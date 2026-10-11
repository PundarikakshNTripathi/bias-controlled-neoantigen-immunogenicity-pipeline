"""Module 1: Data Harmonisation Engine (FR-01 to FR-06).

This pass ingests ONE source: the IEDB T cell assay export (`tcell_full_v3`).
TESLA-derived and HiTIDE corpora are declared in the contract and are not
ingested yet; adding them is a matter of producing another IEDB-shaped frame
and letting the same dedup/union path merge it, not of changing this logic.

Harmonisation is a substantive step, not a clerical one: errors introduced
here propagate silently into every downstream measurement.

Required behaviour
------------------
FR-01  Ingest IEDB, TESLA-derived and HiTIDE corpora into one canonical frame.
       (IEDB T cell export only in this pass.)
FR-02  Reconcile allele nomenclature; assign 34-residue contact pseudo-sequences.
FR-03  Filter peptides to 8-14 residues; deduplicate across overlapping sources.
FR-04  Retain source dataset, assay type and negative-sample provenance per record.
FR-05  Emit a content-hashed, frozen canonical dataset.
FR-06  Reproduce published dataset statistics where a reference exists; document
       any discrepancy rather than quietly matching.

The unit of analysis
--------------------
One canonical row is one **(peptide, allele) pair**, not one assay row. Rows
that share a pair are repeated measurements of one observation (analysis plan,
section 2), so every pair carries `n_assay_records`,
`n_positive_assay_records` and `positive_fraction` alongside the union of its
records' provenance.

A pair whose assay records disagree has to be labelled by a rule, and the rule
changes the corpus. The three pre-registered rules are selected by
`label_rule`, which is part of the hashed criterion, so each rule is a
separate artefact:

* `majority` (the pre-registered PRIMARY rule): positive when more than half of
  the pair's assay records are positive, negative when more than half are
  negative, dropped under `pair_label_tie` otherwise.
* `strict` (sensitivity analysis (a)): every pair whose assay records disagree
  is dropped under `duplicate_label_conflict`.
* `any_positive` (sensitivity analysis (c)): positive when any assay record is
  positive. No pair is dropped.

Choosing between them on the basis of a downstream result is exactly the
failure this project is about; the primary rule is fixed in
`docs/analysis_plan.md` before any model exists.

Hard rules
----------
* A record without provenance is dropped, never defaulted. Positive-to-negative
  ratios are artefacts of experimental design, so an analysis that pools across
  sources without provenance is uninterpretable.
* An allele lacking sequence information is either removed or mapped to a
  same-two-digit-typing representative, and the choice is recorded PER RECORD
  in `allele_mapping_note`.
* Record counts reported anywhere must come from the frame, never from a
  constant in the source. This module contains no corpus statistic, and
  `tests/test_harmonise.py` asserts that by parsing it.
* Degenerate input (an empty frame, or a frame whose every record is dropped)
  raises. It never returns an empty corpus.

Filter order
------------
The order below is fixed, because it decides which bucket a record that fails
several checks is counted in, and therefore what every drop count MEANS.

    missing_provenance -> epitope object type -> modified residues ->
    MHC class -> allele reconciliation -> pseudo-sequence availability ->
    peptide alphabet/length -> label mappability -> pair collapse ->
    pair label rule

Two choices in that order are deliberate:

* Provenance is checked FIRST. A record that cannot say where it came from is
  not a record of anything, so it is never counted under a scientific reason.
* Scope (class I, human HLA, resolvable allele) is checked BEFORE peptide
  length. A first implementation checked length first, and the class II
  peptides in the T cell export - which are mostly longer than 14 residues -
  landed in `peptide_too_long`, which made the reported length-exclusion
  count a mixture of "too long for class I" and "not class I at all". With
  this order, `peptide_too_short` / `peptide_too_long` count only in-scope
  human class I records, which is the quantity worth reporting.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd

from nbx.data.alleles import (
    ALLELE_DROP_REASONS,
    AlleleResolution,
    pseudo_sequence_key,
    reconcile_allele,
)
from nbx.evaluate.gate import UNMEASURED, render_metric

__all__ = [
    "PEPTIDE_MIN_LENGTH",
    "PEPTIDE_MAX_LENGTH",
    "PSEUDO_SEQUENCE_LENGTH",
    "QUALITATIVE_LABEL_MAP",
    "RAW_COLUMNS",
    "CANONICAL_COLUMNS",
    "DROP_REASONS",
    "SOURCE_IEDB_TCELL",
    "LABEL_RULE_MAJORITY",
    "LABEL_RULE_STRICT",
    "LABEL_RULE_ANY_POSITIVE",
    "LABEL_RULES",
    "DEFAULT_LABEL_RULE",
    "LABEL_RULE_DESCRIPTIONS",
    "Provenance",
    "PseudoSequenceTable",
    "PseudoSequenceRegistry",
    "HarmonisedCorpus",
    "ReferenceStatistic",
    "REFERENCE_STATISTICS",
    "HarmonisationError",
    "DegenerateCorpusError",
    "ArtefactHashConflictError",
    "ArtefactIntegrityError",
    "UnsupportedSourceError",
    "UnknownLabelRuleError",
    "compare_to_references",
    "read_iedb_tcell",
    "harmonise_frame",
    "write_corpus",
    "harmonise",
    "default_artefact_name",
]

PEPTIDE_MIN_LENGTH = 8
PEPTIDE_MAX_LENGTH = 14
PSEUDO_SEQUENCE_LENGTH = 34   # NetMHCpan convention: polymorphic residues
                              # within 4.0 A of the peptide

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

SOURCE_IEDB_TCELL = "IEDB:tcell_full_v3"
"""Identity of the only source ingested in this pass."""

NEGATIVE_ORIGIN_IEDB_ASSAY = "iedb_assay_negative"
"""IEDB negatives are assay-reported non-responses, not sampled decoys.

This distinction is the reason `negative_sample_origin` exists: a decoy-derived
negative set and an assay-reported negative set do not have comparable
positive-to-negative ratios, and pooling them without the label makes the
ratio uninterpretable.
"""

NEGATIVE_ORIGIN_NOT_APPLICABLE = "not_applicable:positive_record"
"""Explicit marker on positive records. Never an empty cell, never a default."""

PMID_UNREPORTED = "pmid_unreported"
"""Explicit marker for a record IEDB ships without a PubMed identifier.

The PMID is auxiliary provenance: the mandatory triple is source dataset,
assay type and negative-sample origin. Records with no PMID are retained and
counted separately in the report (`records_without_reference_pmid`) rather
than dropped or silently blanked.
"""

# --------------------------------------------------------------------------
# The qualitative-measurement mapping, stated once, in the open.
# IEDB grades a positive T cell response by magnitude. This project's label is
# binary recognition, so every graded positive maps to 1. The grades are NOT
# averaged, thresholded or weighted anywhere; if a later pass wants to use
# them, it must read them from the raw export, not from this frame.
# --------------------------------------------------------------------------
QUALITATIVE_LABEL_MAP: dict[str, int] = {
    "Positive": 1,
    "Positive-High": 1,
    "Positive-Intermediate": 1,
    "Positive-Low": 1,
    "Negative": 0,
}

# Raw IEDB columns this module reads, addressed as 'Group|Field' because the
# export carries two header rows (group row, then field row).
RAW_COLUMNS: tuple[str, ...] = (
    "Assay ID|IEDB IRI",
    "Reference|PMID",
    "Epitope|Object Type",
    "Epitope|Name",
    "Epitope|Modified residues",
    "Assay|Method",
    "Assay|Response measured",
    "Assay|Qualitative Measurement",
    "MHC Restriction|Name",
    "MHC Restriction|Class",
)

# --------------------------------------------------------------------------
# The pre-registered label rules (docs/analysis_plan.md, section 2).
# `majority` is the primary; the other two are the declared sensitivity
# analyses. The rule is part of the hashed criterion, so switching rules
# produces a different artefact rather than quietly replacing one.
# --------------------------------------------------------------------------
LABEL_RULE_MAJORITY = "majority"
LABEL_RULE_STRICT = "strict"
LABEL_RULE_ANY_POSITIVE = "any_positive"

LABEL_RULES: tuple[str, ...] = (
    LABEL_RULE_MAJORITY,
    LABEL_RULE_STRICT,
    LABEL_RULE_ANY_POSITIVE,
)

DEFAULT_LABEL_RULE = LABEL_RULE_MAJORITY
"""The primary rule, fixed in the pre-registration before any model existed."""

LABEL_RULE_DESCRIPTIONS: dict[str, str] = {
    LABEL_RULE_MAJORITY: (
        "majority vote over the pair's assay records; positive when more than "
        "half are positive, negative when more than half are negative, pair "
        "dropped as a tie otherwise (PRIMARY, pre-registered)"
    ),
    LABEL_RULE_STRICT: (
        "every pair whose assay records disagree is dropped (sensitivity "
        "analysis (a))"
    ),
    LABEL_RULE_ANY_POSITIVE: (
        "a pair is positive when any of its assay records is positive "
        "(sensitivity analysis (c))"
    ),
}

PAIR_DROP_REASON_BY_RULE: dict[str, str] = {
    LABEL_RULE_MAJORITY: "pair_label_tie",
    LABEL_RULE_STRICT: "duplicate_label_conflict",
}
"""The reason each rule may give for dropping a pair. `any_positive` drops none."""

CANONICAL_COLUMNS: tuple[str, ...] = (
    "peptide",
    "peptide_length",
    "allele",
    "allele_pseudoseq",
    "label",
    "n_assay_records",
    "n_positive_assay_records",
    "positive_fraction",
    "source_dataset",
    "assay_type",
    "negative_sample_origin",
    "allele_mapping_note",
    "assay_response_measured",
    "reference_pmids",
    "contributing_assay_ids",
    "contributing_sources",
    "n_contributing_sources",
    "n_contributing_records",
)

EPITOPE_OBJECT_TYPE_KEPT = "Linear peptide"
MHC_CLASS_KEPT = "I"

DROP_REASONS: tuple[str, ...] = (
    "missing_provenance",
    "epitope_not_linear_peptide",
    "epitope_has_modified_residues",
    "peptide_missing",
    "peptide_non_standard_residues",
    "peptide_too_short",
    "peptide_too_long",
    "label_not_mappable",
    "mhc_class_missing",
    "mhc_class_not_i",
    *ALLELE_DROP_REASONS,
    "allele_two_digit_mapping_ambiguous",
    "allele_not_in_pseudo_sequence_table",
    "allele_pseudo_sequence_unavailable",
    "duplicate_label_conflict",
    "pair_label_tie",
)

FILTER_ORDER: tuple[str, ...] = (
    "provenance",
    "epitope_object_type",
    "modified_residues",
    "mhc_class",
    "allele_reconciliation",
    "pseudo_sequence_availability",
    "peptide_alphabet_and_length",
    "label_mappability",
    "pair_collapse",
    "pair_label_rule",
)

MANDATORY_PROVENANCE_FIELDS: tuple[str, ...] = ("assay_iri", "assay_type")
"""Raw fields whose absence means the record has lost its provenance.

`source_dataset` is supplied by the caller and `negative_sample_origin` is
derived from the label, so the two raw fields that can actually go missing are
the assay identifier and the assay method.
"""

_DROPPED_PAIR_EXAMPLE_LIMIT = 20
_ABSENT_ALLELE_LIST_LIMIT = 25


class HarmonisationError(RuntimeError):
    """Base class for every refusal this module makes."""


class DegenerateCorpusError(HarmonisationError):
    """The input produced no usable records. Raised, never returned."""


class ArtefactHashConflictError(HarmonisationError):
    """An artefact of the same name exists with different content."""


class ArtefactIntegrityError(HarmonisationError):
    """A written artefact no longer hashes to the value its report records.

    Raised rather than repaired. A report is a description of a specific
    parquet; if the parquet has changed underneath it, refreshing the report
    would attach a true-looking description to an unknown file.
    """


class UnsupportedSourceError(HarmonisationError):
    """A source was requested that this pass does not ingest."""


class UnknownLabelRuleError(HarmonisationError):
    """A label rule was requested that is not one of the pre-registered three."""


def _pair_label(rule: str, n_positive: int, n_records: int) -> tuple[int | None, str]:
    """Apply one label rule to one pair's assay records.

    Returns `(label, drop_reason)`; the label is None exactly when the rule
    refuses to label the pair, and `drop_reason` is then the named reason.
    """
    n_negative = n_records - n_positive
    if rule == LABEL_RULE_MAJORITY:
        if n_positive > n_negative:
            return 1, ""
        if n_negative > n_positive:
            return 0, ""
        return None, PAIR_DROP_REASON_BY_RULE[LABEL_RULE_MAJORITY]
    if rule == LABEL_RULE_STRICT:
        if n_negative == 0:
            return 1, ""
        if n_positive == 0:
            return 0, ""
        return None, PAIR_DROP_REASON_BY_RULE[LABEL_RULE_STRICT]
    if rule == LABEL_RULE_ANY_POSITIVE:
        return (1 if n_positive > 0 else 0), ""
    raise UnknownLabelRuleError(
        f"unknown label rule {rule!r}; the pre-registered rules are "
        f"{list(LABEL_RULES)}"
    )


@dataclass(frozen=True)
class Provenance:
    """Where a record came from. Mandatory on every record."""

    source_dataset: str
    assay_type: str
    negative_sample_origin: str
    allele_mapping_note: str = ""

    def complete(self) -> bool:
        """True when every mandatory field carries real content."""
        return all(
            bool(str(value).strip())
            for value in (
                self.source_dataset,
                self.assay_type,
                self.negative_sample_origin,
            )
        )


# --------------------------------------------------------------------------
# Pseudo-sequences (FR-02)
# --------------------------------------------------------------------------
class PseudoSequenceTable:
    """One allele -> 34-residue contact pseudo-sequence table.

    Deliberately a small interface rather than a dict, so that a second
    source (a NetMHCpan table, an IMGT-derived table) can be added behind the
    same calls and so that the table that resolved each allele can be named
    on the record.
    """

    def __init__(self, entries: Mapping[str, str], source_name: str) -> None:
        if not source_name.strip():
            raise HarmonisationError("a pseudo-sequence table must name its source")
        self.source_name = source_name
        self._by_key: dict[str, str] = {}
        self._by_group: dict[tuple[str, str], dict[str, str]] = defaultdict(dict)
        for raw_name, sequence in entries.items():
            seq = sequence.strip().upper()
            if len(seq) != PSEUDO_SEQUENCE_LENGTH:
                continue
            key = raw_name.strip().upper().replace("*", "")
            self._by_key[key] = seq
            resolution = reconcile_allele(raw_name)
            if resolution.resolved and resolution.protein is not None:
                group = (str(resolution.locus), str(resolution.group))
                self._by_group[group][str(resolution.canonical)] = seq

    def __len__(self) -> int:
        return len(self._by_key)

    def lookup(self, locus: str, group: str, protein: str) -> str | None:
        """Exact two-field lookup. None when the table does not carry it."""
        return self._by_key.get(pseudo_sequence_key(locus, group, protein).upper())

    def group_members(self, locus: str, group: str) -> dict[str, str]:
        """Every two-field allele in one two-digit group, with its sequence."""
        return dict(self._by_group.get((locus, group), {}))


class PseudoSequenceRegistry:
    """An ordered priority chain of pseudo-sequence tables.

    The first table that can resolve an allele wins, and its `source_name` is
    recorded on the record. Later tables are fallbacks, not overrides.
    """

    def __init__(self, tables: Sequence[PseudoSequenceTable]) -> None:
        if not tables:
            raise HarmonisationError(
                "at least one pseudo-sequence table is required; alleles cannot "
                "be assigned contact residues without one"
            )
        self.tables = tuple(tables)

    @property
    def source_names(self) -> list[str]:
        return [table.source_name for table in self.tables]

    def resolve(self, resolution: AlleleResolution) -> tuple[str | None, str, str | None]:
        """Resolve one reconciled allele to a pseudo-sequence.

        Returns `(sequence, note, drop_reason)`. Exactly one of `sequence` and
        `drop_reason` is not None. `note` always describes what was done, and
        is carried on the record.

        The two failure reasons are kept apart deliberately, because an
        earlier version conflated them and so reported a table coverage gap as
        an ambiguity in the input:

        * `allele_not_in_pseudo_sequence_table` - the input WAS resolved to a
          two-field allele; no table carries it, and its two-digit group does
          not admit an unambiguous representative. The limit is the table.
        * `allele_two_digit_mapping_ambiguous` - the input carried only
          two-digit resolution, and the group it names covers molecules with
          different contact residues. The limit is the input.
        """
        locus, group = str(resolution.locus), str(resolution.group)
        protein = resolution.protein
        fully_resolved = protein is not None

        if protein is not None:
            for table in self.tables:
                found = table.lookup(locus, group, protein)
                if found is not None:
                    return found, f"pseudoseq=direct_lookup:{table.source_name}", None

        # No exact entry: may we use a same-two-digit representative? Only if
        # every two-field allele in the group that the table knows about has
        # the same contact residues. Otherwise the mapping would be a guess.
        for table in self.tables:
            members = table.group_members(locus, group)
            if not members:
                continue
            distinct = sorted(set(members.values()))
            if len(distinct) == 1:
                representatives = ",".join(sorted(members))
                note = (
                    f"pseudoseq=two_digit_representative:{table.source_name}"
                    f"; mapped={resolution.canonical}->group {locus}*{group}"
                    f"; representatives={representatives}"
                )
                return distinct[0], note, None
            if fully_resolved:
                return (
                    None,
                    (
                        f"pseudoseq=unresolved; {resolution.canonical} is absent "
                        f"from {table.source_name} and its group {locus}*{group} "
                        f"carries {len(distinct)} distinct pseudo-sequences there, "
                        "so no representative may be chosen"
                    ),
                    "allele_not_in_pseudo_sequence_table",
                )
            return (
                None,
                (
                    f"pseudoseq=unresolved; group {locus}*{group} carries "
                    f"{len(distinct)} distinct pseudo-sequences in "
                    f"{table.source_name}"
                ),
                "allele_two_digit_mapping_ambiguous",
            )

        unseen = (
            f"pseudoseq=unresolved; no entry for {resolution.canonical} in "
            + ",".join(self.source_names)
        )
        if fully_resolved:
            return None, unseen, "allele_not_in_pseudo_sequence_table"
        return None, unseen, "allele_pseudo_sequence_unavailable"


def read_pseudo_sequence_dat(path: Path, source_name: str) -> PseudoSequenceTable:
    """Read a whitespace-separated `NAME SEQUENCE` pseudo-sequence file.

    The file is read from `data/raw/`; it is never vendored into this
    repository, because its origin and redistribution terms are unverified
    (see docs/progress.md, M1 Step 0).
    """
    entries: dict[str, str] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        entries[parts[0]] = parts[1]
    if not entries:
        raise HarmonisationError(f"no usable pseudo-sequence entries in {path}")
    return PseudoSequenceTable(entries=entries, source_name=source_name)


# --------------------------------------------------------------------------
# FR-06 reference comparison
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class ReferenceStatistic:
    """A published dataset statistic to compare the built frame against.

    `locator` is mandatory and must point at something a reader can check
    (publication plus section/table, or an artefact path plus hash). A
    reference without a locator is not a reference.
    """

    name: str
    value: float | None
    locator: str

    def __post_init__(self) -> None:
        if not self.locator.strip():
            raise HarmonisationError(
                f"reference statistic {self.name!r} has no locator; an "
                "uncitable number may not be used as a reference"
            )


NO_REFERENCE_NOTE = (
    "No published reference statistic has been located for this subset of the "
    "IEDB T cell export (linear class I human peptides, 8-14 residues, with "
    "the allele reconciliation and deduplication rules applied here). The "
    "comparison therefore reports every reference as unmeasured. This is a "
    "gap in the literature search, not a match: nothing here has been checked "
    "against an external number."
)

REFERENCE_STATISTICS: tuple[ReferenceStatistic, ...] = ()
"""Empty on purpose. See NO_REFERENCE_NOTE.

Populate this tuple only with a value that has a locator. The comparison
machinery is implemented and tested, so adding a reference requires no code
change.
"""


def _render_count(value: float | None) -> str:
    if value is None:
        return UNMEASURED
    numeric = float(value)
    return render_metric(numeric, 0) if numeric.is_integer() else render_metric(numeric)


def compare_to_references(
    observed: Mapping[str, float | None],
    references: Sequence[ReferenceStatistic] = REFERENCE_STATISTICS,
    note: str = NO_REFERENCE_NOTE,
) -> dict:
    """Compare observed statistics against published ones (FR-06).

    Where no reference has been located the reference, the discrepancy and the
    locator all render as an em dash. A difference is reported as a
    discrepancy; it is never narrowed, rounded or explained away here.
    """
    by_name = {ref.name: ref for ref in references}
    located = sum(ref.value is not None for ref in references)

    statistics = []
    for name in sorted(observed):
        value = observed[name]
        ref = by_name.get(name)
        if ref is None or ref.value is None:
            statistics.append(
                {
                    "name": name,
                    "observed": _render_count(value),
                    "reference": UNMEASURED,
                    "discrepancy": UNMEASURED,
                    "locator": ref.locator if ref is not None else UNMEASURED,
                    "status": "no_reference_located",
                }
            )
            continue
        if value is None:
            statistics.append(
                {
                    "name": name,
                    "observed": UNMEASURED,
                    "reference": _render_count(ref.value),
                    "discrepancy": UNMEASURED,
                    "locator": ref.locator,
                    "status": "not_measured",
                }
            )
            continue
        delta = float(value) - float(ref.value)
        statistics.append(
            {
                "name": name,
                "observed": _render_count(value),
                "reference": _render_count(ref.value),
                "discrepancy": _render_count(delta),
                "locator": ref.locator,
                "status": "reproduced" if delta == 0 else "discrepancy",
            }
        )

    return {
        "references_located": located,
        "note": note if located == 0 else "",
        "statistics": statistics,
    }


# --------------------------------------------------------------------------
# Reading the raw export
# --------------------------------------------------------------------------
def _iedb_header_names(reader: Iterable[list[str]]) -> list[str]:
    rows = iter(reader)
    group_row = next(rows)
    field_row = next(rows)
    if len(group_row) != len(field_row):
        raise HarmonisationError(
            "the two IEDB header rows have different widths; the export is not "
            "the expected shape"
        )
    names: list[str] = []
    last_group = ""
    for group, field_name in zip(group_row, field_row):
        group = group.strip() or last_group
        last_group = group
        names.append(f"{group}|{field_name.strip()}")
    return names


def read_iedb_tcell(path: Path, columns: Sequence[str] = RAW_COLUMNS) -> pd.DataFrame:
    """Read the needed columns of an IEDB T cell export.

    The export carries TWO header rows (group, then field), which pandas
    cannot combine with `usecols`. So the header is read with `csv` first, the
    needed columns are located by position, and the body is then read
    positionally.
    """
    path = Path(path)

    def _open():
        if path.suffix.lower() == ".zip":
            archive = zipfile.ZipFile(path)
            members = [n for n in archive.namelist() if n.lower().endswith(".csv")]
            if len(members) != 1:
                raise HarmonisationError(
                    f"expected exactly one CSV inside {path.name}, found {members}"
                )
            return io.TextIOWrapper(
                archive.open(members[0]), encoding="utf-8-sig", newline=""
            )
        return path.open("r", encoding="utf-8-sig", newline="")

    with _open() as handle:
        names = _iedb_header_names(csv.reader(handle))

    missing = [c for c in columns if c not in names]
    if missing:
        raise HarmonisationError(
            f"{path.name} does not carry the expected column(s): {missing}"
        )
    indices = [names.index(c) for c in columns]
    ordered = sorted(indices)

    with _open() as handle:
        frame = pd.read_csv(
            handle,
            header=None,
            skiprows=2,
            usecols=ordered,
            dtype=str,
            keep_default_na=False,
            na_values=[],
            low_memory=False,
        )
    frame.columns = [names[i] for i in ordered]
    return frame[list(columns)]


# --------------------------------------------------------------------------
# The canonical frame
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class HarmonisedCorpus:
    """A built canonical frame, its report and its content hash."""

    frame: pd.DataFrame
    report: dict
    content_hash: str = field(default="", compare=False)

    @property
    def short_hash(self) -> str:
        return f"sha256:{self.content_hash[:6]}..{self.content_hash[-3:]}"


def _clean(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def harmonisation_criterion(
    source_dataset: str,
    pseudo_sources: Sequence[str],
    label_rule: str = DEFAULT_LABEL_RULE,
) -> dict:
    """The rule set that defines a harmonised corpus, as a hashable dict.

    Analogous to `nbx.splits.freeze.SplitCriterion`: it is part of the
    artefact's identity, because two different rule sets can produce the same
    rows while meaning different things. In particular the filter order
    decides what every drop count means, so changing it changes the artefact
    even when the retained frame is byte-identical.

    `label_rule` is in here for the same reason: it decides what a pair whose
    assays disagree is labelled, so each rule must be its own artefact with
    its own hash.
    """
    if label_rule not in LABEL_RULES:
        raise UnknownLabelRuleError(
            f"unknown label rule {label_rule!r}; the pre-registered rules are "
            f"{list(LABEL_RULES)}"
        )
    return {
        "filter_order": list(FILTER_ORDER),
        "unit_of_analysis": "peptide_allele_pair",
        "label_rule": label_rule,
        "label_rule_description": LABEL_RULE_DESCRIPTIONS[label_rule],
        "label_map": {str(k): int(v) for k, v in QUALITATIVE_LABEL_MAP.items()},
        "peptide_length_window": [PEPTIDE_MIN_LENGTH, PEPTIDE_MAX_LENGTH],
        "pseudo_sequence_length": PSEUDO_SEQUENCE_LENGTH,
        "pseudo_sequence_sources": list(pseudo_sources),
        "epitope_object_type_kept": EPITOPE_OBJECT_TYPE_KEPT,
        "mhc_class_kept": MHC_CLASS_KEPT,
        "mandatory_provenance_fields": list(MANDATORY_PROVENANCE_FIELDS),
        "pair_drop_reason": PAIR_DROP_REASON_BY_RULE.get(label_rule, ""),
        "source_dataset": source_dataset,
        "schema": list(CANONICAL_COLUMNS),
    }


def content_hash(frame: pd.DataFrame, criterion: Mapping[str, object]) -> str:
    """sha256 over the canonical JSON of the criterion plus the frame content.

    Same pattern as `nbx.splits.freeze.FrozenSplit.compute_hash`: canonical
    JSON with sorted keys and no insignificant whitespace, hashed. Rows are
    hashed in the frame's (deterministically sorted) order, one at a time, so
    this works on a large frame without building one giant string.
    """
    digest = hashlib.sha256()
    digest.update(
        json.dumps(
            {
                "criterion": dict(criterion),
                "schema": list(frame.columns),
                "n_records": int(len(frame)),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )
    for record in frame.to_dict("records"):
        digest.update(
            json.dumps(
                {k: (int(v) if isinstance(v, bool) else v) for k, v in record.items()},
                sort_keys=True,
                separators=(",", ":"),
                default=str,
            ).encode("utf-8")
        )
        digest.update(b"\n")
    return digest.hexdigest()


def harmonise_frame(
    raw: pd.DataFrame,
    pseudo: PseudoSequenceRegistry,
    source_dataset: str = SOURCE_IEDB_TCELL,
    input_files: Mapping[str, dict] | None = None,
    label_rule: str = DEFAULT_LABEL_RULE,
) -> HarmonisedCorpus:
    """Build the canonical frame from one IEDB-shaped raw frame.

    One output row is one (peptide, allele) pair, labelled by `label_rule`.
    Every count in the returned report is computed from the frame that was
    just built or from the counters incremented while building it.
    """
    if label_rule not in LABEL_RULES:
        raise UnknownLabelRuleError(
            f"unknown label rule {label_rule!r}; the pre-registered rules are "
            f"{list(LABEL_RULES)}"
        )
    missing = [c for c in RAW_COLUMNS if c not in raw.columns]
    if missing:
        raise HarmonisationError(f"raw frame is missing column(s): {missing}")
    if len(raw) == 0:
        raise DegenerateCorpusError(
            "the raw frame carries no rows; there is nothing to harmonise"
        )

    dropped: Counter[str] = Counter()
    missing_provenance_fields: Counter[str] = Counter()
    qualitative_seen: Counter[str] = Counter()
    unresolved_allele_strings: Counter[str] = Counter()
    absent_from_pseudo_table: Counter[str] = Counter()
    two_digit_mappings = 0
    kept: list[dict] = []

    for row in raw.to_dict("records"):
        assay_iri = _clean(row["Assay ID|IEDB IRI"])
        assay_type = _clean(row["Assay|Method"])

        # 1. Provenance, first. A record that cannot say where it came from is
        #    dropped and counted. Nothing is defaulted in.
        provenance_gaps = [
            name
            for name, value in (("assay_iri", assay_iri), ("assay_type", assay_type))
            if not value
        ]
        if provenance_gaps:
            dropped["missing_provenance"] += 1
            for name in provenance_gaps:
                missing_provenance_fields[name] += 1
            continue

        # 2. Epitope must be an unmodified linear peptide.
        if _clean(row["Epitope|Object Type"]) != EPITOPE_OBJECT_TYPE_KEPT:
            dropped["epitope_not_linear_peptide"] += 1
            continue
        if _clean(row["Epitope|Modified residues"]):
            dropped["epitope_has_modified_residues"] += 1
            continue

        # 3. MHC class.
        mhc_class = _clean(row["MHC Restriction|Class"])
        if not mhc_class:
            dropped["mhc_class_missing"] += 1
            continue
        if mhc_class.upper() != MHC_CLASS_KEPT:
            dropped["mhc_class_not_i"] += 1
            continue

        # 4. Allele nomenclature.
        raw_allele = _clean(row["MHC Restriction|Name"])
        resolution = reconcile_allele(raw_allele)
        if not resolution.resolved:
            dropped[str(resolution.drop_reason)] += 1
            unresolved_allele_strings[raw_allele or UNMEASURED] += 1
            continue

        # 5. Pseudo-sequence.
        sequence, pseudo_note, pseudo_drop = pseudo.resolve(resolution)
        if sequence is None:
            dropped[str(pseudo_drop)] += 1
            unresolved_allele_strings[raw_allele] += 1
            if pseudo_drop == "allele_not_in_pseudo_sequence_table":
                absent_from_pseudo_table[str(resolution.canonical)] += 1
            continue

        # 6. Peptide alphabet and length. Checked after scope, so that these
        #    counts describe in-scope human class I records only.
        peptide = _clean(row["Epitope|Name"]).upper()
        if not peptide:
            dropped["peptide_missing"] += 1
            continue
        if not set(peptide) <= STANDARD_AMINO_ACIDS:
            dropped["peptide_non_standard_residues"] += 1
            continue
        if len(peptide) < PEPTIDE_MIN_LENGTH:
            dropped["peptide_too_short"] += 1
            continue
        if len(peptide) > PEPTIDE_MAX_LENGTH:
            dropped["peptide_too_long"] += 1
            continue

        # 7. Label.
        qualitative = _clean(row["Assay|Qualitative Measurement"])
        qualitative_seen[qualitative or UNMEASURED] += 1
        if qualitative not in QUALITATIVE_LABEL_MAP:
            dropped["label_not_mappable"] += 1
            continue
        label = QUALITATIVE_LABEL_MAP[qualitative]

        if "two_digit_representative" in pseudo_note:
            two_digit_mappings += 1

        note_parts = [f"raw={resolution.raw}", f"canonical={resolution.canonical}"]
        note_parts.extend(resolution.notes)
        note_parts.append(pseudo_note)
        provenance = Provenance(
            source_dataset=source_dataset,
            assay_type=assay_type,
            negative_sample_origin=(
                NEGATIVE_ORIGIN_IEDB_ASSAY
                if label == 0
                else NEGATIVE_ORIGIN_NOT_APPLICABLE
            ),
            allele_mapping_note="; ".join(note_parts),
        )
        if not provenance.complete():
            # Defensive: should be unreachable, but a silent default here would
            # be exactly the failure this module exists to prevent.
            dropped["missing_provenance"] += 1
            continue

        pmid = _clean(row["Reference|PMID"]) or PMID_UNREPORTED
        kept.append(
            {
                "peptide": peptide,
                "peptide_length": len(peptide),
                "allele": str(resolution.canonical),
                "allele_pseudoseq": sequence,
                "label": int(label),
                "assay_response_measured": _clean(row["Assay|Response measured"])
                or UNMEASURED,
                "reference_pmid": pmid,
                "assay_id": assay_iri.rstrip("/").rsplit("/", maxsplit=1)[-1],
                **asdict(provenance),
            }
        )

    # 8. Collapse to one row per (peptide, allele) pair, keeping the union of
    #    provenance, then label the pair with the selected rule. Rows sharing a
    #    pair are repeated measurements of one observation, so the pair - not
    #    the assay row - is the unit of analysis from here on.
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for record in kept:
        groups[(record["peptide"], record["allele"])].append(record)

    record_level_dropped = int(sum(dropped.values()))

    merged: list[dict] = []
    dropped_pairs: list[tuple[str, str]] = []
    dropped_pair_records = 0
    dropped_pairs_by_reason: Counter[str] = Counter()
    dropped_pair_records_by_label: Counter[int] = Counter()
    pattern: Counter[str] = Counter()
    pre_dedup_by_label: Counter[int] = Counter(record["label"] for record in kept)
    for pair in sorted(groups):
        members = groups[pair]
        n_assay_records = len(members)
        n_positive = int(sum(int(m["label"]) for m in members))
        n_negative = n_assay_records - n_positive

        # The disagreement pattern is a property of the data, counted over
        # every pair that reached the rule - including the ones it drops.
        if n_negative == 0:
            pattern["unanimous_positive"] += 1
        elif n_positive == 0:
            pattern["unanimous_negative"] += 1
        else:
            pattern["mixed"] += 1
        if n_positive == n_negative:
            pattern["tied"] += 1

        pair_label, pair_drop_reason = _pair_label(
            label_rule, n_positive, n_assay_records
        )
        if pair_label is None:
            dropped_pairs.append(pair)
            dropped_pair_records += n_assay_records
            dropped_pairs_by_reason[pair_drop_reason] += 1
            dropped_pair_records_by_label.update(m["label"] for m in members)
            dropped[pair_drop_reason] += n_assay_records
            continue

        first = members[0]
        sources = sorted({m["source_dataset"] for m in members})
        merged.append(
            {
                "peptide": first["peptide"],
                "peptide_length": first["peptide_length"],
                "allele": first["allele"],
                "allele_pseudoseq": first["allele_pseudoseq"],
                "label": int(pair_label),
                "n_assay_records": n_assay_records,
                "n_positive_assay_records": n_positive,
                "positive_fraction": n_positive / n_assay_records,
                "source_dataset": ";".join(sources),
                "assay_type": ";".join(sorted({m["assay_type"] for m in members})),
                "negative_sample_origin": ";".join(
                    sorted({m["negative_sample_origin"] for m in members})
                ),
                "allele_mapping_note": ";".join(
                    sorted({m["allele_mapping_note"] for m in members})
                ),
                "assay_response_measured": ";".join(
                    sorted({m["assay_response_measured"] for m in members})
                ),
                "reference_pmids": ";".join(
                    sorted({m["reference_pmid"] for m in members})
                ),
                "contributing_assay_ids": ";".join(
                    sorted({m["assay_id"] for m in members})
                ),
                "contributing_sources": ";".join(sources),
                "n_contributing_sources": len(sources),
                "n_contributing_records": len(members),
            }
        )

    records_in_labelled_pairs = len(kept) - dropped_pair_records
    duplicates_merged = records_in_labelled_pairs - len(merged)

    if not merged:
        raise DegenerateCorpusError(
            "every record was dropped; the harmonised corpus would be empty. "
            f"drop reasons: {dict(sorted(dropped.items()))}"
        )

    frame = pd.DataFrame(merged, columns=list(CANONICAL_COLUMNS))
    frame = frame.sort_values(["peptide", "allele"]).reset_index(drop=True)
    for column in (
        "label",
        "peptide_length",
        "n_assay_records",
        "n_positive_assay_records",
        "n_contributing_sources",
        "n_contributing_records",
    ):
        frame[column] = frame[column].astype("int64")
    frame["positive_fraction"] = frame["positive_fraction"].astype("float64")

    criterion = harmonisation_criterion(
        source_dataset, pseudo.source_names, label_rule=label_rule
    )
    digest = content_hash(frame, criterion)

    # --- every number below is read back off the frame just built ----------
    label_counts = {str(k): int(v) for k, v in frame["label"].value_counts().items()}
    length_distribution = {
        str(k): int(v) for k, v in frame["peptide_length"].value_counts().items()
    }
    per_source: dict[str, dict] = {}
    for source, subset in frame.groupby("source_dataset", sort=True):
        per_source[str(source)] = {
            "record_count": int(len(subset)),
            "allele_count": int(subset["allele"].nunique()),
            "n_positive": int((subset["label"] == 1).sum()),
            "n_negative": int((subset["label"] == 0).sum()),
        }
    allele_record_counts = {
        str(k): int(v) for k, v in frame["allele"].value_counts().items()
    }
    observed = {
        "record_count": float(len(frame)),
        "allele_count": float(frame["allele"].nunique()),
        "n_positive": float((frame["label"] == 1).sum()),
        "n_negative": float((frame["label"] == 0).sum()),
    }

    absent_ranked = sorted(
        absent_from_pseudo_table.items(), key=lambda kv: (-kv[1], kv[0])
    )
    pairs_total = int(len(merged) + len(dropped_pairs))
    reconciliation = {
        "note": (
            "one canonical row is one (peptide, allele) pair. Record-level "
            "drops, pair-level drops, records collapsed into a surviving pair "
            "and the surviving pairs themselves account for every raw row."
        ),
        "raw_row_count": int(len(raw)),
        "records_dropped_record_level": int(record_level_dropped),
        "records_dropped_pair_level": int(dropped_pair_records),
        "records_in_labelled_pairs": int(records_in_labelled_pairs),
        "records_collapsed_into_pairs": int(duplicates_merged),
        "pairs_labelled": int(len(merged)),
        "pairs_dropped": int(len(dropped_pairs)),
        "pairs_total": pairs_total,
        "record_level_balances": bool(
            record_level_dropped
            + dropped_pair_records
            + duplicates_merged
            + len(merged)
            == len(raw)
        ),
        "pair_level_balances": bool(
            len(merged) + len(dropped_pairs) == pairs_total
            and records_in_labelled_pairs == int(frame["n_assay_records"].sum())
        ),
    }

    report = {
        "module": "M1 nbx.data.harmonise",
        "source_dataset": source_dataset,
        "sources_ingested": [source_dataset],
        "sources_declared_not_ingested": ["TESLA-derived", "HiTIDE/NCI"],
        "input_files": dict(input_files or {}),
        "criterion": criterion,
        "content_hash": digest,
        "schema": list(frame.columns),
        "unit_of_analysis": "peptide_allele_pair",
        "label_rule": label_rule,
        "label_rule_description": LABEL_RULE_DESCRIPTIONS[label_rule],
        "label_rules_available": list(LABEL_RULES),
        "primary_label_rule": DEFAULT_LABEL_RULE,
        "record_count": int(len(frame)),
        "pair_count": int(len(frame)),
        "allele_count": int(frame["allele"].nunique()),
        "peptide_count": int(frame["peptide"].nunique()),
        "label_counts": label_counts,
        "length_distribution": length_distribution,
        "per_source_breakdown": per_source,
        "records_per_allele": allele_record_counts,
        "raw_row_count": int(len(raw)),
        "records_before_deduplication": int(len(kept)),
        "duplicate_records_merged": int(duplicates_merged),
        "reconciliation": reconciliation,
        "pair_pattern_summary": {
            "note": (
                "counted over every (peptide, allele) pair that reached the "
                "label rule, including the pairs the rule then dropped. Tied "
                "pairs are a SUBSET of mixed pairs; unanimous_positive + "
                "unanimous_negative + mixed equals pairs_total."
            ),
            "pairs_total": pairs_total,
            "unanimous_positive": int(pattern["unanimous_positive"]),
            "unanimous_negative": int(pattern["unanimous_negative"]),
            "mixed": int(pattern["mixed"]),
            "tied": int(pattern["tied"]),
        },
        "pairs_dropped_by_label_rule": int(len(dropped_pairs)),
        "pairs_dropped_by_reason": {
            str(k): int(v) for k, v in sorted(dropped_pairs_by_reason.items())
        },
        # The label rule is not label-neutral: dropping a pair removes
        # positives and negatives in whatever ratio the pair happened to have,
        # and labelling a mixed pair reassigns the minority's records. Either
        # way the corpus-level positive-to-negative ratio moves, and that ratio
        # is itself a design artefact, so the shift is reported rather than
        # left for someone to discover.
        "label_counts_before_deduplication": {
            str(k): int(v) for k, v in sorted(pre_dedup_by_label.items())
        },
        "records_in_dropped_pairs_by_label": {
            str(k): int(v) for k, v in sorted(dropped_pair_records_by_label.items())
        },
        "dropped_pair_examples": [
            list(pair) for pair in dropped_pairs[:_DROPPED_PAIR_EXAMPLE_LIMIT]
        ],
        "labelled_pairs_with_mixed_negative_sample_origin": int(
            frame["negative_sample_origin"].str.contains(";", regex=False).sum()
        ),
        "alleles_absent_from_pseudo_sequence_table": {
            "drop_reason": "allele_not_in_pseudo_sequence_table",
            "note": (
                "alleles that WERE resolved to a two-field name but that no "
                "pseudo-sequence table carries, and whose two-digit group does "
                "not admit an unambiguous representative. Their records are "
                "dropped, never approximated. This is a coverage gap in the "
                "table, not an ambiguity in the input."
            ),
            "pseudo_sequence_sources": pseudo.source_names,
            "n_alleles": int(len(absent_from_pseudo_table)),
            "n_records": int(sum(absent_from_pseudo_table.values())),
            "listed_limit": _ABSENT_ALLELE_LIST_LIMIT,
            "n_listed": int(min(len(absent_ranked), _ABSENT_ALLELE_LIST_LIMIT)),
            "alleles": [
                {"allele": str(allele), "records": int(count)}
                for allele, count in absent_ranked[:_ABSENT_ALLELE_LIST_LIMIT]
            ],
        },
        "dropped_by_reason": {str(k): int(v) for k, v in sorted(dropped.items())},
        "dropped_total": int(sum(dropped.values())),
        "dropped_for_missing_provenance": int(dropped["missing_provenance"]),
        "missing_provenance_fields": {
            str(k): int(v) for k, v in sorted(missing_provenance_fields.items())
        },
        "records_without_reference_pmid": int(
            frame["reference_pmids"].eq(PMID_UNREPORTED).sum()
        ),
        "two_digit_representative_mappings": int(two_digit_mappings),
        "peptide_length_window": [PEPTIDE_MIN_LENGTH, PEPTIDE_MAX_LENGTH],
        "peptides_excluded_for_length": int(
            dropped["peptide_too_short"] + dropped["peptide_too_long"]
        ),
        "label_map": {str(k): int(v) for k, v in QUALITATIVE_LABEL_MAP.items()},
        "qualitative_measurements_seen": {
            str(k): int(v) for k, v in sorted(qualitative_seen.items())
        },
        "unresolved_allele_strings": {
            str(k): int(v)
            for k, v in sorted(
                unresolved_allele_strings.items(), key=lambda kv: (-kv[1], kv[0])
            )
        },
        "pseudo_sequence_sources": pseudo.source_names,
        "filter_order": list(FILTER_ORDER),
        "negative_sample_origins": sorted(
            set(frame["negative_sample_origin"].unique())
        ),
        "reference_comparison": compare_to_references(observed),
    }

    # Defensive invariants. A violation here means a downstream measurement
    # would be uninterpretable, so it raises rather than warns.
    for column in ("source_dataset", "assay_type", "negative_sample_origin",
                   "allele_mapping_note"):
        if not frame[column].astype(str).str.strip().all():
            raise HarmonisationError(
                f"column {column!r} has an empty value on at least one retained "
                "record; provenance may not be blank"
            )
    if not frame["allele_pseudoseq"].str.len().eq(PSEUDO_SEQUENCE_LENGTH).all():
        raise HarmonisationError("a retained record has a wrong-length pseudo-sequence")
    if frame.duplicated(["peptide", "allele"]).any():
        raise HarmonisationError("deduplication left a duplicate (peptide, allele) pair")
    if not reconciliation["record_level_balances"]:
        raise HarmonisationError(
            "record reconciliation does not balance: kept pairs, dropped "
            "records and collapsed records do not account for every raw row. "
            f"{reconciliation}"
        )
    if not reconciliation["pair_level_balances"]:
        raise HarmonisationError(
            f"pair reconciliation does not balance: {reconciliation}"
        )
    patterns: dict = report["pair_pattern_summary"]  # type: ignore[assignment]
    if (
        patterns["unanimous_positive"]
        + patterns["unanimous_negative"]
        + patterns["mixed"]
        != patterns["pairs_total"]
    ):
        raise HarmonisationError(
            f"pair disagreement patterns do not partition the pairs: {patterns}"
        )
    if (frame["n_positive_assay_records"] > frame["n_assay_records"]).any():
        raise HarmonisationError("a pair reports more positive assays than assays")

    return HarmonisedCorpus(frame=frame, report=report, content_hash=digest)


# --------------------------------------------------------------------------
# Writing (FR-05)
# --------------------------------------------------------------------------
def write_corpus(
    corpus: HarmonisedCorpus, out_dir: Path, name: str = "iedb_tcell"
) -> dict[str, Path]:
    """Write the dataset and its report, refusing to overwrite a different one.

    Re-writing byte-identical content is a no-op. Writing content with a
    different hash under an existing name raises: a changed criterion is a new
    artefact, not an edit of the old one.
    """
    out_dir = Path(out_dir).resolve()
    if "splits" in out_dir.parts:
        raise HarmonisationError(
            "refusing to write into a splits directory; frozen splits are "
            "immutable (hard rule 4)"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    dataset_path = out_dir / f"{name}.parquet"
    report_path = out_dir / f"{name}.report.json"

    if report_path.exists():
        existing = json.loads(report_path.read_text(encoding="utf-8"))
        if existing.get("content_hash") != corpus.content_hash:
            raise ArtefactHashConflictError(
                f"{report_path.name} exists with content hash "
                f"{existing.get('content_hash')} but the corpus just built has "
                f"{corpus.content_hash}. A different corpus is a different "
                "artefact: write it under a different name."
            )
        # The hash matches, so the criterion and the rows are the same corpus.
        # The report is a regenerable DESCRIPTION of that corpus, so it is
        # refreshed rather than refused: adding a diagnostic to the report is
        # not a new experiment. A changed criterion, by contrast, changes
        # `harmonisation_criterion` and therefore the hash, and lands on the
        # refusal above - which is how the stale-report failure that motivated
        # putting the criterion in the hash is now caught.
        if dataset_path.exists():
            # Before describing the parquet again, check that it is still the
            # parquet the report describes. A report is only as trustworthy as
            # the file it points at.
            recorded_criterion = existing.get("criterion")
            if not isinstance(recorded_criterion, dict):
                raise ArtefactIntegrityError(
                    f"{report_path.name} records no criterion, so "
                    f"{dataset_path.name} cannot be re-hashed against it. "
                    "Refusing to refresh a report over an unverifiable artefact."
                )
            on_disk = content_hash(
                pd.read_parquet(dataset_path), recorded_criterion
            )
            if on_disk != existing.get("content_hash"):
                raise ArtefactIntegrityError(
                    f"{dataset_path.name} hashes to {on_disk} but "
                    f"{report_path.name} records "
                    f"{existing.get('content_hash')}. The parquet has changed "
                    "since its report was written, so the report may not be "
                    "refreshed over it. Investigate the file; do not re-run "
                    "over it."
                )
            report_path.write_text(
                json.dumps(
                    corpus.report, indent=2, sort_keys=True, ensure_ascii=False
                ),
                encoding="utf-8",
            )
            return {"dataset": dataset_path, "report": report_path}

    corpus.frame.to_parquet(dataset_path, engine="pyarrow", index=False)
    report_path.write_text(
        json.dumps(corpus.report, indent=2, sort_keys=True, ensure_ascii=False),
        encoding="utf-8",
    )
    return {"dataset": dataset_path, "report": report_path}


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------
SOURCE_KEY_IEDB_TCELL = "iedb_tcell"
SOURCE_KEY_PSEUDO = "mhc_pseudo"
SUPPORTED_SOURCE_KEYS: tuple[str, ...] = (SOURCE_KEY_IEDB_TCELL, SOURCE_KEY_PSEUDO)


def _file_digest(path: Path) -> dict:
    digest = hashlib.sha256()
    size = 0
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(io.DEFAULT_BUFFER_SIZE), b""):
            digest.update(chunk)
            size += len(chunk)
    return {"path": str(path), "sha256": digest.hexdigest(), "bytes": size}


def default_artefact_name(label_rule: str = DEFAULT_LABEL_RULE) -> str:
    """The artefact base name for a rule: one artefact per rule, by name."""
    if label_rule not in LABEL_RULES:
        raise UnknownLabelRuleError(
            f"unknown label rule {label_rule!r}; the pre-registered rules are "
            f"{list(LABEL_RULES)}"
        )
    return f"iedb_tcell_{label_rule}"


def harmonise(
    sources: dict[str, Path],
    out_dir: Path,
    name: str | None = None,
    label_rule: str = DEFAULT_LABEL_RULE,
) -> tuple[HarmonisedCorpus, dict[str, Path]]:
    """Build and write the canonical dataset. See FR-01 to FR-06.

    `sources` keys: `iedb_tcell` (the T cell export, .zip or .csv, required)
    and `mhc_pseudo` (an allele -> pseudo-sequence table, required). Any other
    key is refused: this pass ingests the IEDB T cell export only, and
    half-ingesting a corpus is worse than not ingesting it.

    `label_rule` selects one of the pre-registered pair-label rules and is
    part of the artefact's hashed criterion. When `name` is not given, the
    artefact is named after the rule, so one rule can never overwrite another.
    """
    unknown = sorted(set(sources) - set(SUPPORTED_SOURCE_KEYS))
    if unknown:
        raise UnsupportedSourceError(
            f"source(s) {unknown} are not ingested in this pass. Supported: "
            f"{list(SUPPORTED_SOURCE_KEYS)}. TESLA-derived and HiTIDE corpora "
            "are declared in the contract and still to be added."
        )
    for required in SUPPORTED_SOURCE_KEYS:
        if required not in sources:
            raise UnsupportedSourceError(
                f"source {required!r} is required; got {sorted(sources)}"
            )

    tcell_path = Path(sources[SOURCE_KEY_IEDB_TCELL])
    pseudo_path = Path(sources[SOURCE_KEY_PSEUDO])
    for path in (tcell_path, pseudo_path):
        if not path.exists():
            raise HarmonisationError(f"input file does not exist: {path}")

    registry = PseudoSequenceRegistry(
        [read_pseudo_sequence_dat(pseudo_path, source_name=f"file:{pseudo_path.name}")]
    )
    raw = read_iedb_tcell(tcell_path)
    corpus = harmonise_frame(
        raw,
        pseudo=registry,
        source_dataset=SOURCE_IEDB_TCELL,
        input_files={
            SOURCE_KEY_IEDB_TCELL: _file_digest(tcell_path),
            SOURCE_KEY_PSEUDO: _file_digest(pseudo_path),
        },
        label_rule=label_rule,
    )
    paths = write_corpus(
        corpus, out_dir=out_dir, name=name or default_artefact_name(label_rule)
    )
    return corpus, paths
