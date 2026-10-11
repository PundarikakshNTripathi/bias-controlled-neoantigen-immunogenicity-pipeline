"""Tests for Module 1, data harmonisation (FR-01 to FR-06).

Every fixture in this file is SYNTHETIC AND ILLUSTRATIVE. The peptides,
alleles and pseudo-sequences are hand-written to exercise one branch each;
none of them is a measured quantity and none is copied from a real export.
Row counts asserted here are properties of these illustrative fixtures, not
statistics about any corpus.

The real counts live only in the harmonisation report produced by running
`nbx data harmonise` on the real export.
"""

from __future__ import annotations

import ast
import json
from pathlib import Path

import pandas as pd
import pytest

from nbx.data.harmonise import (
    DEFAULT_LABEL_RULE,
    LABEL_RULES,
    PEPTIDE_MAX_LENGTH,
    PEPTIDE_MIN_LENGTH,
    PSEUDO_SEQUENCE_LENGTH,
    QUALITATIVE_LABEL_MAP,
    RAW_COLUMNS,
    ArtefactHashConflictError,
    ArtefactIntegrityError,
    DegenerateCorpusError,
    PseudoSequenceRegistry,
    PseudoSequenceTable,
    ReferenceStatistic,
    UnknownLabelRuleError,
    UnsupportedSourceError,
    compare_to_references,
    harmonise,
    harmonise_frame,
    write_corpus,
)
from nbx.evaluate.gate import UNMEASURED

# --------------------------------------------------------------------------
# Illustrative fixtures
# --------------------------------------------------------------------------

# Illustrative 34-residue strings. They are the right LENGTH on purpose and
# are otherwise meaningless; no real pseudo-sequence is reproduced here.
_SEQ_A0201 = "A" * PSEUDO_SEQUENCE_LENGTH
_SEQ_B0702 = "B" * PSEUDO_SEQUENCE_LENGTH
_SEQ_B35_SHARED = "C" * PSEUDO_SEQUENCE_LENGTH
_SEQ_B4402 = "D" * PSEUDO_SEQUENCE_LENGTH
_SEQ_B4403 = "E" * PSEUDO_SEQUENCE_LENGTH
_SEQ_C0701 = "F" * PSEUDO_SEQUENCE_LENGTH

ILLUSTRATIVE_PSEUDO_ENTRIES = {
    # Names follow the real file's convention: no asterisk.
    "HLA-A02:01": _SEQ_A0201,
    "HLA-B07:02": _SEQ_B0702,
    # Same two-digit group, identical sequence -> an unambiguous representative.
    "HLA-B35:01": _SEQ_B35_SHARED,
    "HLA-B35:02": _SEQ_B35_SHARED,
    # Same two-digit group, different sequences -> ambiguous, must not be guessed.
    "HLA-B44:02": _SEQ_B4402,
    "HLA-B44:03": _SEQ_B4403,
    "HLA-C07:01": _SEQ_C0701,
    "H2-Kb": "G" * PSEUDO_SEQUENCE_LENGTH,  # non-human entry, must not be admitted
}


def pseudo_registry() -> PseudoSequenceRegistry:
    table = PseudoSequenceTable(
        entries=ILLUSTRATIVE_PSEUDO_ENTRIES, source_name="illustrative_table"
    )
    return PseudoSequenceRegistry([table])


def raw_row(**overrides) -> dict[str, str]:
    """One IEDB-shaped row. Valid by default; override to break one thing."""
    row = {
        "Assay ID|IEDB IRI": "http://www.iedb.org/assay/illustrative-1",
        "Reference|PMID": "1",
        "Epitope|Object Type": "Linear peptide",
        "Epitope|Name": "SIINFEKLV",
        "Epitope|Modified residues": "",
        "Assay|Method": "illustrative IFNg release",
        "Assay|Response measured": "IFNg",
        "Assay|Qualitative Measurement": "Positive",
        "MHC Restriction|Name": "HLA-A*02:01",
        "MHC Restriction|Class": "I",
    }
    row.update(overrides)
    return row


def raw_frame(rows: list[dict[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=list(RAW_COLUMNS)).astype("string").fillna("")


def harmonise_rows(rows: list[dict[str, str]], label_rule: str = DEFAULT_LABEL_RULE):
    return harmonise_frame(
        raw_frame(rows), pseudo=pseudo_registry(), label_rule=label_rule
    )


def pair_labels(corpus) -> set[tuple[str, str, int]]:
    """The (peptide, allele, label) set a corpus asserts. The unit of analysis."""
    return {
        (str(r.peptide), str(r.allele), int(r.label))
        for r in corpus.frame.itertuples(index=False)
    }


# --------------------------------------------------------------------------
# FR-01  Ingest an IEDB-shaped export into one canonical frame.
# --------------------------------------------------------------------------
def test_fr01_canonical_frame_has_the_declared_schema_and_one_row_per_record():
    # FR-01: IEDB-shaped input becomes a canonical frame with a fixed schema.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLV", "MHC Restriction|Name": "HLA-A*02:01"}),
        raw_row(
            **{
                "Epitope|Name": "APRGPHGGAASGL",
                "MHC Restriction|Name": "HLA-B*07:02",
                "Assay|Qualitative Measurement": "Negative",
            }
        ),
    ]
    corpus = harmonise_rows(rows)

    for column in (
        "peptide",
        "allele",
        "allele_pseudoseq",
        "peptide_length",
        "label",
        "source_dataset",
        "assay_type",
        "negative_sample_origin",
        "allele_mapping_note",
        "contributing_sources",
        "n_contributing_sources",
        "n_contributing_records",
    ):
        assert column in corpus.frame.columns, column

    assert len(corpus.frame) == len(rows)
    assert set(corpus.frame["label"]) == {0, 1}
    assert (corpus.frame["peptide_length"] == corpus.frame["peptide"].str.len()).all()


def test_fr01_positive_levels_all_map_to_one_and_negative_maps_to_zero():
    # FR-01: the qualitative-measurement mapping is explicit and documented.
    assert set(QUALITATIVE_LABEL_MAP.values()) == {0, 1}
    assert QUALITATIVE_LABEL_MAP["Negative"] == 0
    positives = [k for k in QUALITATIVE_LABEL_MAP if k.startswith("Positive")]
    assert positives, "the mapping table must enumerate the Positive* levels"
    assert all(QUALITATIVE_LABEL_MAP[k] == 1 for k in positives)

    # Distinct peptides so that dedup does not merge them.
    peptides = ["SIINFEKLA", "SIINFEKLC", "SIINFEKLD", "SIINFEKLE", "SIINFEKLF"]
    rows = [
        raw_row(
            **{
                "Epitope|Name": peptide,
                "Assay|Qualitative Measurement": level,
            }
        )
        for peptide, level in zip(peptides, QUALITATIVE_LABEL_MAP)
    ]
    corpus = harmonise_rows(rows)
    observed = dict(
        zip(corpus.frame["peptide"], corpus.frame["label"], strict=False)
    )
    expected = {
        peptide: QUALITATIVE_LABEL_MAP[level]
        for peptide, level in zip(peptides, QUALITATIVE_LABEL_MAP)
    }
    assert observed == expected
    # The mapping must be visible in the report, not buried in the code.
    assert corpus.report["label_map"] == {
        k: int(v) for k, v in QUALITATIVE_LABEL_MAP.items()
    }


def test_fr01_non_peptidic_modified_and_non_class_i_records_are_dropped_by_reason():
    # FR-01: only linear, unmodified, class I human peptide records are admitted.
    rows = [
        raw_row(),
        raw_row(**{"Epitope|Object Type": "Non-peptidic"}),
        raw_row(**{"Epitope|Object Type": "Discontinuous protein residues"}),
        raw_row(**{"Epitope|Modified residues": "S2 phosphorylation"}),
        raw_row(**{"MHC Restriction|Class": "II"}),
        raw_row(**{"Assay|Qualitative Measurement": ""}),
        raw_row(**{"Epitope|Name": "SIINFEKLX"}),  # non-standard residue
    ]
    corpus = harmonise_rows(rows)
    drops = corpus.report["dropped_by_reason"]

    assert len(corpus.frame) == 1
    assert drops["epitope_not_linear_peptide"] == 2
    assert drops["epitope_has_modified_residues"] == 1
    assert drops["mhc_class_not_i"] == 1
    assert drops["label_not_mappable"] == 1
    assert drops["peptide_non_standard_residues"] == 1
    assert corpus.report["dropped_total"] == len(rows) - len(corpus.frame)


# --------------------------------------------------------------------------
# FR-02  Reconcile allele nomenclature; assign 34-residue pseudo-sequences.
# --------------------------------------------------------------------------
def test_fr02_nomenclature_variants_reconcile_to_one_canonical_allele():
    # FR-02: differing spellings of the same molecule become one allele string.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLA", "MHC Restriction|Name": "HLA-A*02:01"}),
        raw_row(**{"Epitope|Name": "SIINFEKLC", "MHC Restriction|Name": "hla-a*02:01"}),
        raw_row(**{"Epitope|Name": "SIINFEKLD", "MHC Restriction|Name": "HLA-A02:01"}),
        raw_row(
            **{"Epitope|Name": "SIINFEKLE", "MHC Restriction|Name": "HLA-A*02:01:01:02"}
        ),
    ]
    corpus = harmonise_rows(rows)
    assert set(corpus.frame["allele"]) == {"HLA-A*02:01"}
    assert len(corpus.frame) == len(rows)

    # The field-truncation decision is recorded on the record that needed it.
    note = corpus.frame.loc[
        corpus.frame["peptide"] == "SIINFEKLE", "allele_mapping_note"
    ].iloc[0]
    assert "truncated_to_two_field" in note


def test_fr02_pseudo_sequence_is_assigned_and_its_source_is_recorded_per_record():
    # FR-02: every retained record carries a 34-residue pseudo-sequence and
    # names the table that resolved it.
    corpus = harmonise_rows([raw_row()])
    row = corpus.frame.iloc[0]
    assert row["allele_pseudoseq"] == _SEQ_A0201
    assert len(row["allele_pseudoseq"]) == PSEUDO_SEQUENCE_LENGTH
    assert "illustrative_table" in row["allele_mapping_note"]
    assert corpus.frame["allele_pseudoseq"].str.len().eq(PSEUDO_SEQUENCE_LENGTH).all()
    assert corpus.report["pseudo_sequence_sources"] == ["illustrative_table"]


def test_fr02_unresolvable_and_out_of_scope_alleles_are_removed_with_counted_reasons():
    # FR-02: an allele that cannot be resolved is removed, and the reason is
    # counted. Nothing is guessed.
    rows = [
        raw_row(),
        raw_row(**{"MHC Restriction|Name": ""}),
        raw_row(**{"MHC Restriction|Name": "HLA class I"}),
        raw_row(**{"MHC Restriction|Name": "HLA-A2"}),  # serotype, many alleles
        raw_row(**{"MHC Restriction|Name": "HLA-Cw12"}),  # serological
        raw_row(**{"MHC Restriction|Name": "HLA-DR7"}),  # class II
        raw_row(**{"MHC Restriction|Name": "H2-Kb"}),  # mouse
        raw_row(**{"MHC Restriction|Name": "HLA-A*02:01 T163A mutant"}),
        raw_row(**{"MHC Restriction|Name": "HLA-A*02:01N"}),  # null allele
        raw_row(**{"MHC Restriction|Name": "HLA-E*01:01"}),  # non-classical
        raw_row(**{"MHC Restriction|Name": "HLA-C*18:01"}),  # whole group absent
    ]
    corpus = harmonise_rows(rows)
    drops = corpus.report["dropped_by_reason"]

    assert len(corpus.frame) == 1
    assert drops["allele_missing"] == 1
    assert drops["allele_unspecified_resolution"] == 1  # 'HLA class I'
    # 'HLA-A2' and 'HLA-Cw12' are serotypes: one serotype covers many
    # molecules with different contact residues, so neither can be resolved.
    assert drops["allele_serological_notation"] == 2
    assert drops["allele_not_class_i"] == 1
    assert drops["allele_non_human"] == 1
    assert drops["allele_engineered_variant"] == 1
    assert drops["allele_expression_suffix"] == 1
    assert drops["allele_non_classical_class_i"] == 1
    # 'HLA-C*18:01' resolves to a two-field name; the table simply does not
    # carry it. That is a coverage gap in the table, not an ambiguity.
    assert drops["allele_not_in_pseudo_sequence_table"] == 1
    assert sum(drops.values()) == len(rows) - len(corpus.frame)


def test_fr02_two_digit_representative_is_used_only_when_unambiguous():
    # FR-02: an allele with no sequence may take a same-two-digit representative
    # only when every peer in the group agrees; the choice is recorded per record.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLA", "MHC Restriction|Name": "HLA-B*35:43"}),
        raw_row(**{"Epitope|Name": "SIINFEKLC", "MHC Restriction|Name": "HLA-B*44:99"}),
    ]
    corpus = harmonise_rows(rows)

    assert len(corpus.frame) == 1
    kept = corpus.frame.iloc[0]
    assert kept["peptide"] == "SIINFEKLA"
    assert kept["allele"] == "HLA-B*35:43"  # identity is preserved
    assert kept["allele_pseudoseq"] == _SEQ_B35_SHARED
    note = kept["allele_mapping_note"]
    assert "two_digit_representative" in note
    assert "HLA-B*35:43" in note and "35" in note
    assert "illustrative_table" in note

    drops = corpus.report["dropped_by_reason"]
    # 'HLA-B*44:99' is fully resolved and absent from the table. Its group
    # carries two distinct contact sets, so no representative may be chosen -
    # but the record's problem is the table's coverage, not an ambiguous input.
    assert drops["allele_not_in_pseudo_sequence_table"] == 1
    assert "allele_two_digit_mapping_ambiguous" not in drops
    assert corpus.report["two_digit_representative_mappings"] == 1


def test_fr02_two_digit_input_in_an_ambiguous_group_is_its_own_drop_reason():
    # FR-02: the two reasons are distinguished. An input that carries only
    # two-digit resolution and whose group has several distinct contact sets is
    # genuinely ambiguous; an absent two-field allele is a table coverage gap.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLA", "MHC Restriction|Name": "HLA-B*44"}),
        raw_row(**{"Epitope|Name": "SIINFEKLC", "MHC Restriction|Name": "HLA-A*11"}),
        raw_row(**{"Epitope|Name": "SIINFEKLD", "MHC Restriction|Name": "HLA-B*35"}),
    ]
    corpus = harmonise_rows(rows)
    drops = corpus.report["dropped_by_reason"]

    assert drops["allele_two_digit_mapping_ambiguous"] == 1  # HLA-B*44
    assert drops["allele_pseudo_sequence_unavailable"] == 1  # HLA-A*11, no group
    assert "allele_not_in_pseudo_sequence_table" not in drops
    # HLA-B*35: two-digit input, group unanimous -> representative, retained.
    assert set(corpus.frame["allele"]) == {"HLA-B*35"}
    assert corpus.frame.iloc[0]["allele_pseudoseq"] == _SEQ_B35_SHARED


def test_fr02_report_lists_the_alleles_absent_from_the_pseudo_sequence_table():
    # FR-02: the table's coverage gap is visible, per allele, with counts.
    rows = [
        raw_row(),
        raw_row(**{"Epitope|Name": "SIINFEKLC", "MHC Restriction|Name": "HLA-C*18:01"}),
        raw_row(**{"Epitope|Name": "SIINFEKLD", "MHC Restriction|Name": "HLA-C*18:01"}),
        raw_row(**{"Epitope|Name": "SIINFEKLE", "MHC Restriction|Name": "HLA-B*44:99"}),
    ]
    corpus = harmonise_rows(rows)
    section = corpus.report["alleles_absent_from_pseudo_sequence_table"]

    assert section["n_alleles"] == 2
    assert section["n_records"] == 3
    by_allele = {entry["allele"]: entry["records"] for entry in section["alleles"]}
    assert by_allele == {"HLA-C*18:01": 2, "HLA-B*44:99": 1}
    # Listed most-affected first, so a coverage gap cannot hide in the tail.
    assert [entry["records"] for entry in section["alleles"]] == sorted(
        by_allele.values(), reverse=True
    )
    assert section["pseudo_sequence_sources"] == ["illustrative_table"]


# --------------------------------------------------------------------------
# FR-03  Length filter; deduplicate across overlapping sources.
# --------------------------------------------------------------------------
def test_fr03_peptides_outside_the_class_i_length_window_are_excluded_and_counted():
    # FR-03: 8-14 residues inclusive; exclusions are counted, not silent.
    rows = []
    for length in range(PEPTIDE_MIN_LENGTH - 2, PEPTIDE_MAX_LENGTH + 3):
        rows.append(
            raw_row(**{"Epitope|Name": "A" * (length - 1) + "K"})
        )
    corpus = harmonise_rows(rows)

    kept_lengths = sorted(corpus.frame["peptide_length"])
    assert kept_lengths == list(range(PEPTIDE_MIN_LENGTH, PEPTIDE_MAX_LENGTH + 1))

    n_too_short = PEPTIDE_MIN_LENGTH - (PEPTIDE_MIN_LENGTH - 2)
    n_too_long = (PEPTIDE_MAX_LENGTH + 2) - PEPTIDE_MAX_LENGTH
    drops = corpus.report["dropped_by_reason"]
    assert drops["peptide_too_short"] == n_too_short
    assert drops["peptide_too_long"] == n_too_long
    assert corpus.report["length_distribution"] == {
        str(length): 1 for length in kept_lengths
    }


def test_fr03_duplicate_peptide_allele_pairs_merge_keeping_the_union_of_provenance():
    # FR-03: deduplicate on (peptide, allele), keeping the union of provenance.
    rows = [
        raw_row(**{"Assay|Method": "assay alpha"}),
        raw_row(**{"Assay|Method": "assay beta", "Reference|PMID": "2"}),
        raw_row(**{"Assay|Method": "assay alpha"}),
    ]
    corpus = harmonise_rows(rows)

    assert len(corpus.frame) == 1
    row = corpus.frame.iloc[0]
    assert row["n_contributing_records"] == len(rows)
    assert "assay alpha" in row["assay_type"] and "assay beta" in row["assay_type"]
    assert row["n_contributing_sources"] == 1
    assert corpus.report["duplicate_records_merged"] == len(rows) - len(corpus.frame)


def test_fr03_strict_rule_drops_duplicates_that_disagree_on_label_and_counts_them():
    # FR-03: under the `strict` rule a (peptide, allele) pair with conflicting
    # labels is NOT silently resolved: the pair is dropped and counted.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLA"}),
        raw_row(
            **{"Epitope|Name": "SIINFEKLA", "Assay|Qualitative Measurement": "Negative"}
        ),
        raw_row(**{"Epitope|Name": "SIINFEKLC"}),
    ]
    corpus = harmonise_rows(rows, label_rule="strict")

    assert set(corpus.frame["peptide"]) == {"SIINFEKLC"}
    assert corpus.report["pairs_dropped_by_label_rule"] == 1
    assert corpus.report["dropped_by_reason"]["duplicate_label_conflict"] == 2
    dropped_pairs = corpus.report["dropped_pair_examples"]
    assert ["SIINFEKLA", "HLA-A*02:01"] in [list(c) for c in dropped_pairs]

    # The policy is not label-neutral, so its effect on the positive-to-negative
    # ratio is reported rather than left implicit.
    assert corpus.report["records_in_dropped_pairs_by_label"] == {"0": 1, "1": 1}
    assert corpus.report["label_counts_before_deduplication"] == {"0": 1, "1": 2}


# --------------------------------------------------------------------------
# FR-03 / analysis plan section 2: the (peptide, allele) PAIR is the unit of
# analysis, and the pre-registered label rules decide what a pair's label is.
# Every fixture below is synthetic and illustrative.
# --------------------------------------------------------------------------
_PEP_MAJORITY_POS = "MAYPASSPV"   # 2 positive, 1 negative
_PEP_MAJORITY_NEG = "MAYFAILVK"   # 1 positive, 2 negative
_PEP_TIED = "TIEDPAIRK"           # 1 positive, 1 negative
_PEP_UNANIMOUS_POS = "AGREEPSVK"  # 2 positive
_PEP_UNANIMOUS_NEG = "AGREENEGK"  # 1 negative
_ALLELE = "HLA-A*02:01"


def disagreeing_rows() -> list[dict[str, str]]:
    """Pairs with a deliberate assay pattern each. Synthetic, illustrative."""
    pos = "Positive"
    neg = "Negative"

    def row(peptide: str, level: str, assay: str) -> dict[str, str]:
        return raw_row(
            **{
                "Epitope|Name": peptide,
                "Assay|Qualitative Measurement": level,
                "Assay|Method": assay,
            }
        )

    return [
        row(_PEP_MAJORITY_POS, pos, "assay alpha"),
        row(_PEP_MAJORITY_POS, "Positive-Low", "assay beta"),
        row(_PEP_MAJORITY_POS, neg, "assay gamma"),
        row(_PEP_MAJORITY_NEG, pos, "assay alpha"),
        row(_PEP_MAJORITY_NEG, neg, "assay beta"),
        row(_PEP_MAJORITY_NEG, neg, "assay gamma"),
        row(_PEP_TIED, pos, "assay alpha"),
        row(_PEP_TIED, neg, "assay beta"),
        row(_PEP_UNANIMOUS_POS, pos, "assay alpha"),
        row(_PEP_UNANIMOUS_POS, "Positive-High", "assay beta"),
        row(_PEP_UNANIMOUS_NEG, neg, "assay alpha"),
    ]


def test_pair_is_the_unit_and_assay_counts_are_retained_per_pair():
    # Analysis plan section 2: rows sharing a pair are repeated measurements of
    # one observation, and the pair keeps its assay count and positive fraction.
    corpus = harmonise_rows(disagreeing_rows(), label_rule="any_positive")
    by_peptide = corpus.frame.set_index("peptide")

    assert len(corpus.frame) == 5  # five distinct pairs in the fixture
    assert not corpus.frame.duplicated(["peptide", "allele"]).any()
    assert by_peptide.loc[_PEP_MAJORITY_POS, "n_assay_records"] == 3
    assert by_peptide.loc[_PEP_MAJORITY_POS, "n_positive_assay_records"] == 2
    assert by_peptide.loc[_PEP_MAJORITY_POS, "positive_fraction"] == pytest.approx(
        2 / 3
    )
    assert by_peptide.loc[_PEP_UNANIMOUS_NEG, "n_assay_records"] == 1
    assert by_peptide.loc[_PEP_UNANIMOUS_NEG, "n_positive_assay_records"] == 0
    assert by_peptide.loc[_PEP_UNANIMOUS_NEG, "positive_fraction"] == pytest.approx(
        0.0
    )
    # The existing contributing-record count and the new assay count are the
    # same quantity seen from two sides; they may never disagree.
    assert (
        corpus.frame["n_assay_records"] == corpus.frame["n_contributing_records"]
    ).all()
    assert (
        corpus.frame["positive_fraction"]
        == corpus.frame["n_positive_assay_records"] / corpus.frame["n_assay_records"]
    ).all()


def test_majority_rule_labels_by_majority_and_drops_ties():
    # Analysis plan section 2: majority vote over the pair's assays, ties
    # dropped under a named reason. All Positive* levels count as positive.
    corpus = harmonise_rows(disagreeing_rows(), label_rule="majority")

    assert pair_labels(corpus) == {
        (_PEP_MAJORITY_POS, _ALLELE, 1),
        (_PEP_MAJORITY_NEG, _ALLELE, 0),
        (_PEP_UNANIMOUS_POS, _ALLELE, 1),
        (_PEP_UNANIMOUS_NEG, _ALLELE, 0),
    }
    assert _PEP_TIED not in set(corpus.frame["peptide"])
    assert corpus.report["dropped_by_reason"]["pair_label_tie"] == 2
    assert corpus.report["pairs_dropped_by_reason"] == {"pair_label_tie": 1}
    assert "duplicate_label_conflict" not in corpus.report["dropped_by_reason"]
    assert corpus.report["label_rule"] == "majority"


def test_strict_rule_drops_every_pair_whose_assay_records_disagree():
    # Analysis plan section 2(a): the strict sensitivity analysis. Recomputed
    # here independently of the module's own grouping.
    rows = disagreeing_rows()
    corpus = harmonise_rows(rows, label_rule="strict")

    observed: dict[tuple[str, str], set[int]] = {}
    for row in rows:
        key = (row["Epitope|Name"], _ALLELE)
        observed.setdefault(key, set()).add(
            QUALITATIVE_LABEL_MAP[row["Assay|Qualitative Measurement"]]
        )
    expected = {
        (peptide, allele, next(iter(labels)))
        for (peptide, allele), labels in observed.items()
        if len(labels) == 1
    }
    assert pair_labels(corpus) == expected
    assert corpus.report["pairs_dropped_by_reason"] == {"duplicate_label_conflict": 3}
    assert corpus.report["dropped_by_reason"]["duplicate_label_conflict"] == 8


def test_any_positive_rule_labels_a_pair_positive_if_any_assay_is_positive():
    # Analysis plan section 2(c): the permissive sensitivity analysis. No pair
    # is dropped by this rule.
    corpus = harmonise_rows(disagreeing_rows(), label_rule="any_positive")

    assert pair_labels(corpus) == {
        (_PEP_MAJORITY_POS, _ALLELE, 1),
        (_PEP_MAJORITY_NEG, _ALLELE, 1),
        (_PEP_TIED, _ALLELE, 1),
        (_PEP_UNANIMOUS_POS, _ALLELE, 1),
        (_PEP_UNANIMOUS_NEG, _ALLELE, 0),
    }
    assert corpus.report["pairs_dropped_by_label_rule"] == 0
    assert corpus.report["pairs_dropped_by_reason"] == {}


def test_strict_rule_reproduces_the_previous_conflict_dropping_behaviour():
    # Regression guard, synthetic twin of the real-corpus comparison: `strict`
    # keeps exactly the pairs whose assay records all agree, which is what the
    # module did before the rules existed.
    rows = disagreeing_rows()
    strict = harmonise_rows(rows, label_rule="strict")
    majority = harmonise_rows(rows, label_rule="majority")
    permissive = harmonise_rows(rows, label_rule="any_positive")

    assert pair_labels(strict) <= pair_labels(permissive)
    assert pair_labels(strict) != pair_labels(majority)
    # Unanimous pairs are labelled identically by all three rules.
    unanimous = {
        (_PEP_UNANIMOUS_POS, _ALLELE, 1),
        (_PEP_UNANIMOUS_NEG, _ALLELE, 0),
    }
    for corpus in (strict, majority, permissive):
        assert unanimous <= pair_labels(corpus)


def test_label_rule_is_part_of_the_hashed_criterion():
    # Each rule is a separate artefact with its own content hash (section 2).
    rows = disagreeing_rows()
    hashes = {}
    for rule in LABEL_RULES:
        corpus = harmonise_rows(rows, label_rule=rule)
        assert corpus.report["criterion"]["label_rule"] == rule
        hashes[rule] = corpus.content_hash
    assert len(set(hashes.values())) == len(LABEL_RULES)
    assert DEFAULT_LABEL_RULE == "majority"


def test_unknown_label_rule_is_refused():
    with pytest.raises(UnknownLabelRuleError):
        harmonise_rows([raw_row()], label_rule="whichever_wins")


def test_pair_pattern_summary_counts_come_from_the_pairs():
    # The disagreement pattern is reported, per rule, from the pairs - never
    # from a constant. Under `any_positive` no pair is dropped, so the summary
    # must be reproducible from the frame alone.
    corpus = harmonise_rows(disagreeing_rows(), label_rule="any_positive")
    frame = corpus.frame
    summary = corpus.report["pair_pattern_summary"]

    n_pos = frame["n_positive_assay_records"]
    n_neg = frame["n_assay_records"] - frame["n_positive_assay_records"]
    assert summary["pairs_total"] == len(frame)
    assert summary["unanimous_positive"] == int((n_neg == 0).sum())
    assert summary["unanimous_negative"] == int((n_pos == 0).sum())
    assert summary["mixed"] == int(((n_pos > 0) & (n_neg > 0)).sum())
    assert summary["tied"] == int((n_pos == n_neg).sum())
    assert (
        summary["unanimous_positive"]
        + summary["unanimous_negative"]
        + summary["mixed"]
        == summary["pairs_total"]
    )
    assert summary["tied"] <= summary["mixed"]
    assert summary["note"]


def test_reconciliation_balances_at_record_and_pair_level_for_every_rule():
    # Kept + dropped (by reason) + collapsed-into-pair must equal the raw row
    # count, and the pair counts must balance too.
    rows = disagreeing_rows() + [
        raw_row(**{"Epitope|Name": "SIINFEKLX"}),  # non-standard residue
        raw_row(**{"Epitope|Name": "SIINFEKLA", "Assay|Method": ""}),  # no provenance
        raw_row(**{"MHC Restriction|Class": "II"}),  # out of scope
    ]
    for rule in LABEL_RULES:
        corpus = harmonise_rows(rows, label_rule=rule)
        rec = corpus.report["reconciliation"]

        assert rec["raw_row_count"] == len(rows)
        assert (
            rec["records_dropped_record_level"]
            + rec["records_dropped_pair_level"]
            + rec["records_collapsed_into_pairs"]
            + rec["pairs_labelled"]
            == rec["raw_row_count"]
        ), rule
        assert rec["pairs_labelled"] == len(corpus.frame)
        assert rec["pairs_total"] == rec["pairs_labelled"] + rec["pairs_dropped"]
        assert rec["record_level_balances"] is True
        assert rec["pair_level_balances"] is True
        assert rec["records_in_labelled_pairs"] == int(
            corpus.frame["n_assay_records"].sum()
        )
        # The record-level drop total still accounts for every dropped row.
        assert corpus.report["dropped_total"] == (
            rec["records_dropped_record_level"] + rec["records_dropped_pair_level"]
        )


# --------------------------------------------------------------------------
# FR-04  Provenance is mandatory per record; degenerate input must raise.
# --------------------------------------------------------------------------
def test_fr04_record_missing_provenance_is_dropped_not_defaulted():
    # FR-04 (load-bearing): a blank assay method or assay identifier means the
    # record has lost its provenance. It is dropped and counted; no default is
    # substituted.
    rows = [
        raw_row(),
        raw_row(**{"Epitope|Name": "SIINFEKLC", "Assay|Method": ""}),
        raw_row(**{"Epitope|Name": "SIINFEKLD", "Assay|Method": "   "}),
        raw_row(**{"Epitope|Name": "SIINFEKLE", "Assay ID|IEDB IRI": ""}),
    ]
    corpus = harmonise_rows(rows)

    assert len(corpus.frame) == 1
    assert corpus.report["dropped_by_reason"]["missing_provenance"] == 3
    assert corpus.report["missing_provenance_fields"] == {
        "assay_type": 2,
        "assay_iri": 1,
    }
    # Nothing that looks like a filled-in default may appear in the frame.
    for column in ("assay_type", "source_dataset", "negative_sample_origin"):
        values = corpus.frame[column].tolist()
        assert all(v.strip() for v in values)
        assert not any(
            v.strip().lower() in {"unknown", "n/a", "na", "none", "null", "0"}
            for v in values
        )


def test_fr04_every_retained_record_carries_the_full_provenance_triple():
    # FR-04: source dataset, assay type and negative-sample origin on every row.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLA"}),
        raw_row(
            **{"Epitope|Name": "SIINFEKLC", "Assay|Qualitative Measurement": "Negative"}
        ),
    ]
    corpus = harmonise_rows(rows)

    assert set(corpus.frame["source_dataset"]) == {corpus.report["source_dataset"]}
    positives = corpus.frame[corpus.frame["label"] == 1]
    negatives = corpus.frame[corpus.frame["label"] == 0]
    # Negatives say where they came from; positives carry an explicit
    # non-applicable marker rather than an empty cell.
    assert negatives["negative_sample_origin"].eq("iedb_assay_negative").all()
    assert positives["negative_sample_origin"].str.startswith("not_applicable").all()
    assert corpus.frame["allele_mapping_note"].str.len().gt(0).all()


def test_fr04_empty_input_raises_rather_than_returning_an_empty_corpus():
    # FR-04: degenerate input must raise, never return silently.
    with pytest.raises(DegenerateCorpusError):
        harmonise_rows([])


def test_fr04_input_whose_every_record_is_dropped_raises():
    # FR-04: an all-dropped corpus is degenerate too.
    with pytest.raises(DegenerateCorpusError):
        harmonise_rows([raw_row(**{"MHC Restriction|Class": "II"})])


def test_fr04_unknown_source_key_is_refused():
    # FR-04: this pass ingests the IEDB T cell export only; anything else is
    # refused loudly rather than half-ingested.
    with pytest.raises(UnsupportedSourceError):
        harmonise({"tesla": Path("does-not-exist.csv")}, out_dir=Path("."))


# --------------------------------------------------------------------------
# FR-05  Content-hashed, frozen canonical dataset.
# --------------------------------------------------------------------------
def test_fr05_identical_content_gives_an_identical_hash():
    # FR-05: the hash is over content, so it is reproducible.
    rows = [raw_row(**{"Epitope|Name": "SIINFEKLA"}), raw_row()]
    first = harmonise_rows(rows)
    second = harmonise_rows(list(reversed(rows)))
    assert first.content_hash == second.content_hash
    assert first.content_hash.startswith("sha256:") is False  # bare hex digest
    assert len(first.content_hash) == 64


def test_fr05_changed_content_gives_a_different_hash():
    # FR-05: a changed input must not reuse an existing identity.
    base = harmonise_rows([raw_row()])
    changed_label = harmonise_rows(
        [raw_row(**{"Assay|Qualitative Measurement": "Negative"})]
    )
    changed_peptide = harmonise_rows([raw_row(**{"Epitope|Name": "SIINFEKLC"})])
    assert len({base.content_hash, changed_label.content_hash,
                changed_peptide.content_hash}) == 3


def test_fr05_write_refuses_to_overwrite_an_artefact_with_different_content(tmp_path):
    # FR-05: a frozen artefact is never silently replaced.
    first = harmonise_rows([raw_row()])
    paths = write_corpus(first, out_dir=tmp_path, name="illustrative")
    assert paths["dataset"].exists() and paths["report"].exists()
    report = json.loads(paths["report"].read_text(encoding="utf-8"))
    assert report["content_hash"] == first.content_hash

    # Re-writing identical content is a no-op, not an error.
    again = write_corpus(first, out_dir=tmp_path, name="illustrative")
    assert again["dataset"] == paths["dataset"]

    different = harmonise_rows([raw_row(**{"Epitope|Name": "SIINFEKLC"})])
    with pytest.raises(ArtefactHashConflictError):
        write_corpus(different, out_dir=tmp_path, name="illustrative")


def test_fr05_written_dataset_round_trips_with_the_same_content_hash(tmp_path):
    # FR-05: the hash identifies what was actually written.
    corpus = harmonise_rows([raw_row(), raw_row(**{"Epitope|Name": "SIINFEKLC"})])
    paths = write_corpus(corpus, out_dir=tmp_path, name="illustrative")
    reloaded = pd.read_parquet(paths["dataset"])
    assert len(reloaded) == len(corpus.frame)
    assert list(reloaded.columns) == list(corpus.frame.columns)

    # Re-writing the same corpus re-hashes the parquet on disk and agrees.
    again = write_corpus(corpus, out_dir=tmp_path, name="illustrative")
    assert again["dataset"] == paths["dataset"]


def test_fr05_write_refuses_when_the_parquet_no_longer_matches_its_report(tmp_path):
    # FR-05: refreshing a report is only legitimate if the parquet is still the
    # corpus the report describes. An edited parquet is caught, not described.
    corpus = harmonise_rows([raw_row(), raw_row(**{"Epitope|Name": "SIINFEKLC"})])
    paths = write_corpus(corpus, out_dir=tmp_path, name="illustrative")

    tampered = pd.read_parquet(paths["dataset"])
    tampered.loc[0, "label"] = 1 - int(tampered.loc[0, "label"])
    tampered.to_parquet(paths["dataset"], engine="pyarrow", index=False)

    with pytest.raises(ArtefactIntegrityError):
        write_corpus(corpus, out_dir=tmp_path, name="illustrative")


# --------------------------------------------------------------------------
# FR-06  Statistics come from the frame; references are compared, not assumed.
# --------------------------------------------------------------------------
def test_fr06_every_reported_statistic_is_recomputed_from_the_frame():
    # FR-06: no count in the report may come from anywhere but the frame.
    rows = [
        raw_row(**{"Epitope|Name": "SIINFEKLA"}),
        raw_row(
            **{
                "Epitope|Name": "APRGPHGGAASGL",
                "MHC Restriction|Name": "HLA-B*07:02",
                "Assay|Qualitative Measurement": "Negative",
            }
        ),
        raw_row(**{"Epitope|Name": "SIINFEKLLL", "MHC Restriction|Name": "HLA-C*07:01"}),
    ]
    corpus = harmonise_rows(rows)
    frame = corpus.frame
    report = corpus.report

    assert report["record_count"] == len(frame)
    assert report["allele_count"] == frame["allele"].nunique()
    assert report["label_counts"] == {
        str(k): int(v) for k, v in frame["label"].value_counts().items()
    }
    assert report["length_distribution"] == {
        str(k): int(v) for k, v in frame["peptide_length"].value_counts().items()
    }
    per_source = report["per_source_breakdown"]
    assert set(per_source) == set(frame["source_dataset"].unique())
    for source, stats in per_source.items():
        subset = frame[frame["source_dataset"] == source]
        assert stats["record_count"] == len(subset)
        assert stats["allele_count"] == subset["allele"].nunique()
        assert stats["n_positive"] == int((subset["label"] == 1).sum())
        assert stats["n_negative"] == int((subset["label"] == 0).sum())


def test_fr06_a_statistic_with_no_located_reference_renders_as_an_em_dash():
    # FR-06: where no published reference exists, the comparison reports an em
    # dash. It never invents a reference and never quietly claims a match.
    corpus = harmonise_rows([raw_row()])
    comparison = corpus.report["reference_comparison"]
    assert comparison["references_located"] == 0
    assert comparison["note"]
    for entry in comparison["statistics"]:
        assert entry["observed"] != UNMEASURED
        assert entry["reference"] == UNMEASURED
        assert entry["discrepancy"] == UNMEASURED
        assert entry["status"] == "no_reference_located"
        assert entry["locator"] == UNMEASURED


def test_fr06_a_located_reference_is_compared_and_a_discrepancy_is_documented():
    # FR-06: the comparison mechanism works when a reference with a locator
    # exists; a mismatch is reported as a discrepancy, not rounded away.
    # Both values below are illustrative, not measured.
    observed = {"record_count": 10.0, "allele_count": 2.0}
    references = (
        ReferenceStatistic(
            name="record_count", value=12.0, locator="illustrative source, p. 1"
        ),
        ReferenceStatistic(
            name="allele_count", value=2.0, locator="illustrative source, p. 2"
        ),
    )
    result = compare_to_references(observed, references)
    by_name = {entry["name"]: entry for entry in result["statistics"]}
    assert result["references_located"] == len(references)
    assert by_name["record_count"]["status"] == "discrepancy"
    assert by_name["record_count"]["discrepancy"] != UNMEASURED
    assert by_name["allele_count"]["status"] == "reproduced"
    assert by_name["allele_count"]["locator"] == "illustrative source, p. 2"


def test_fr06_no_published_statistic_constant_lives_in_the_source():
    # FR-06: corpus-scale numbers must never be hard-coded. Any integer that
    # could plausibly be a record or allele count is forbidden in the module.
    source_path = Path(__file__).resolve().parents[1] / "src" / "nbx" / "data"
    for module_path in sorted(source_path.glob("*.py")):
        tree = ast.parse(module_path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, int):
                assert abs(node.value) < 1000, (
                    f"{module_path.name}: literal {node.value} looks like a "
                    "corpus statistic; counts must come from the frame"
                )
            if isinstance(node, ast.Constant) and isinstance(node.value, float):
                assert abs(node.value) < 1000, module_path.name
