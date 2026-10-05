"""HLA class I nomenclature reconciliation (FR-02).

IEDB MHC restriction strings are heterogeneous: full molecular typings
(`HLA-A*02:01`), asterisk-free variants (`HLA-A02:01`), over-specified
typings (`HLA-A*02:01:01:02`), serological notations (`HLA-A2`, `HLA-Cw12`),
locus-only or class-only placeholders (`HLA class I`), class II molecules,
engineered mutants (`HLA-A*02:01 T163A mutant`) and non-human molecules
(`H2-Kb`).

This module turns one raw string into one explicit decision. There are only
three possible outcomes:

1.  a canonical two-field allele (`HLA-A*02:01`),
2.  a canonical allele group (`HLA-A*02`) when the input carries only
    two-digit resolution, which the caller may resolve to a representative
    only if that is unambiguous, or
3.  a drop, with a named reason.

Nothing is guessed. A string whose molecule cannot be determined is dropped,
never approximated, and every decision taken on the way is recorded in
`notes` so that it can be carried on the record itself.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

__all__ = [
    "AlleleResolution",
    "CLASSICAL_CLASS_I_LOCI",
    "NON_CLASSICAL_CLASS_I_LOCI",
    "ALLELE_DROP_REASONS",
    "reconcile_allele",
    "pseudo_sequence_key",
]

CLASSICAL_CLASS_I_LOCI = ("A", "B", "C")
"""The classical class I loci. Scope of this project."""

NON_CLASSICAL_CLASS_I_LOCI = ("E", "F", "G")
"""Non-classical class I loci: out of scope, dropped with their own reason."""

ALLELE_DROP_REASONS = (
    "allele_missing",
    "allele_engineered_variant",
    "allele_unspecified_resolution",
    "allele_serological_notation",
    "allele_non_human",
    "allele_not_class_i",
    "allele_non_classical_class_i",
    "allele_expression_suffix",
    "allele_unparseable",
)
"""Every reason this module may give for removing an allele."""

_ENGINEERED_MARKERS = ("mutant", "mutation", "chimer", "transgen", "variant of")
_CLASS_PLACEHOLDER = re.compile(r"(?i)^hla[\s-]*class\s*(i+|1|2)\b")
_LOCUS_SPLIT = re.compile(r"^([A-Za-z]+)(.*)$")
_DIGITS_ONLY = re.compile(r"^\d+$")


@dataclass(frozen=True)
class AlleleResolution:
    """The outcome of reconciling one raw MHC restriction string."""

    raw: str
    canonical: str | None = None
    locus: str | None = None
    group: str | None = None
    protein: str | None = None
    drop_reason: str | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def resolved(self) -> bool:
        return self.canonical is not None and self.drop_reason is None

    @property
    def group_level_only(self) -> bool:
        """True when the input carried only two-digit (allele group) resolution."""
        return self.resolved and self.protein is None


def pseudo_sequence_key(locus: str, group: str, protein: str) -> str:
    """The lookup key used by pseudo-sequence tables: `HLA-A02:01`, no asterisk."""
    return f"HLA-{locus}{group}:{protein}"


def _drop(raw: str, reason: str, *notes: str) -> AlleleResolution:
    return AlleleResolution(raw=raw, drop_reason=reason, notes=tuple(notes))


def reconcile_allele(raw: str | None) -> AlleleResolution:
    """Reconcile one raw MHC restriction string to a canonical class I allele.

    The checks are applied in a fixed, documented order so that the drop
    reason attached to a record is reproducible.
    """
    text = " ".join((raw or "").split())
    if not text:
        return _drop("", "allele_missing")

    lowered = text.lower()

    # 1. Engineered or otherwise annotated molecules are not the wild-type
    #    molecule and are not pooled with it.
    if any(marker in lowered for marker in _ENGINEERED_MARKERS):
        return _drop(text, "allele_engineered_variant")

    # 2. Class-level placeholders carry no typing information at all.
    if _CLASS_PLACEHOLDER.match(text):
        return _drop(text, "allele_unspecified_resolution")

    # 3. Human only.
    if not lowered.startswith("hla"):
        return _drop(text, "allele_non_human")

    remainder = text[len("hla"):].lstrip("-").strip()
    if not remainder:
        return _drop(text, "allele_unspecified_resolution")

    # An annotation in parentheses or after a space is not part of the typing.
    if " " in remainder or "(" in remainder:
        return _drop(text, "allele_unspecified_resolution")

    split = _LOCUS_SPLIT.match(remainder)
    if split is None:
        return _drop(text, "allele_unparseable")

    locus_raw, tail = split.group(1).upper(), split.group(2)
    notes: list[str] = []

    # 4. Locus. `Cw` is the superseded serological spelling of the C locus.
    if locus_raw in ("C", "CW"):
        locus = "C"
        if locus_raw == "CW":
            notes.append("locus_cw_normalised_to_c")
    elif locus_raw.startswith("D"):
        return _drop(text, "allele_not_class_i")
    elif locus_raw in NON_CLASSICAL_CLASS_I_LOCI:
        return _drop(text, "allele_non_classical_class_i")
    elif locus_raw in CLASSICAL_CLASS_I_LOCI:
        locus = locus_raw
    else:
        return _drop(text, "allele_unparseable")

    # 5. Typing fields.
    if not tail:
        return _drop(text, "allele_unspecified_resolution")

    if "*" not in tail and ":" not in tail:
        # Serological / low-resolution typing: `HLA-A2`, `HLA-Cw12`. One
        # serotype covers many molecules with different contact residues, so
        # this cannot be mapped to a representative without guessing.
        return _drop(text, "allele_serological_notation")

    fields = tail.lstrip("*").split(":")
    if any(not f for f in fields):
        return _drop(text, "allele_unparseable")

    # An expression suffix (N, Q, S, C, A) marks a molecule that is not a
    # normally expressed surface heterodimer.
    if fields[-1][-1].isalpha():
        return _drop(text, "allele_expression_suffix")

    if not all(_DIGITS_ONLY.match(f) for f in fields):
        return _drop(text, "allele_unparseable")

    group = fields[0].zfill(2)
    if group != fields[0]:
        notes.append("group_zero_padded")

    if len(fields) == 1:
        return AlleleResolution(
            raw=text,
            canonical=f"HLA-{locus}*{group}",
            locus=locus,
            group=group,
            protein=None,
            notes=tuple([*notes, "two_digit_resolution_only"]),
        )

    protein = fields[1].zfill(2)
    if len(fields) > 2:
        notes.append("truncated_to_two_field")

    return AlleleResolution(
        raw=text,
        canonical=f"HLA-{locus}*{group}:{protein}",
        locus=locus,
        group=group,
        protein=protein,
        notes=tuple(notes),
    )
