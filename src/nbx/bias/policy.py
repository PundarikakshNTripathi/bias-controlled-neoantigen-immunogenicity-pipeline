"""Adaptive, allele-aware masking policy (Section 7.6.2, objective O3).

The published mitigation masks peptide position 2 and the C-terminus
uniformly across all alleles, and its authors state this as a simplification
because anchor positions vary by allele. This module generalises it:

    M(a) = { j : I( x_j ; y | a ) > tau }                              (7.9)

The reduction property is the whole point of the comparison, and is tested:
when the informative position set is identical across alleles, the adaptive
policy produces exactly the uniform policy. That makes the ablation clean --
a difference in results can then be attributed to per-allele adaptation and
to nothing else.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Hashable, Mapping, Sequence
from dataclasses import dataclass

from nbx.bias.diagnostics import mutual_information

__all__ = [
    "MaskingPolicy",
    "uniform_policy",
    "per_position_mi",
    "adaptive_policy",
    "select_threshold",
]


@dataclass(frozen=True)
class MaskingPolicy:
    """Which peptide positions are masked, per allele.

    Positions are zero-indexed. `default` applies to any allele absent from
    `per_allele`, so a policy is total over the allele space rather than
    silently masking nothing for an unseen allele.
    """

    per_allele: Mapping[str, frozenset[int]]
    default: frozenset[int]
    peptide_length: int
    name: str = "adaptive"

    def mask_for(self, allele: Hashable) -> frozenset[int]:
        return self.per_allele.get(str(allele), self.default)

    def is_uniform(self) -> bool:
        """True when every allele receives the same mask."""
        masks = set(self.per_allele.values())
        if not masks:
            return True
        return len(masks) == 1 and next(iter(masks)) == self.default

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "peptide_length": self.peptide_length,
            "default": sorted(self.default),
            "per_allele": {a: sorted(m) for a, m in sorted(self.per_allele.items())},
            "is_uniform": self.is_uniform(),
        }


def uniform_policy(peptide_length: int = 9, alleles: Sequence[str] = ()) -> MaskingPolicy:
    """The published scheme: mask position 2 and the C-terminus, everywhere.

    Position 2 is index 1; the C-terminus is the last index. These are the
    primary anchor positions for most HLA class I alleles, which is exactly
    why masking them uniformly is an approximation rather than a rule.
    """
    if peptide_length < 3:
        raise ValueError("peptide_length must be at least 3")
    mask = frozenset({1, peptide_length - 1})
    return MaskingPolicy(
        per_allele={str(a): mask for a in alleles},
        default=mask,
        peptide_length=peptide_length,
        name="uniform",
    )


def per_position_mi(
    peptides: Sequence[Sequence[Hashable]],
    labels: Sequence[int],
    alleles: Sequence[Hashable],
    peptide_length: int,
) -> dict[str, list[float]]:
    """I(x_j ; y | a) for every position j, computed within each allele.

    Returns a mapping from allele to a list of per-position mutual
    information values, one per peptide position. Alleles whose subset has a
    constant label contribute a list of zeros: with no label variation there
    is no information for any position to carry.
    """
    if not (len(peptides) == len(labels) == len(alleles)):
        raise ValueError("peptides, labels and alleles must be parallel sequences")

    groups: dict[str, list[int]] = defaultdict(list)
    for i, a in enumerate(alleles):
        groups[str(a)].append(i)

    out: dict[str, list[float]] = {}
    for allele, idx in groups.items():
        ys = [int(labels[i]) for i in idx]
        if len(set(ys)) < 2:
            out[allele] = [0.0] * peptide_length
            continue
        values = []
        for j in range(peptide_length):
            residues = [peptides[i][j] for i in idx]
            values.append(mutual_information(ys, residues))
        out[allele] = values
    return out


def select_threshold(
    mi_by_allele: Mapping[str, Sequence[float]],
    quantile: float = 0.90,
) -> float:
    """Pick tau as a quantile of the pooled per-position MI distribution.

    The specification says tau is selected on validation data. This provides
    a defensible default rule; the caller is expected to pass validation-set
    MI, not training-set MI, and the choice is recorded in the policy
    artefact so it can be audited later.
    """
    if not 0.0 < quantile < 1.0:
        raise ValueError("quantile must lie strictly between 0 and 1")
    pooled = sorted(v for vals in mi_by_allele.values() for v in vals)
    if not pooled:
        raise ValueError("no mutual-information values to threshold")
    k = int(quantile * (len(pooled) - 1))
    return pooled[k]


def adaptive_policy(
    mi_by_allele: Mapping[str, Sequence[float]],
    tau: float,
    peptide_length: int,
    restrict_to: frozenset[int] | None = None,
) -> MaskingPolicy:
    """Build M(a) = { j : I(x_j ; y | a) > tau } for each allele.

    Parameters
    ----------
    restrict_to
        When given, only positions in this set may be masked. The intended
        use is to restrict masking to candidate anchor positions, so that
        positions known to carry T cell contact signal are never masked away
        merely because they happen to be informative in a biased corpus.
        Section 7.6.2 flags this as a real cost of masking.
    """
    per_allele: dict[str, frozenset[int]] = {}
    for allele, values in mi_by_allele.items():
        if len(values) != peptide_length:
            raise ValueError(
                f"allele {allele!r} has {len(values)} MI values, "
                f"expected {peptide_length}"
            )
        selected = {j for j, v in enumerate(values) if v > tau}
        if restrict_to is not None:
            selected &= set(restrict_to)
        per_allele[allele] = frozenset(selected)

    fallback = frozenset({1, peptide_length - 1})
    if restrict_to is not None:
        fallback &= set(restrict_to)

    return MaskingPolicy(
        per_allele=per_allele,
        default=frozenset(fallback),
        peptide_length=peptide_length,
        name="adaptive",
    )
