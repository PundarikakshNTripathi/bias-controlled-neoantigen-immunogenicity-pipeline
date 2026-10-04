"""Bias diagnostics: the integrity core of this project.

Implements, exactly as specified in Section 7.4 of the synopsis:

    I(b ; y) = H(y) - H(y | b)                                        (7.1)
    H(y)     = -sum_y P(y) log P(y)                                   (7.2)
    H(y | b) = -sum_b P(b) sum_y P(y|b) log P(y|b)                    (7.3)
    rho(a)   = log2( N+(a) / N-(a) )                                  (7.4)
    kappa(f) = spearman( rho(a), E[f(x) | a] )                        (7.5)

These four quantities are the project's headline reported measurements. They
are implemented here, with tests that pin them against analytically known
cases, so that every downstream claim rests on arithmetic that has been
checked rather than on arithmetic that merely ran.

Conventions fixed by the specification
--------------------------------------
* Logarithms are natural (nats) for entropy and mutual information. The unit
  is reported alongside every value; `log_base` lets a caller ask for bits.
* rho(a) saturates at +/- RHO_SATURATION for alleles possessing only positive
  or only negative examples, rather than returning infinity.
* Alleles represented by a single sample are excluded from rho and from kappa.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from collections.abc import Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field

__all__ = [
    "RHO_SATURATION",
    "entropy",
    "conditional_entropy",
    "mutual_information",
    "log_imbalance_ratio",
    "per_allele_imbalance",
    "spearman",
    "shortcut_dependence",
    "BiasProfile",
    "profile_dataset",
]

# Alleles with only one label class saturate here instead of diverging.
# log2(10) = 3.3219... is the ">10:1" threshold referenced in the documents,
# so a saturation well above it keeps saturated alleles clearly separated
# from merely-imbalanced ones.
RHO_SATURATION: float = 10.0

_EXCLUDE_SINGLE_SAMPLE_ALLELES = True


def _as_list(xs: Iterable) -> list:
    return xs if isinstance(xs, list) else list(xs)


# --------------------------------------------------------------------------
# Information theory
# --------------------------------------------------------------------------
def entropy(labels: Sequence[Hashable], log_base: float | None = None) -> float:
    """Shannon entropy H(y) of an observed label sequence.

    Equation (7.2). Empty input is an error rather than zero: an entropy of
    zero is a measurement ("the label never varies"), and silently returning
    it for "no data" is exactly the sort of placeholder this project forbids.
    """
    xs = _as_list(labels)
    if not xs:
        raise ValueError("entropy() requires at least one observation")

    n = len(xs)
    h = 0.0
    for count in Counter(xs).values():
        p = count / n
        h -= p * math.log(p)
    return h if log_base is None else h / math.log(log_base)


def conditional_entropy(
    labels: Sequence[Hashable],
    by: Sequence[Hashable],
    log_base: float | None = None,
) -> float:
    """Conditional entropy H(y | b).

    Equation (7.3). `labels` and `by` are parallel sequences: `by[i]` is the
    value of the conditioning variable for observation `labels[i]`.
    """
    ys = _as_list(labels)
    bs = _as_list(by)
    if len(ys) != len(bs):
        raise ValueError(
            f"labels and by must be parallel sequences; got {len(ys)} and {len(bs)}"
        )
    if not ys:
        raise ValueError("conditional_entropy() requires at least one observation")

    n = len(ys)
    groups: dict[Hashable, list] = defaultdict(list)
    for b, y in zip(bs, ys):
        groups[b].append(y)

    h = 0.0
    for group in groups.values():
        p_b = len(group) / n
        h += p_b * entropy(group)
    return h if log_base is None else h / math.log(log_base)


def mutual_information(
    labels: Sequence[Hashable],
    by: Sequence[Hashable],
    log_base: float | None = None,
) -> float:
    """Mutual information I(b ; y) between a candidate bias variable and the label.

    Equation (7.1), in the marginal-conditional form the specification uses.

    A high value for a variable with no direct biological causal relationship
    to the label marks that variable as a bias candidate. The function does
    not make that judgement: it reports the number.

    The result is clamped at zero from below. Mutual information is
    non-negative in exact arithmetic, so a small negative value can only be
    floating-point noise; a large one would be a bug and is left to surface.
    """
    mi = entropy(labels, log_base) - conditional_entropy(labels, by, log_base)
    return 0.0 if -1e-12 < mi < 0.0 else mi


# --------------------------------------------------------------------------
# Intra-allelic imbalance
# --------------------------------------------------------------------------
def log_imbalance_ratio(n_pos: int, n_neg: int) -> float:
    """rho(a) = log2(N+ / N-), saturating when one class is absent.

    Equation (7.4). Saturation, rather than infinity, keeps the value usable
    in the rank correlation of `shortcut_dependence` without special-casing.
    """
    if n_pos < 0 or n_neg < 0:
        raise ValueError("counts must be non-negative")
    if n_pos == 0 and n_neg == 0:
        raise ValueError("allele has no samples; it should have been excluded")
    if n_neg == 0:
        return RHO_SATURATION
    if n_pos == 0:
        return -RHO_SATURATION
    return math.log2(n_pos / n_neg)


def per_allele_imbalance(
    alleles: Sequence[Hashable],
    labels: Sequence[int],
    exclude_single_sample: bool = _EXCLUDE_SINGLE_SAMPLE_ALLELES,
) -> dict[Hashable, float]:
    """rho(a) for every allele in the corpus.

    `labels` are binary, with 1 meaning positive. Alleles represented by a
    single sample are excluded by default, as the specification requires:
    their rho is always saturated and carries no information about imbalance.
    """
    als = _as_list(alleles)
    ys = _as_list(labels)
    if len(als) != len(ys):
        raise ValueError(
            f"alleles and labels must be parallel; got {len(als)} and {len(ys)}"
        )

    pos: Counter = Counter()
    neg: Counter = Counter()
    for a, y in zip(als, ys):
        if y not in (0, 1, True, False):
            raise ValueError(f"labels must be binary; saw {y!r}")
        (pos if int(y) == 1 else neg)[a] += 1

    out: dict[Hashable, float] = {}
    for a in set(als):
        p, n = pos[a], neg[a]
        if exclude_single_sample and (p + n) < 2:
            continue
        out[a] = log_imbalance_ratio(p, n)
    return out


# --------------------------------------------------------------------------
# Shortcut dependence
# --------------------------------------------------------------------------
def _rank(values: Sequence[float]) -> list[float]:
    """Ranks with ties averaged, as Spearman's rho requires."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def spearman(xs: Sequence[float], ys: Sequence[float]) -> float:
    """Spearman rank correlation, with ties averaged.

    Implemented here rather than imported so that the project's headline
    diagnostic has no dependency that could change its definition underneath
    it. Verified against scipy.stats.spearmanr in the test suite.
    """
    xs = _as_list(xs)
    ys = _as_list(ys)
    if len(xs) != len(ys):
        raise ValueError("inputs must be the same length")
    if len(xs) < 3:
        raise ValueError(
            "rank correlation over fewer than 3 alleles is not interpretable"
        )

    rx, ry = _rank(xs), _rank(ys)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0.0 or dy == 0.0:
        raise ValueError(
            "a rank correlation is undefined when one variable is constant"
        )
    return num / (dx * dy)


def shortcut_dependence(
    train_imbalance: Mapping[Hashable, float],
    test_alleles: Sequence[Hashable],
    test_scores: Sequence[float],
) -> float:
    """kappa(f): the project's headline diagnostic. Equation (7.5).

    Rank correlation between an allele's imbalance in the TRAINING corpus and
    the mean score the model assigns that allele's TEST peptides.

    A model that has learned peptide-level biology rather than allele label
    statistics should sit near zero. A clearly positive value means the model
    scores an allele according to how that allele was labelled in training.

    Only alleles present in both the training imbalance profile and the test
    set contribute. The asymmetry is deliberate: an allele the model never
    trained on has no training imbalance to correlate against.
    """
    als = _as_list(test_alleles)
    sc = _as_list(test_scores)
    if len(als) != len(sc):
        raise ValueError("test_alleles and test_scores must be parallel")

    sums: dict[Hashable, float] = defaultdict(float)
    counts: Counter = Counter()
    for a, s in zip(als, sc):
        if a in train_imbalance:
            sums[a] += float(s)
            counts[a] += 1

    shared = sorted(counts, key=str)
    if len(shared) < 3:
        raise ValueError(
            f"kappa needs at least 3 alleles shared between train and test; "
            f"found {len(shared)}"
        )

    rho_vals = [train_imbalance[a] for a in shared]
    mean_scores = [sums[a] / counts[a] for a in shared]
    return spearman(rho_vals, mean_scores)


# --------------------------------------------------------------------------
# Dataset-level profile
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class BiasProfile:
    """The per-dataset bias profile produced by layer L2.

    This object is what the modelling and evaluation layers consume. It is
    deliberately a plain dataclass so it serialises to JSON without loss.
    """

    n_records: int
    label_entropy: float
    mutual_information: dict[str, float]
    per_allele_rho: dict[str, float] = field(default_factory=dict)
    n_alleles_excluded_single_sample: int = 0
    log_unit: str = "nats"

    @property
    def imbalanced_alleles(self) -> dict[str, float]:
        """Alleles whose positive-to-negative ratio exceeds 10:1 either way."""
        threshold = math.log2(10.0)
        return {
            a: r for a, r in self.per_allele_rho.items() if abs(r) > threshold
        }

    def bias_candidates(self, threshold_fraction: float = 0.25) -> dict[str, float]:
        """Variables carrying more than `threshold_fraction` of H(y).

        Reported as candidates, never as a verdict: whether a variable is a
        shortcut or a genuine determinant is a judgement about biology, not
        about arithmetic.
        """
        if self.label_entropy <= 0.0:
            return {}
        return {
            k: v
            for k, v in self.mutual_information.items()
            if v / self.label_entropy > threshold_fraction
        }

    def to_dict(self) -> dict:
        return {
            "n_records": self.n_records,
            "label_entropy": self.label_entropy,
            "log_unit": self.log_unit,
            "mutual_information": dict(self.mutual_information),
            "per_allele_rho": dict(self.per_allele_rho),
            "n_alleles_excluded_single_sample": self.n_alleles_excluded_single_sample,
            "n_alleles_imbalanced_beyond_10_to_1": len(self.imbalanced_alleles),
        }


def profile_dataset(
    labels: Sequence[int],
    variables: Mapping[str, Sequence[Hashable]],
    allele_key: str = "hla_allele",
) -> BiasProfile:
    """Compute the full bias profile for one corpus.

    `variables` maps a candidate bias variable name to a parallel sequence of
    its values, for example {"hla_allele": [...], "source_dataset": [...]}.
    The variable named by `allele_key` additionally drives the per-allele
    imbalance report.
    """
    ys = _as_list(labels)
    if not ys:
        raise ValueError("profile_dataset() requires at least one record")
    for name, vals in variables.items():
        if len(_as_list(vals)) != len(ys):
            raise ValueError(f"variable {name!r} is not parallel to labels")

    mi = {
        name: mutual_information(ys, vals) for name, vals in variables.items()
    }

    rho: dict[str, float] = {}
    excluded = 0
    if allele_key in variables:
        alleles = _as_list(variables[allele_key])
        rho_raw = per_allele_imbalance(alleles, ys)
        rho = {str(k): v for k, v in rho_raw.items()}
        excluded = len(set(alleles)) - len(rho_raw)

    return BiasProfile(
        n_records=len(ys),
        label_entropy=entropy(ys),
        mutual_information=mi,
        per_allele_rho=rho,
        n_alleles_excluded_single_sample=excluded,
    )
