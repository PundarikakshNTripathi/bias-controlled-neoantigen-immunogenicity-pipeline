"""Evaluation metrics (Section 7.7.1).

AUPRC is the headline discrimination metric, in preference to AUROC, because
the class imbalance in this problem is extreme and the ROC curve is
correspondingly optimistic about it.

Alongside it the TESLA-derived ranking metrics are reported, so results are
comparable with the consortium literature:

    fraction ranked          recall within the top 100 candidates
    top-20 immunogenic       immunogenic proportion of the top 20

All functions take scores and binary labels and return a float, or raise.
None of them silently return 0.0 for a degenerate input: a metric that
cannot be computed is an error, not a zero.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

__all__ = [
    "auprc",
    "auroc",
    "precision_at_k",
    "recall_at_k",
    "fraction_ranked",
    "top_k_immunogenic_fraction",
    "matthews_corrcoef",
    "allele_stratified",
]


def _check(scores: Sequence[float], labels: Sequence[int]) -> tuple[list, list]:
    s = [float(x) for x in scores]
    y = [int(x) for x in labels]
    if len(s) != len(y):
        raise ValueError("scores and labels must be parallel sequences")
    if not s:
        raise ValueError("metric requires at least one observation")
    if any(v not in (0, 1) for v in y):
        raise ValueError("labels must be binary 0/1")
    return s, y


def _descending(scores: Sequence[float], labels: Sequence[int]) -> list[int]:
    """Labels ordered by descending score, ties broken pessimistically.

    A pessimistic tie break puts negatives ahead of positives at the same
    score. Optimistic tie breaking is a quiet way to inflate ranking metrics
    on a model that assigns many identical scores.
    """
    order = sorted(range(len(scores)), key=lambda i: (-scores[i], labels[i]))
    return [labels[i] for i in order]


def auprc(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Area under the precision-recall curve, by exact step summation.

    Equivalent to sklearn's average_precision_score: the sum over thresholds
    of (recall increment) x (precision at that point). Verified against it in
    the test suite.
    """
    s, y = _check(scores, labels)
    n_pos = sum(y)
    if n_pos == 0:
        raise ValueError("AUPRC is undefined with no positive examples")

    ordered = _descending(s, y)
    tp = 0
    total = 0.0
    for i, label in enumerate(ordered, start=1):
        if label == 1:
            tp += 1
            total += tp / i          # precision at this recall step
    return total / n_pos


def auroc(scores: Sequence[float], labels: Sequence[int]) -> float:
    """Area under the ROC curve, via the rank-sum identity, ties averaged.

    Reported as secondary context only. The documents are explicit that
    AUROC is optimistic under this problem's class imbalance.
    """
    s, y = _check(scores, labels)
    n_pos, n_neg = sum(y), len(y) - sum(y)
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC is undefined when one class is absent")

    order = sorted(range(len(s)), key=lambda i: s[i])
    ranks = [0.0] * len(s)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and s[order[j + 1]] == s[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1

    rank_sum = sum(r for r, label in zip(ranks, y) if label == 1)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def precision_at_k(scores: Sequence[float], labels: Sequence[int], k: int) -> float:
    """PPV@k = TP@k / k. The operative quantity in a clinical setting.

    A manufactured therapy can include only a handful of peptides, so the
    ranking's behaviour in its top k is what determines therapeutic value,
    not its average behaviour over the whole candidate list.
    """
    s, y = _check(scores, labels)
    if k <= 0:
        raise ValueError("k must be positive")
    if k > len(s):
        raise ValueError(f"k={k} exceeds the {len(s)} available candidates")
    return sum(_descending(s, y)[:k]) / k


def recall_at_k(scores: Sequence[float], labels: Sequence[int], k: int) -> float:
    """Recall@k = TP@k / P."""
    s, y = _check(scores, labels)
    n_pos = sum(y)
    if n_pos == 0:
        raise ValueError("recall is undefined with no positive examples")
    if k <= 0:
        raise ValueError("k must be positive")
    k = min(k, len(s))
    return sum(_descending(s, y)[:k]) / n_pos


def fraction_ranked(scores: Sequence[float], labels: Sequence[int]) -> float:
    """TESLA metric: recall within the top 100 candidates."""
    return recall_at_k(scores, labels, 100)


def top_k_immunogenic_fraction(
    scores: Sequence[float], labels: Sequence[int], k: int = 20
) -> float:
    """TESLA metric: the immunogenic proportion of the top k ranked candidates."""
    return precision_at_k(scores, labels, min(k, len(scores)))


def matthews_corrcoef(predictions: Sequence[int], labels: Sequence[int]) -> float:
    """Matthews correlation coefficient, for balanced-split comparisons.

    Takes hard predictions, not scores. Returns 0.0 when the denominator
    vanishes, which is the conventional definition rather than a placeholder:
    a predictor that outputs one class everywhere has no correlation with the
    truth, and that is a measurement.
    """
    p, y = _check([float(v) for v in predictions], labels)
    pred = [int(v) for v in p]
    if any(v not in (0, 1) for v in pred):
        raise ValueError("predictions must be binary 0/1")

    tp = sum(1 for a, b in zip(pred, y) if a == 1 and b == 1)
    tn = sum(1 for a, b in zip(pred, y) if a == 0 and b == 0)
    fp = sum(1 for a, b in zip(pred, y) if a == 1 and b == 0)
    fn = sum(1 for a, b in zip(pred, y) if a == 0 and b == 1)

    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0.0:
        return 0.0
    return (tp * tn - fp * fn) / denom


def allele_stratified(
    metric,
    scores: Sequence[float],
    labels: Sequence[int],
    alleles: Sequence[str],
    weight_by_support: bool = True,
) -> tuple[dict[str, float], float]:
    """Compute a metric per allele and aggregate with explicit weighting (R3).

    Returns (per_allele, aggregate). Alleles for which the metric is
    undefined are omitted from both, and their omission is visible because
    they are absent from the per-allele mapping.

    Pooled metrics on an imbalanced allele distribution conflate coverage
    with accuracy, which is why this function exists and why the aggregate is
    never computed by simply pooling all records together.
    """
    s, y = _check(scores, labels)
    if len(alleles) != len(s):
        raise ValueError("alleles must be parallel to scores and labels")

    groups: dict[str, tuple[list[float], list[int]]] = {}
    for score, label, allele in zip(s, y, alleles):
        a = str(allele)
        groups.setdefault(a, ([], []))
        groups[a][0].append(score)
        groups[a][1].append(label)

    per_allele: dict[str, float] = {}
    supports: dict[str, int] = {}
    for allele, (gs, gy) in sorted(groups.items()):
        try:
            per_allele[allele] = metric(gs, gy)
            supports[allele] = len(gy)
        except ValueError:
            continue          # undefined for this allele; omitted, not zeroed

    if not per_allele:
        raise ValueError("the metric was undefined for every allele")

    if weight_by_support:
        total = sum(supports[a] for a in per_allele)
        aggregate = sum(per_allele[a] * supports[a] for a in per_allele) / total
    else:
        aggregate = sum(per_allele.values()) / len(per_allele)

    return per_allele, aggregate
