"""Resampling statistics for the M6 harness (analysis plan section 5).

Pure numpy, no I/O, no knowledge of runs or reports. Everything here answers
one of two questions about a test set whose records are grouped by allele:

* how large is AUPRC(model) - AUPRC(control), and is it distinguishable from
  zero once the test alleles are resampled (the paired, allele-clustered
  bootstrap), and
* how stable is the shortcut-dependence kappa under the same resampling.

Why the cluster is the allele
-----------------------------
Predictions for peptides presented on one allele are correlated, and the test
set is a sample of alleles. A record-level bootstrap would treat them as
independent and understate the variance. One resample draws `n_alleles` alleles
with replacement; each drawn allele brings all of its test records. A resample
is represented by a vector of multiplicities (`counts`), which is exactly a
multinomial draw. The same `counts` matrix is reused for every comparison, so
every comparison is paired on the same resamples.

Ties (decision D14): seed-free, expectation over random orderings
------------------------------------------------------------------
`nbx.evaluate.metrics` breaks ties against the model (negatives first). That is
the conservative choice when the thing under test is a model. It is the
ANTI-conservative choice when the thing penalised is a CONTROL: an HLA-only
control assigns one score to every peptide of an allele, so it is a block of
ties, and it would be scored below chance relative to any scorer with
continuous output, including a model with no skill at all.

The definition used here is the EXPECTED average precision over uniformly random
orderings within each group of tied scores (`expected_average_precision`). It has
a closed form, involves no random draw and no seed, and for scores without ties it
equals `metrics.auprc` exactly. The ranking metrics (precision, recall and MCC at a
cut) are defined the same way.

Why not threshold-grouped average precision (`grouped_average_precision`, the
semantics of sklearn's average_precision_score, where equal scores form ONE
threshold)? It is seed-free too, and it is kept here as a reference, but it is NOT
the expectation over orderings: on tied blocks it is biased low relative to a
scorer with the same (absent) information and continuous output. Measured in the
suite, a zero-skill model with continuous scores beats an allele-constant control
under it far more often than alpha, which is exactly the hazard it was meant to
remove. The suite pins both facts.

A bootstrap resample repeats a record once per drawn copy of its allele; copies
of one record have equal scores and so form a tie group, handled by the same
rule. `weighted_auprc` evaluates that exactly, without expanding the resample.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from nbx.bias.diagnostics import spearman

__all__ = [
    "ALPHA",
    "MarginEstimate",
    "cluster_counts",
    "grouped_average_precision",
    "expected_average_precision",
    "expected_precision_at_k",
    "expected_recall_at_k",
    "expected_mcc_top_npos",
    "weighted_auprc",
    "bootstrap_p_value",
    "percentile_interval",
    "one_sided_bound",
    "holm_adjust",
    "margin_over_control",
    "kappa_resampled",
]

ALPHA = 0.05
"""Significance level of every one-sided test here (analysis plan section 5)."""

# Upper bound on the number of (resample, record) cells held in memory at once.
_CHUNK_CELLS = 1_500_000


def cluster_counts(n_clusters: int, n_boot: int, seed: int) -> np.ndarray:
    """Multiplicity of each cluster in each of `n_boot` resamples.

    Row b is a draw of `n_clusters` clusters with replacement (a multinomial with
    equal cell probabilities), so every row sums to `n_clusters`.
    """
    if n_clusters < 1 or n_boot < 1:
        raise ValueError("n_clusters and n_boot must be positive")
    rng = np.random.default_rng(seed)
    probs = np.full(n_clusters, 1.0 / n_clusters)
    return rng.multinomial(n_clusters, probs, size=n_boot).astype(np.int64)


# --------------------------------------------------------------------------
# Tie-aware point metrics (seed-free)
# --------------------------------------------------------------------------
def _grouped_cumulative(
    scores: np.ndarray, labels: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per tie group, best first: size, positives, and their cumulative sums."""
    s = np.asarray(scores, dtype=np.float64)
    y = np.asarray(labels, dtype=np.int64)
    if s.shape != y.shape or s.ndim != 1 or s.size == 0:
        raise ValueError("scores and labels must be equal-length, non-empty 1-d arrays")
    if np.any((y != 0) & (y != 1)):
        raise ValueError("labels must be binary 0/1")
    order = np.argsort(-s, kind="stable")
    s_sorted = s[order]
    y_sorted = y[order]
    boundaries = np.flatnonzero(np.diff(s_sorted) != 0) + 1
    starts = np.concatenate(([0], boundaries))
    ends = np.concatenate((boundaries, [s.size]))
    csum = np.concatenate(([0], np.cumsum(y_sorted)))
    size = ends - starts
    pos = csum[ends] - csum[starts]
    return size, pos, np.cumsum(size), np.cumsum(pos)


def grouped_average_precision(scores: np.ndarray, labels: np.ndarray) -> float:
    """Threshold-grouped average precision: tied scores share one threshold.

    sum over thresholds of (positives at the threshold / P) x (precision there),
    the semantics of sklearn's average_precision_score. REFERENCE ONLY: it is not
    the expectation over orderings and is not used by the harness (see the module
    docstring). Raises when there is no positive record rather than returning zero.
    """
    size, pos, cum_size, cum_pos = _grouped_cumulative(scores, labels)
    n_pos = int(pos.sum())
    if n_pos == 0:
        raise ValueError("AUPRC is undefined with no positive examples")
    return float(np.sum(pos * (cum_pos / cum_size)) / n_pos)


def _harmonic_table(max_total: int) -> np.ndarray:
    return np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, max_total + 1))))


def _expected_ap_numerator(
    size: np.ndarray, pos: np.ndarray, harmonic: np.ndarray
) -> np.ndarray:
    """Sum over tie groups of the expected sum of precisions at their positives.

    A group of m records holding p positives, preceded by N0 records of which T0 are
    positive, occupies ranks N0+1..N0+m in a uniformly random order. A positive at
    rank N0+r has (r-1)(p-1)/(m-1) other positives of the group ahead of it in
    expectation, and the precision is linear in that count, so
        E[sum of precisions] = sum_r (p/m) (T0 + 1 + (r-1)(p-1)/(m-1)) / (N0 + r),
    evaluated with harmonic numbers. Works on arrays whose last axis indexes groups.
    """
    n0 = np.cumsum(size, axis=-1) - size
    t0 = np.cumsum(pos, axis=-1) - pos
    s1 = harmonic[n0 + size] - harmonic[n0]               # sum 1/(N0+r)
    s2 = size - n0 * s1                                   # sum r/(N0+r)
    m = np.maximum(size, 1)
    spread = (pos - 1) / np.maximum(size - 1, 1)
    return np.sum((pos / m) * ((t0 + 1) * s1 + spread * (s2 - s1)), axis=-1)


def expected_average_precision(scores: np.ndarray, labels: np.ndarray) -> float:
    """Average precision expected over uniformly random orderings within ties.

    Seed-free and exact. Equals `metrics.auprc` when no two scores are equal.
    Raises when there is no positive record rather than returning zero.
    """
    size, pos, _, _ = _grouped_cumulative(scores, labels)
    n_pos = int(pos.sum())
    if n_pos == 0:
        raise ValueError("AUPRC is undefined with no positive examples")
    harmonic = _harmonic_table(int(size.sum()))
    return float(_expected_ap_numerator(size, pos, harmonic) / n_pos)


def _expected_top_positives(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Expected number of positives in the top k over random orderings inside ties."""
    size, pos, cum_size, cum_pos = _grouped_cumulative(scores, labels)
    n = int(size.sum())
    if k <= 0:
        raise ValueError("k must be positive")
    if k > n:
        raise ValueError(f"k={k} exceeds the {n} available candidates")
    g = int(np.searchsorted(cum_size, k))          # the group holding rank k
    before_size = int(cum_size[g - 1]) if g > 0 else 0
    before_pos = int(cum_pos[g - 1]) if g > 0 else 0
    slots = k - before_size
    return before_pos + slots * float(pos[g]) / float(size[g])


def expected_precision_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """PPV@k, expected over random orderings within tied groups. Raises if k > n."""
    return _expected_top_positives(scores, labels, k) / k


def expected_recall_at_k(scores: np.ndarray, labels: np.ndarray, k: int) -> float:
    """Recall@k with k clamped to n (as metrics.recall_at_k), tie-expected."""
    y = np.asarray(labels, dtype=np.int64)
    n_pos = int(y.sum())
    if n_pos == 0:
        raise ValueError("recall is undefined with no positive examples")
    return _expected_top_positives(scores, labels, min(k, y.size)) / n_pos


def expected_mcc_top_npos(scores: np.ndarray, labels: np.ndarray) -> float:
    """MCC when the n_pos highest scores are predicted positive.

    Computed from the EXPECTED confusion counts under random ordering inside ties
    (a plug-in of the expected counts, not the expectation of MCC). Raises when one
    class is absent: MCC is then not a measurement.
    """
    y = np.asarray(labels, dtype=np.int64)
    n, n_pos = int(y.size), int(y.sum())
    if n_pos == 0 or n_pos == n:
        raise ValueError("MCC is undefined when only one class is present")
    tp = _expected_top_positives(scores, labels, n_pos)
    fp = n_pos - tp
    fn = n_pos - tp
    tn = (n - n_pos) - fp
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denom == 0.0:
        raise ValueError("MCC denominator vanishes")
    return float((tp * tn - fp * fn) / denom)


# --------------------------------------------------------------------------
# The bootstrap
# --------------------------------------------------------------------------
def weighted_auprc(
    scores: np.ndarray,
    labels: np.ndarray,
    allele_idx: np.ndarray,
    counts: np.ndarray,
) -> np.ndarray:
    """Expected average precision of every resample in `counts` (length n_boot).

    Equal to `expected_average_precision` on the expanded resample (a record
    repeated as many times as its allele was drawn). A resample with no positive
    record is undefined and returns NaN, never zero.
    """
    s = np.asarray(scores, dtype=np.float64)
    y_all = np.asarray(labels, dtype=np.int64)
    a_all = np.asarray(allele_idx, dtype=np.int64)
    n = s.size
    n_boot = counts.shape[0]
    order = np.argsort(-s, kind="stable")
    s_sorted = s[order]
    y = y_all[order]
    a = a_all[order]
    starts = np.concatenate(([0], np.flatnonzero(np.diff(s_sorted) != 0) + 1))

    sizes = np.bincount(a_all, minlength=counts.shape[1])
    harmonic = _harmonic_table(int(np.max(counts @ sizes)))

    out = np.empty(n_boot, dtype=np.float64)
    step = max(1, _CHUNK_CELLS // max(n, 1))
    for lo in range(0, n_boot, step):
        hi = min(n_boot, lo + step)
        w = counts[lo:hi][:, a]                         # copies of each record
        wp = w * y                                      # copies of each positive
        size = np.add.reduceat(w, starts, axis=1)       # per tie group
        pos = np.add.reduceat(wp, starts, axis=1)
        n_pos = wp.sum(axis=1)
        total = _expected_ap_numerator(size, pos, harmonic)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[lo:hi] = np.where(n_pos > 0, total / np.maximum(n_pos, 1), np.nan)
    return out


def bootstrap_p_value(deltas: np.ndarray) -> tuple[float | None, int]:
    """One-sided p = (1 + #{b : delta_b <= 0}) / (B + 1) over defined resamples.

    Returns (None, 0) when no resample is defined. NaN resamples are excluded and
    counted out of B; the second element is the number that were used. The opposite
    tail is `bootstrap_p_value(-deltas)`.
    """
    d = np.asarray(deltas, dtype=np.float64)
    d = d[~np.isnan(d)]
    if d.size == 0:
        return None, 0
    return (1.0 + float(np.sum(d <= 0.0))) / (d.size + 1.0), int(d.size)


def percentile_interval(
    values: np.ndarray, level: float = 0.90
) -> tuple[float | None, float | None, int]:
    """Two-sided percentile interval over the defined resamples.

    The default level 0.90 (5th and 95th percentiles) is the interval that matches
    a one-sided test at alpha = 0.05.
    """
    v = np.asarray(values, dtype=np.float64)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return None, None, 0
    tail = (1.0 - level) / 2.0 * 100.0
    lo, hi = np.percentile(v, [tail, 100.0 - tail])
    return float(lo), float(hi), int(v.size)


def one_sided_bound(values: np.ndarray, side: str, level: float = 0.95) -> float | None:
    """One-sided percentile bound: the 5th percentile ("lower") or 95th ("upper")."""
    v = np.asarray(values, dtype=np.float64)
    v = v[~np.isnan(v)]
    if v.size == 0:
        return None
    if side == "lower":
        return float(np.percentile(v, (1.0 - level) * 100.0))
    if side == "upper":
        return float(np.percentile(v, level * 100.0))
    raise ValueError("side must be 'lower' or 'upper'")


def holm_adjust(p_values: Sequence[float]) -> list[float]:
    """Holm step-down adjusted p-values, in the order given."""
    m = len(p_values)
    order = sorted(range(m), key=lambda i: p_values[i])
    adjusted = [0.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (m - rank) * p_values[i]))
        adjusted[i] = running
    return adjusted


@dataclass(frozen=True)
class MarginEstimate:
    """AUPRC(model) - AUPRC(control), seed-averaged, with its bootstrap test.

    `p_value` is one-sided against delta <= 0. `ci90_low`/`ci90_high` are the 5th
    and 95th percentiles of the same bootstrap distribution, the two-sided interval
    that matches that test; `one_sided_lower95` is its lower end by definition.
    `p_value_worse` is the opposite tail, and `significantly_worse` means the point
    delta is negative and that tail is below alpha.
    """

    delta: float
    per_seed_mean: float
    per_seed_sd: float | None
    p_value: float | None
    p_value_worse: float | None
    significantly_worse: bool
    ci90_low: float | None
    ci90_high: float | None
    one_sided_lower95: float | None
    n_valid_resamples: int
    n_boot: int
    n_seeds: int


def margin_over_control(
    model_scores: Sequence[np.ndarray],
    control_scores: Sequence[np.ndarray],
    labels: np.ndarray,
    allele_idx: np.ndarray,
    counts: np.ndarray,
) -> MarginEstimate:
    """Paired, allele-clustered bootstrap of the seed-averaged AUPRC margin.

    `model_scores[i]` and `control_scores[i]` are the predictions of seed i over
    the SAME test records, in the same order as `labels` and `allele_idx`. The
    statistic is mean over seeds of AUPRC(model) - AUPRC(control). Each resample
    draws test alleles (rows of `counts`) and evaluates every seed on it, so the
    pairing is by seed and by resample. The test is one-sided against delta <= 0;
    the opposite tail answers "significantly worse".
    """
    if len(model_scores) != len(control_scores) or not model_scores:
        raise ValueError("model and control must supply the same, non-zero seeds")
    y = np.asarray(labels, dtype=np.int64)
    a = np.asarray(allele_idx, dtype=np.int64)
    if counts.shape[1] != int(a.max()) + 1:
        raise ValueError("counts must have one column per allele cluster")

    models = [np.asarray(m, dtype=np.float64) for m in model_scores]
    controls = [np.asarray(c, dtype=np.float64) for c in control_scores]

    per_seed = np.array(
        [
            expected_average_precision(m, y) - expected_average_precision(c, y)
            for m, c in zip(models, controls)
        ]
    )
    acc = np.zeros(counts.shape[0], dtype=np.float64)
    for m, c in zip(models, controls):
        acc += weighted_auprc(m, y, a, counts) - weighted_auprc(c, y, a, counts)
    boot = acc / len(models)

    p, n_valid = bootstrap_p_value(boot)
    p_worse, _ = bootstrap_p_value(-boot)
    lo, hi, _ = percentile_interval(boot, 0.90)
    delta = float(per_seed.mean())
    sd = float(np.std(per_seed, ddof=1)) if per_seed.size > 1 else None
    return MarginEstimate(
        delta=delta,
        per_seed_mean=delta,
        per_seed_sd=sd,
        p_value=p,
        p_value_worse=p_worse,
        significantly_worse=bool(delta < 0.0 and p_worse is not None and p_worse < ALPHA),
        ci90_low=lo,
        ci90_high=hi,
        one_sided_lower95=one_sided_bound(boot, "lower"),
        n_valid_resamples=n_valid,
        n_boot=int(counts.shape[0]),
        n_seeds=len(models),
    )


def kappa_resampled(
    rho: np.ndarray,
    mean_scores_by_seed: Sequence[np.ndarray],
    counts: np.ndarray,
) -> np.ndarray:
    """Seed-averaged kappa of every allele resample (NaN where undefined).

    `rho[j]` is the TRAINING imbalance of allele j and `mean_scores_by_seed[s][j]`
    its mean predicted test score under seed s, for the alleles shared by train
    and test. `counts[b, j]` is how often allele j is drawn in resample b. kappa is
    the project's own `spearman` of the two vectors, so the definition is the one
    `shortcut_dependence` uses; with unit counts the two agree exactly.
    """
    rho_v = np.asarray(rho, dtype=np.float64)
    means = [np.asarray(m, dtype=np.float64) for m in mean_scores_by_seed]
    out = np.full(counts.shape[0], np.nan)
    positions = np.arange(rho_v.size)
    for b in range(counts.shape[0]):
        idx = np.repeat(positions, counts[b])
        if idx.size < 3:
            continue
        xs = rho_v[idx].tolist()
        try:
            out[b] = float(np.mean([spearman(xs, m[idx].tolist()) for m in means]))
        except ValueError:
            continue                                  # undefined: counted out
    return out
