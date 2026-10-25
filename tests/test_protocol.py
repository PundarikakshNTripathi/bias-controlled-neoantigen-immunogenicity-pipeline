"""Validity tests of the M6 evaluation harness itself.

These do not test a model. They test the instrument: that it can DETECT a real
margin, that it REFUSES a shortcut, that it is calibrated under the null, that
it declines to produce a headline number when its preconditions fail, and that
it is deterministic. Every dataset here is synthetic and is built to a stated
generative rule; no number below is a finding about immunogenicity.

Requirements pinned (analysis_plan.md sections 4, 5, 6, 9; SRS FR-60..FR-68;
decisions D4, D14, D15, D16 of the master plan):
  FR-62  matched controls          -> refusal and structural tests
  FR-63  no metrics w/o control    -> PENDING and refusal tests
  FR-68  mean +/- SD, no best seed -> determinism / reporting tests
  D4     burden of proof on model  -> gate scenarios
  D14    seed-free tie handling    -> stats tests
  D15    coherent test and interval-> margin block tests
  D16    pins and structural checks-> seeds, digest, structure, CLI tests
"""

from __future__ import annotations

import ast
import dataclasses
import inspect
import itertools
import json
import re
import zlib
from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import average_precision_score
from typer.testing import CliRunner

from nbx.bias.diagnostics import per_allele_imbalance, shortcut_dependence
from nbx.cli import app
from nbx.evaluate import protocol, stats
from nbx.evaluate.gate import UNMEASURED, GateStatus
from nbx.evaluate.metrics import auprc, matthews_corrcoef, precision_at_k, recall_at_k
from nbx.evaluate.protocol import (
    DEFAULT_BOOTSTRAP_SEED,
    DEFAULT_RANKING_K,
    EXPECTED_SEEDS,
    MIN_SEEDS,
    PRODUCTION_N_BOOT,
    ProtocolRefusal,
    RunConfig,
    SplitRecord,
    VariantRun,
    evaluate_runs,
    report_to_json,
    run_protocol,
    write_run_dir,
)
from nbx.splits.freeze import SplitCriterion, build_split

SEEDS = tuple(sorted(EXPECTED_SEEDS))
CFG = RunConfig(
    epochs=10, n_parameters=1000, early_stopping_metric="val_auprc", tuning_trials=5
)
N_BOOT = 300          # small on purpose: the suite must stay fast
K = 20                # ranking k, small because the synthetic test set is small
NO_METRIC_TEXT = r"kappa|\d\.\d{3}"


# ==========================================================================
# Synthetic worlds
# ==========================================================================
@dataclasses.dataclass
class World:
    record: SplitRecord
    test_idx: np.ndarray            # sorted corpus indices of the test records
    test_alleles: list[str]
    test_labels: np.ndarray
    x: np.ndarray                   # peptide feature, indexed by corpus index
    train_rho: dict
    train_rate: dict

    @property
    def n(self) -> int:
        return len(self.test_idx)


def make_world(kind: str, n_alleles: int = 16, per_allele: int = 200,
               seed: int = 11, extra_small_allele: bool = False) -> World:
    """kind='shortcut': label depends ONLY on the allele (rate p_a); x is noise.
    kind='signal': label still has an allele-specific rate (so TRAINING is
    imbalanced), and the peptide feature x is informative about the label with the
    same class-conditional distribution in every allele. Both are split with the
    real peptide-grouped splitter, so the test set is allele-balanced by criterion."""
    rng = np.random.default_rng(seed)
    p_a = np.linspace(0.12, 0.88, n_alleles)
    rng.shuffle(p_a)
    alleles, labels, xs = [], [], []
    for a in range(n_alleles):
        y = (rng.random(per_allele) < p_a[a]).astype(int)
        if kind == "shortcut":
            x = rng.normal(size=per_allele)              # carries no information
        else:
            # The peptide feature depends on the label alone, identically for every
            # allele, so the signal is a property of the peptide and not the allele.
            x = 0.75 * (2 * y - 1) + rng.normal(size=per_allele)
        alleles += [f"A{a:02d}"] * per_allele
        labels += y.tolist()
        xs += x.tolist()
    if extra_small_allele:                # below min_samples_per_allele
        for _ in range(10):
            alleles.append("SMALL")
            labels.append(int(rng.random() < 0.5))
            xs.append(float(rng.normal()))
    peptides = [f"PEP{i:05d}" for i in range(len(alleles))]      # one per record
    split = build_split(alleles, labels, SplitCriterion(seed=3), groups=peptides)
    record = SplitRecord.from_frozen(split, alleles, labels, peptides)
    test_idx = np.array(split.test_idx)
    train_alleles = [alleles[i] for i in split.train_idx]
    train_labels = [labels[i] for i in split.train_idx]
    rho = per_allele_imbalance(train_alleles, train_labels)
    rate: dict = {}
    for a in set(train_alleles):
        ys = [y for al, y in zip(train_alleles, train_labels) if al == a]
        rate[a] = sum(ys) / len(ys)
    return World(
        record=record,
        test_idx=test_idx,
        test_alleles=[alleles[i] for i in test_idx],
        test_labels=np.array([labels[i] for i in test_idx]),
        x=np.array(xs),
        train_rho=rho,
        train_rate=rate,
    )


def run_of(world: World, variant: str, role: str, seed: int, scores: np.ndarray,
           matched_to: str | None = None, config: RunConfig = CFG) -> VariantRun:
    assert len(scores) == world.n
    return VariantRun(
        variant=variant, role=role, seed=seed,
        scores={int(i): float(s) for i, s in zip(world.test_idx, scores)},
        config=config, matched_to=matched_to,
    )


def rho_of(world: World) -> np.ndarray:
    return np.array([world.train_rho[a] for a in world.test_alleles])


def rate_of(world: World) -> np.ndarray:
    return np.array([world.train_rate[a] for a in world.test_alleles])


def name_key(name: str) -> int:
    return zlib.crc32(name.encode("utf-8")) % 1000


def noise(rng, n, sd=1.0):
    return rng.normal(scale=sd, size=n)


def assemble(world: World, model_fn, hla_fn, pep_fn, base_fn=None, base_hla_fn=None,
             base_pep_fn=None, seeds=SEEDS, model="model") -> list[VariantRun]:
    """model + matched controls (+ a baseline with its own matched controls).
    Each *_fn(seed, rng) -> scores over the world's test records."""
    runs: list[VariantRun] = []

    def add(name, role, fn, matched_to=None):
        for s in seeds:
            rng = np.random.default_rng([s, name_key(name)])
            runs.append(run_of(world, name, role, s, fn(s, rng), matched_to))

    add(model, "model", model_fn)
    add(f"{model}__hla", "hla_only", hla_fn, model)
    add(f"{model}__pep", "peptide_only", pep_fn, model)
    if base_fn is not None:
        add("base", "baseline", base_fn)
        add("base__hla", "hla_only", base_hla_fn or hla_fn, "base")
        add("base__pep", "peptide_only", base_pep_fn or pep_fn, "base")
    return runs


def evaluate(world: World, runs, record: SplitRecord | None = None, **kw):
    params = dict(n_boot=N_BOOT, bootstrap_seed=7, ranking_k=K)
    params.update(kw)
    rec = record or world.record
    return evaluate_runs(rec, runs, split_hash=rec.content_hash, **params)


@pytest.fixture(scope="module")
def shortcut_world() -> World:
    return make_world("shortcut")


@pytest.fixture(scope="module")
def signal_world() -> World:
    return make_world("signal")


# ==========================================================================
# Statistics building blocks (D14: seed-free, threshold-grouped AUPRC)
# ==========================================================================
def _random_case(rng, n, tied):
    labels = (rng.random(n) < rng.uniform(0.2, 0.7)).astype(int)
    if labels.sum() in (0,):
        labels[0] = 1
    scores = rng.normal(size=n)
    if tied:
        scores = np.round(scores, rng.integers(0, 2))      # many repeated values
    return scores, labels


def test_grouped_average_precision_matches_sklearn_on_tied_and_untied_cases():
    rng = np.random.default_rng(0)
    worst = 0.0
    for case in range(400):
        scores, labels = _random_case(rng, int(rng.integers(5, 80)), tied=case % 2 == 0)
        got = stats.grouped_average_precision(scores, labels)
        want = average_precision_score(labels, scores)
        worst = max(worst, abs(got - want))
    assert worst < 1e-12, worst


def test_untied_scores_give_exactly_metrics_auprc():
    rng = np.random.default_rng(1)
    for _ in range(100):
        scores, labels = _random_case(rng, int(rng.integers(5, 80)), tied=False)
        assert len(set(scores.tolist())) == len(scores)
        got = stats.expected_average_precision(scores, labels)
        assert got == pytest.approx(auprc(scores.tolist(), labels.tolist()), abs=1e-12)
        # and the threshold-grouped reference agrees with both when nothing ties
        assert stats.grouped_average_precision(scores, labels) == pytest.approx(
            got, abs=1e-12)


def _expected_ap_by_enumeration(scores, labels):
    """Average precision averaged over every ordering consistent with the scores."""
    n = len(scores)
    total, count = 0.0, 0
    for perm in itertools.permutations(range(n)):
        ordered = [scores[i] for i in perm]
        if any(ordered[i] < ordered[i + 1] for i in range(n - 1)):
            continue                                   # not a valid descending order
        tp, ap = 0, 0.0
        for rank, i in enumerate(perm, start=1):
            if labels[i] == 1:
                tp += 1
                ap += tp / rank
        total += ap / sum(labels)
        count += 1
    return total / count


def test_expected_average_precision_is_the_exact_mean_over_tie_orderings():
    rng = np.random.default_rng(10)
    worst = 0.0
    for _ in range(40):
        scores = np.round(rng.normal(size=7), 0)           # heavy ties
        labels = (rng.random(7) < 0.5).astype(int)
        labels[0] = 1
        want = _expected_ap_by_enumeration(scores.tolist(), labels.tolist())
        worst = max(worst, abs(stats.expected_average_precision(scores, labels) - want))
    assert worst < 1e-12, worst


def test_tied_block_is_not_ordered_by_position_or_label():
    labels = np.array([1, 0, 1, 0, 1, 0])
    same = np.zeros(6)
    want = _expected_ap_by_enumeration(same.tolist(), labels.tolist())
    assert stats.expected_average_precision(same, labels) == pytest.approx(want)
    assert stats.expected_average_precision(same, labels) == (
        stats.expected_average_precision(same[::-1], labels[::-1])
    )
    # the pessimistic convention of metrics.auprc is strictly lower
    assert auprc(same.tolist(), labels.tolist()) < want


def test_threshold_grouped_ap_is_biased_against_a_constant_score_control():
    """Why the harness uses the expectation over orderings and not sklearn's
    threshold-grouped definition. Zero-skill scorers with continuous output have
    mean AP at the chance level; an allele-constant control scores exactly its
    expectation under the harness definition but below it under the grouped one."""
    n_alleles, per = 16, 30
    labels = np.tile([1, 0], n_alleles * per // 2)
    lookup = np.repeat(np.linspace(0.0, 1.0, n_alleles), per)
    rng = np.random.default_rng(11)
    random_ap = np.array([
        stats.expected_average_precision(rng.normal(size=labels.size), labels)
        for _ in range(300)
    ])
    se = random_ap.std(ddof=1) / np.sqrt(random_ap.size)
    expected_control = stats.expected_average_precision(lookup, labels)
    grouped_control = stats.grouped_average_precision(lookup, labels)
    assert abs(random_ap.mean() - expected_control) < 4 * se
    assert random_ap.mean() - grouped_control > 4 * se
    assert grouped_control == pytest.approx(average_precision_score(labels, lookup))


def test_weighted_auprc_equals_expected_ap_of_the_expanded_resample_with_ties():
    rng = np.random.default_rng(2)
    n_alleles, per = 7, 15
    allele_idx = np.repeat(np.arange(n_alleles), per)
    labels = (rng.random(n_alleles * per) < 0.4).astype(int)
    for tied in (False, True):
        scores = rng.normal(size=n_alleles * per)
        if tied:
            scores = np.round(scores, 0)
        counts = rng.multinomial(n_alleles, np.full(n_alleles, 1 / n_alleles), size=25)
        got = stats.weighted_auprc(scores, labels, allele_idx, counts)
        for b in range(counts.shape[0]):
            expanded = np.repeat(np.arange(len(labels)), counts[b, allele_idx])
            if labels[expanded].sum() == 0:
                assert np.isnan(got[b])
                continue
            want = stats.expected_average_precision(
                scores[expanded], labels[expanded])
            assert got[b] == pytest.approx(want, abs=1e-12)


def test_weighted_auprc_with_unit_counts_is_the_plain_auprc():
    rng = np.random.default_rng(3)
    allele_idx = np.repeat(np.arange(4), 25)
    labels = (rng.random(100) < 0.5).astype(int)
    scores = rng.normal(size=100)
    ones = np.ones((1, 4), dtype=np.int64)
    got = stats.weighted_auprc(scores, labels, allele_idx, ones)[0]
    assert got == pytest.approx(auprc(scores.tolist(), labels.tolist()), abs=1e-12)


def test_there_is_no_tie_seed_anywhere_in_the_statistics_module():
    assert not hasattr(stats, "break_ties")
    for name in ("margin_over_control", "weighted_auprc",
                 "expected_average_precision"):
        params = inspect.signature(getattr(stats, name)).parameters
        assert not any("tie" in p for p in params)


def test_bootstrap_seed_changes_only_the_resampling_not_the_point_estimates():
    rng = np.random.default_rng(4)
    n_alleles, per = 12, 20
    labels = np.tile([1, 0], n_alleles * per // 2)
    allele = np.repeat(np.arange(n_alleles), per)
    lookup = np.repeat(np.linspace(0, 1, n_alleles), per)          # heavy ties
    model = [rng.normal(size=labels.size) for _ in SEEDS]
    ctrl = [lookup for _ in SEEDS]
    a = stats.margin_over_control(model, ctrl, labels, allele,
                                  stats.cluster_counts(n_alleles, 100, 1))
    b = stats.margin_over_control(model, ctrl, labels, allele,
                                  stats.cluster_counts(n_alleles, 100, 2))
    assert a.delta == b.delta and a.per_seed_sd == b.per_seed_sd
    assert (a.p_value, a.ci90_low) != (b.p_value, b.ci90_low)      # resampling differs


def _expected_topk_positives_by_enumeration(scores, labels, k):
    """Average top-k positives over every ordering consistent with the scores."""
    n = len(scores)
    totals, count = 0, 0
    for perm in itertools.permutations(range(n)):
        ordered = [scores[i] for i in perm]
        if any(ordered[i] < ordered[i + 1] for i in range(n - 1)):
            continue                                   # not a valid descending order
        totals += sum(labels[i] for i in perm[:k])
        count += 1
    return totals / count


def test_expected_precision_and_recall_average_over_orderings_within_ties():
    rng = np.random.default_rng(5)
    for _ in range(25):
        n = 6
        scores = np.round(rng.normal(size=n), 0)
        labels = (rng.random(n) < 0.5).astype(int)
        if labels.sum() == 0:
            labels[0] = 1
        for k in (1, 2, 3, 5):
            exp_tp = _expected_topk_positives_by_enumeration(
                scores.tolist(), labels.tolist(), k
            )
            assert stats.expected_precision_at_k(scores, labels, k) == pytest.approx(
                exp_tp / k, abs=1e-12)
            assert stats.expected_recall_at_k(scores, labels, k) == pytest.approx(
                exp_tp / labels.sum(), abs=1e-12)


def test_untied_ranking_metrics_equal_metrics_py():
    rng = np.random.default_rng(6)
    for _ in range(60):
        scores, labels = _random_case(rng, 40, tied=False)
        for k in (1, 7, 20):
            assert stats.expected_precision_at_k(scores, labels, k) == pytest.approx(
                precision_at_k(scores.tolist(), labels.tolist(), k), abs=1e-12)
            assert stats.expected_recall_at_k(scores, labels, k) == pytest.approx(
                recall_at_k(scores.tolist(), labels.tolist(), k), abs=1e-12)
        n_pos = int(labels.sum())
        pred = np.zeros(40, dtype=int)
        pred[np.argsort(-scores, kind="stable")[:n_pos]] = 1
        if 0 < n_pos < 40:
            assert stats.expected_mcc_top_npos(scores, labels) == pytest.approx(
                matthews_corrcoef(pred.tolist(), labels.tolist()), abs=1e-12)


def test_ranking_metrics_refuse_degenerate_inputs_instead_of_returning_zero():
    with pytest.raises(ValueError):
        stats.expected_precision_at_k(np.arange(3.0), np.array([1, 0, 1]), 5)
    with pytest.raises(ValueError):
        stats.expected_mcc_top_npos(np.arange(3.0), np.array([1, 1, 1]))
    with pytest.raises(ValueError):
        stats.grouped_average_precision(np.arange(3.0), np.array([0, 0, 0]))


def test_p_value_formula_and_exclusion_of_undefined_resamples():
    deltas = np.array([-0.1, 0.0, 0.2, 0.3, np.nan])
    p, n_valid = stats.bootstrap_p_value(deltas)
    assert n_valid == 4
    assert p == pytest.approx((1 + 2) / (4 + 1))        # -0.1 and 0.0 are <= 0
    p_worse, _ = stats.bootstrap_p_value(-deltas)       # the opposite tail
    assert p_worse == pytest.approx((1 + 3) / (4 + 1))    # -0.0, -0.2, -0.3


def test_interval_helpers_are_the_matching_90_percent_pair():
    x = np.random.default_rng(7).normal(size=500)
    lo, hi, n = stats.percentile_interval(x, level=0.90)
    assert n == 500
    assert lo == pytest.approx(np.percentile(x, 5))
    assert hi == pytest.approx(np.percentile(x, 95))
    assert stats.one_sided_bound(x, "lower") == pytest.approx(np.percentile(x, 5))
    assert stats.one_sided_bound(x, "upper") == pytest.approx(np.percentile(x, 95))


def test_holm_adjustment():
    adj = stats.holm_adjust([0.01, 0.04, 0.03])
    assert adj == pytest.approx([0.03, 0.06, 0.06])
    assert stats.holm_adjust([0.5]) == pytest.approx([0.5])
    assert stats.holm_adjust([0.4, 0.9]) == pytest.approx([0.8, 0.9])
    assert max(stats.holm_adjust([0.9, 0.99])) <= 1.0


def test_cluster_counts_are_resamples_with_replacement_and_seeded():
    c = stats.cluster_counts(n_clusters=9, n_boot=50, seed=5)
    assert c.shape == (50, 9)
    assert (c.sum(axis=1) == 9).all() and (c >= 0).all()
    assert np.array_equal(c, stats.cluster_counts(9, 50, 5))
    assert not np.array_equal(c, stats.cluster_counts(9, 50, 6))


def test_kappa_bootstrap_with_unit_counts_equals_shortcut_dependence():
    rng = np.random.default_rng(8)
    alleles = [f"A{i}" for i in range(12)]
    rho = {a: float(rng.normal()) for a in alleles}
    rows_a = [a for a in alleles for _ in range(5)]
    rows_s = rng.normal(size=len(rows_a))
    want = shortcut_dependence(rho, rows_a, rows_s.tolist())
    mean_by_allele = np.array(
        [rows_s[[i for i, a in enumerate(rows_a) if a == al]].mean() for al in alleles]
    )
    got = stats.kappa_resampled(
        np.array([rho[a] for a in alleles]), [mean_by_allele],
        np.ones((1, len(alleles)), dtype=np.int64),
    )
    assert got[0] == pytest.approx(want, abs=1e-12)


# ==========================================================================
# Scenarios
# ==========================================================================
def _shortcut_fns(world: World):
    shortcut = rate_of(world)

    def model(s, rng):                      # exploits the shortcut, no noise
        return shortcut.copy()

    def hla(s, rng):                        # a pure per-allele lookup
        return shortcut.copy()

    def pep(s, rng):
        return noise(rng, world.n)

    def base(s, rng):                       # an undebiased reference, same shortcut
        return shortcut + noise(rng, world.n, 0.02)

    return model, hla, pep, base


def _signal_fns(world: World, model_shortcut: float, base_shortcut: float,
                model_noise: float = 0.3, pep_noise: float = 0.3):
    x_test = world.x[world.test_idx]
    rho = rho_of(world)

    def model(s, rng):
        return x_test + model_shortcut * rho + noise(rng, world.n, model_noise)

    def hla(s, rng):
        return rate_of(world)

    def pep(s, rng):                         # peptide-only sees x, as the model does
        return x_test + noise(rng, world.n, pep_noise)

    def base(s, rng):
        return x_test + base_shortcut * rho + noise(rng, world.n, 0.3)

    return model, hla, pep, base


# (i) Planted shortcut: label depends on the allele only -------------------
def test_planted_shortcut_with_scores_identical_to_the_control_is_never_passed(
    shortcut_world,
):
    w = shortcut_world
    model, hla, pep, base = _shortcut_fns(w)
    report = evaluate(w, assemble(w, model, hla, pep, base))
    v = report["variants"]["model"]
    assert v["gate"]["status"] == GateStatus.NOT_SIGNIFICANT.value, v["gate"]
    assert v["gate"]["rankable"] is False
    assert v["rank"] == UNMEASURED and report["ranking"] == []
    for name in ("model__hla", "model__pep"):
        assert report["variants"][name]["gate"]["status"] == GateStatus.CONTROL.value
        assert report["variants"][name]["rank"] == UNMEASURED
        assert report["variants"][name]["gate"]["rankable"] is False


def test_zero_skill_model_vs_lookup_control_false_pass_rate_is_bounded(shortcut_world):
    """The shortcut scenario with an HLA-only control that is a pure per-allele
    lookup (a block of ties) and a model that is the same lookup plus noise, i.e.
    it has no peptide skill. PASSED needs a significant margin AND the kappa and
    peptide-only conditions, so the rate of significant margins is an upper bound on
    the false-pass rate. Bound 0.10, as in the null-calibration test; observed rate
    in the message. A single dataset is not evidence either way, hence replicates."""
    w = shortcut_world
    allele_idx = np.unique(w.test_alleles, return_inverse=True)[1]
    counts = stats.cluster_counts(int(allele_idx.max()) + 1, 100, seed=17)
    lookup = rate_of(w)
    n_rep, hits = 200, 0
    for rep in range(n_rep):
        rng = np.random.default_rng(5000 + rep)
        model = [lookup + noise(rng, w.n, 0.02) for _ in SEEDS]
        ctrl = [lookup for _ in SEEDS]
        est = stats.margin_over_control(model, ctrl, w.test_labels, allele_idx, counts)
        hits += int(est.p_value is not None and est.p_value < 0.05)
    rate = hits / n_rep
    print(f"\nzero-skill vs lookup control: P(p<0.05) = {rate:.4f} ({hits}/{n_rep})")
    assert rate <= 0.10, f"observed false-significance rate {rate:.4f} ({hits}/{n_rep})"


# (ii) Planted signal ---------------------------------------------------------
def test_planted_signal_without_shortcut_is_detected_and_passed(signal_world):
    w = signal_world
    model, hla, pep, base = _signal_fns(w, model_shortcut=0.0, base_shortcut=0.3)
    report = evaluate(w, assemble(w, model, hla, pep, base))
    v = report["variants"]["model"]
    assert v["margin_vs_hla_only"]["significant"] is True
    assert v["gate"]["status"] == GateStatus.PASSED.value, v["gate"]
    assert v["gate"]["rankable"] is True
    assert v["rank"] == 1 and report["ranking"] == ["model"]
    change = v["shortcut_dependence"]["abs_kappa_change"]
    assert float(change["upper95_one_sided"]) <= 0.0
    # the peptide-only comparison is reported and labelled; it gates only when the
    # model is significantly WORSE
    assert v["margin_vs_peptide_only"]["significantly_worse"] is False


def test_planted_signal_with_higher_kappa_than_baseline_fails_the_shortcut_rule(
    signal_world,
):
    w = signal_world
    model, hla, pep, base = _signal_fns(w, model_shortcut=0.3, base_shortcut=0.0)
    report = evaluate(w, assemble(w, model, hla, pep, base))
    v = report["variants"]["model"]
    k_model = abs(float(v["shortcut_dependence"]["kappa_mean"]))
    k_base = abs(float(v["shortcut_dependence"]["baseline_kappa"]))
    assert k_model > k_base, "scenario precondition: the model must be the more biased"
    assert v["margin_vs_hla_only"]["significant"] is True
    assert v["gate"]["status"] == GateStatus.FAILED_SHORTCUT.value, v["gate"]
    assert "(ii)" in v["gate"]["reason"]
    assert v["rank"] == UNMEASURED and report["ranking"] == []


def test_a_model_as_shortcut_dependent_as_its_baseline_does_not_pass(signal_world):
    """The audit's measured failure, at harness level: both score alleles by
    their training rate on top of a real signal, so kappa is high for both and the
    interval of the change straddles zero. The burden of proof is on the model."""
    w = signal_world
    model, hla, pep, base = _signal_fns(w, model_shortcut=0.3, base_shortcut=0.3)
    report = evaluate(w, assemble(w, model, hla, pep, base))
    v = report["variants"]["model"]
    assert v["margin_vs_hla_only"]["significant"] is True
    assert abs(float(v["shortcut_dependence"]["kappa_mean"])) > 0.8
    assert v["gate"]["status"] != GateStatus.PASSED.value, v["gate"]
    assert float(v["shortcut_dependence"]["abs_kappa_change"]["upper95_one_sided"]) > 0
    assert report["ranking"] == []


def test_a_model_significantly_worse_than_its_peptide_only_control_is_not_passed(
    signal_world,
):
    w = signal_world
    model, hla, pep, base = _signal_fns(
        w, model_shortcut=0.0, base_shortcut=0.3, model_noise=3.0, pep_noise=0.3
    )
    report = evaluate(w, assemble(w, model, hla, pep, base))
    v = report["variants"]["model"]
    assert v["margin_vs_hla_only"]["significant"] is True
    assert v["margin_vs_peptide_only"]["significantly_worse"] is True
    assert float(v["shortcut_dependence"]["abs_kappa_change"]["upper95_one_sided"]) <= 0
    assert v["gate"]["status"] == GateStatus.FAILED_SHORTCUT.value
    assert "(iv)" in v["gate"]["reason"] and "(ii)" not in v["gate"]["reason"]
    assert v["rank"] == UNMEASURED


def test_the_kappa_change_is_a_gating_quantity_and_is_no_longer_labelled_informational(
    signal_world,
):
    w = signal_world
    report = evaluate(w, assemble(w, *_signal_fns(w, 0.0, 0.3)))
    v = report["variants"]["model"]
    assert "informational_not_used_by_gate" not in v
    change = v["shortcut_dependence"]["abs_kappa_change"]
    assert set(change) >= {"upper95_one_sided", "ci90_low", "ci90_high",
                           "noninferiority_margin", "n_valid_resamples"}
    assert "non-inferiority" in report["kappa_gate_rule"]
    rule = report["kappa_gate_rule"]
    assert "(iii)" in rule and "HLA-only control" in rule    # decision R1
    assert "ceiling" in rule and "not adopted" in rule        # the absolute ceiling


def test_the_hla_only_kappa_interval_is_reported_and_the_model_that_beats_it_passes(
    signal_world,
):
    w = signal_world
    report = evaluate(w, assemble(w, *_signal_fns(w, 0.0, 0.3)))
    v = report["variants"]["model"]
    block = v["shortcut_dependence"]["abs_kappa_change_vs_hla_only"]
    assert set(block) >= {"upper95_one_sided", "ci90_low", "ci90_high",
                          "n_valid_resamples", "required", "definition"}
    assert float(block["upper95_one_sided"]) < 0
    assert v["gate"]["status"] == GateStatus.PASSED.value, v["gate"]
    assert "HLA-only control" in v["gate"]["reason"]


def test_a_model_as_biased_as_a_saturated_baseline_is_stopped_by_condition_iii(
    signal_world,
):
    """The weakness R1 closes. Model and baseline both score alleles almost entirely by
    their training rate, so |kappa| is at its maximum for both and the change against
    the baseline has a degenerate interval at zero, which satisfies (ii). The HLA-only
    control is shortcut-dependent by construction (kappa at its maximum), so the change
    against it is also zero and (iii) fails: the model is not shown to be LESS
    shortcut-dependent than the control."""
    w = signal_world
    model, hla, pep, base = _signal_fns(w, model_shortcut=50.0, base_shortcut=50.0)
    report = evaluate(w, assemble(w, model, hla, pep, base))
    v = report["variants"]["model"]
    assert v["margin_vs_hla_only"]["significant"] is True
    change = v["shortcut_dependence"]["abs_kappa_change"]
    assert float(change["upper95_one_sided"]) <= 0.0, "scenario precondition: (ii) is met"
    assert v["gate"]["status"] == GateStatus.FAILED_SHORTCUT.value, v["gate"]
    assert "(iii)" in v["gate"]["reason"] and "(ii)" not in v["gate"]["reason"]
    assert v["rank"] == UNMEASURED and report["ranking"] == []


def test_an_undefined_hla_only_kappa_withholds_the_pass(signal_world):
    """An HLA-only control that scores every record the same has no defined kappa; the
    interval for (iii) is then unmeasured and the model is PENDING, not passed."""
    w = signal_world
    model, _, pep, base = _signal_fns(w, 0.0, 0.3)

    def flat(s, rng):
        return np.zeros(w.n)

    report = evaluate(w, assemble(w, model, flat, pep, base))
    v = report["variants"]["model"]
    assert v["gate"]["status"] == GateStatus.PENDING.value, v["gate"]
    assert v["gate"]["rankable"] is False and report["ranking"] == []
    assert "HLA-only" in v["harness_note"]


def test_missing_baseline_withholds_the_pass_instead_of_comparing_to_nothing(
    signal_world,
):
    w = signal_world
    model, hla, pep, _ = _signal_fns(w, 0.0, 0.0)
    report = evaluate(w, assemble(w, model, hla, pep))          # no baseline run
    v = report["variants"]["model"]
    assert v["gate"]["status"] == GateStatus.PENDING.value
    assert v["gate"]["rankable"] is False and report["ranking"] == []
    assert v["metrics"]["displayed"] is False                    # FR-63
    assert "pooled" not in v["metrics"]
    assert "baseline" in v["harness_note"].lower()


# (iii) Null calibration ---------------------------------------------------------
def test_null_calibration_false_significance_rate_is_near_alpha():
    """Independent noise for model and control, so there is no true margin.

    Bound: the fraction of replicates with p < 0.05 must be <= 0.10, not 0.05.
    Reasons: with R replicates the observed rate has standard error
    sqrt(.05*.95/R) (about 0.015 at R=200), and a percentile-bootstrap p-value over a
    modest number of clusters is only approximately uniform under the null (it is
    known to be somewhat liberal for few clusters). 0.10 is the bound the project
    set in advance; it is not tuned to this result. If it fails, the observed rate
    is the finding. The observed rate is printed (pytest -s) and in the message."""
    n_alleles, per, n_rep, n_boot = 20, 20, 240, 200
    labels = np.tile([1, 0], n_alleles * per // 2)
    allele_idx = np.repeat(np.arange(n_alleles), per)
    counts = stats.cluster_counts(n_alleles, n_boot, seed=2026)
    hits = 0
    for rep in range(n_rep):
        rng = np.random.default_rng(10_000 + rep)
        model = [rng.normal(size=labels.size) for _ in SEEDS]
        ctrl = [rng.normal(size=labels.size) for _ in SEEDS]
        est = stats.margin_over_control(model, ctrl, labels, allele_idx, counts)
        hits += int(est.p_value is not None and est.p_value < 0.05)
    rate = hits / n_rep
    print(f"\nnull calibration: P(p<0.05) = {rate:.4f} ({hits}/{n_rep} replicates)")
    assert rate <= 0.10, f"observed false-pass rate {rate:.4f} ({hits}/{n_rep})"


def test_margin_power_is_not_trivially_zero():
    """The converse of calibration: a real margin is found."""
    n_alleles, per = 20, 20
    labels = np.tile([1, 0], n_alleles * per // 2)
    allele_idx = np.repeat(np.arange(n_alleles), per)
    counts = stats.cluster_counts(n_alleles, 200, seed=1)
    rng = np.random.default_rng(9)
    model = [labels * 1.0 + rng.normal(size=labels.size) for _ in SEEDS]
    ctrl = [rng.normal(size=labels.size) for _ in SEEDS]
    est = stats.margin_over_control(model, ctrl, labels, allele_idx, counts)
    assert est.delta > 0 and est.p_value is not None and est.p_value < 0.05
    assert est.ci90_low > 0
    assert est.one_sided_lower95 == pytest.approx(est.ci90_low)
    assert est.significantly_worse is False
    worse = stats.margin_over_control(ctrl, model, labels, allele_idx, counts)
    assert worse.significantly_worse is True and worse.p_value > 0.5


# ==========================================================================
# Reporting rules
# ==========================================================================
@pytest.fixture(scope="module")
def passed_report(signal_world):
    w = signal_world
    fns = _signal_fns(w, 0.0, 0.3)
    return evaluate(w, assemble(w, *fns))


def test_every_model_metric_carries_its_margin_and_controls_are_labelled(
    passed_report,
):
    m = passed_report["variants"]["model"]["metrics"]
    assert m["displayed"] is True
    for group in ("pooled", "allele_stratified_support_weighted",
                  "allele_stratified_macro"):
        for name, entry in m[group].items():
            assert re.fullmatch(r"[+-]\d\.\d{4} \((significant|n\.s\.)\)",
                                entry["margin_over_hla_only"]), (group, name)
    ctrl = passed_report["variants"]["model__hla"]
    assert ctrl["gate"]["status"] == "control" and ctrl["rank"] == UNMEASURED
    assert ctrl["metrics"]["pooled"]["auprc"]["margin_over_hla_only"] == "control row"


def test_the_test_and_the_interval_are_a_matching_one_sided_pair(passed_report):
    block = passed_report["variants"]["model"]["margin_vs_hla_only"]
    assert "ci95_low" not in block and "ci95_high" not in block
    assert {"p_value", "ci90_low", "ci90_high", "one_sided_lower95"} <= set(block)
    assert block["one_sided_lower95"] == block["ci90_low"]
    text = passed_report["bootstrap"]["interval"]
    assert "90%" in text and "5th" in text and "95th" in text
    assert "mismatch" not in text and "can co-occur" not in text
    assert "0.950000" not in json.dumps(passed_report["variants"]["model"]
                                        ["margin_vs_hla_only"])


def test_headline_is_mean_and_sd_across_seeds_never_the_best_seed(signal_world):
    w = signal_world
    fns = _signal_fns(w, 0.0, 0.3)
    runs = assemble(w, *fns)
    report = evaluate(w, runs)
    per_seed = []
    for r in runs:
        if r.variant == "model":
            sc = np.array([r.scores[int(i)] for i in w.test_idx])
            per_seed.append(stats.expected_average_precision(sc, w.test_labels))
    entry = report["variants"]["model"]["metrics"]["pooled"]["auprc"]
    assert float(entry["mean"]) == pytest.approx(np.mean(per_seed), abs=5e-6)
    assert float(entry["sd"]) == pytest.approx(np.std(per_seed, ddof=1), abs=5e-6)
    text = report_to_json(report)
    assert not re.search(r'"(best|max|min)[^"]*"', text)
    assert '"per_seed"' not in text          # only the aggregate, never a seed list


def test_report_has_no_null_and_unmeasured_is_an_em_dash(signal_world):
    w = signal_world
    fns = _signal_fns(w, 0.0, 0.3)
    report = evaluate(w, assemble(w, *fns), ranking_k=w.n + 1)   # PPV@k undefined
    text = report_to_json(report)

    def no_none(node):
        if isinstance(node, dict):
            return all(no_none(v) for v in node.values())
        if isinstance(node, list):
            return all(no_none(v) for v in node)
        return node is not None

    assert no_none(json.loads(text)) and ": null" not in text
    assert report["variants"]["model"]["metrics"]["pooled"]["ppv_at_k"]["mean"] == (
        UNMEASURED
    )
    assert UNMEASURED in text


def test_report_is_byte_reproducible_and_content_hashed(signal_world):
    w = signal_world
    fns = _signal_fns(w, 0.0, 0.3)
    runs = assemble(w, *fns)
    a, b = evaluate(w, runs), evaluate(w, runs)
    assert report_to_json(a) == report_to_json(b)
    assert report_to_json(a).encode("utf-8") == report_to_json(b).encode("utf-8")
    assert a["report_hash"] == b["report_hash"] and len(a["report_hash"]) == 64
    c = evaluate(w, runs, bootstrap_seed=8)
    assert c["report_hash"] != a["report_hash"]            # the seed is really used
    assert a["bootstrap"]["seed"] == 7 and a["bootstrap"]["B"] == N_BOOT
    assert "tie" in a["bootstrap"]["tie_handling"] and "seed" in a["bootstrap"][
        "tie_handling"
    ] and "no seed" in a["bootstrap"]["tie_handling"]


def test_report_states_weighting_split_hash_seeds_and_config_hashes(passed_report):
    r = passed_report
    assert r["seeds"] == list(SEEDS) and r["split_hash"]
    assert "support" in r["weighting"]["allele_stratified_support_weighted"].lower()
    assert "equal" in r["weighting"]["allele_stratified_macro"].lower()
    h = {v["config_hash"] for k, v in r["variants"].items() if k.startswith("model")}
    assert len(h) == 1                                     # matched budget, visibly
    assert r["bootstrap"]["unit"].startswith("test alleles")
    assert any("FR-65" in s for s in r["not_done"])
    assert r["test_set"]["n_records"] > 0 and r["test_set"]["n_alleles"] > 0


def test_alleles_outside_the_admitted_test_set_are_reported_separately_never_pooled():
    w = make_world("signal", extra_small_allele=True)
    fns = _signal_fns(w, 0.0, 0.3)
    report = evaluate(w, assemble(w, *fns))
    sep = report["not_admitted"]
    assert any(e["allele"] == "SMALL" for e in sep["alleles"])
    assert sep["metrics"] == UNMEASURED
    assert report["test_set"]["n_records"] == w.n
    assert "SMALL" not in report["test_set"]["alleles"]


def test_matching_block_separates_declared_metadata_from_structurally_verified(
    passed_report,
):
    m = passed_report["matching"]
    assert set(m["declared_metadata"]) >= {"epochs", "n_parameters",
                                           "early_stopping_metric", "tuning_trials"}
    verified = m["structurally_verified"]
    assert "constant within every test allele" in verified["hla_only"]["rule"]
    assert verified["hla_only"]["n_runs_checked"] > 0
    assert "peptide" in verified["peptide_only"]["rule"]
    assert "not verified" in m["not_verified"].lower()
    assert m["n_parameters_rel_tol"] == "0.000000"


def test_report_states_the_number_of_model_variants_and_holm_adjusted_p(signal_world):
    w = signal_world
    runs = assemble(w, *_signal_fns(w, 0.0, 0.3))
    # a second model variant with its own matched controls
    more = assemble(w, *_signal_fns(w, 0.0, 0.3, model_noise=0.8)[:3], model="model2")
    report = evaluate(w, runs + more)
    assert report["n_model_variants_compared"] == 2
    p1 = float(report["variants"]["model"]["margin_vs_hla_only"]["p_value"])
    p2 = float(report["variants"]["model2"]["margin_vs_hla_only"]["p_value"])
    adj = [float(report["variants"][n]["margin_vs_hla_only"]["p_value_holm"])
           for n in ("model", "model2")]
    assert adj == pytest.approx(stats.holm_adjust([p1, p2]), abs=2e-6)
    assert all(a >= p - 1e-9 for a, p in zip(adj, (p1, p2)))
    assert "Holm" in report["multiplicity"]
    # the baseline is not a compared model variant
    assert "p_value_holm" not in report["variants"]["base"]["margin_vs_hla_only"]


def test_preregistered_settings_flag_is_false_whenever_n_boot_is_lowered(
    signal_world, monkeypatch
):
    w = signal_world
    runs = assemble(w, *_signal_fns(w, 0.0, 0.3))
    assert evaluate(w, runs)["preregistered_settings_used"] is False
    monkeypatch.setattr(protocol, "PRODUCTION_N_BOOT", N_BOOT)
    ok = evaluate(w, runs, bootstrap_seed=DEFAULT_BOOTSTRAP_SEED,
                  ranking_k=DEFAULT_RANKING_K)
    assert ok["preregistered_settings_used"] is True
    assert evaluate(w, runs, bootstrap_seed=1,
                    ranking_k=DEFAULT_RANKING_K)["preregistered_settings_used"] is False


def test_the_production_constants_are_the_analysis_plan_values():
    assert PRODUCTION_N_BOOT == 10_000
    assert set(EXPECTED_SEEDS) == {1, 2, 3, 4, 5} and MIN_SEEDS == 5


# ==========================================================================
# (iv) Refusals
# ==========================================================================
def _good_runs(world: World) -> list[VariantRun]:
    return assemble(world, *_signal_fns(world, 0.0, 0.3))


def _write(tmp_path: Path, world: World, runs) -> Path:
    write_run_dir(tmp_path / "run", world.record, runs)
    return tmp_path / "run"


def _refuse(tmp_path, world, runs, match, hash_override=None, capsys=None):
    run_dir = _write(tmp_path, world, runs)
    h = hash_override or world.record.content_hash
    with pytest.raises(ProtocolRefusal, match=match) as exc:
        run_protocol(run_dir, h, n_boot=N_BOOT, bootstrap_seed=7, ranking_k=K)
    assert not re.search(NO_METRIC_TEXT, str(exc.value), re.I)
    out = capsys.readouterr() if capsys else None
    if out is not None:
        assert out.out == "" and out.err == ""
    assert not list(run_dir.glob("report*"))             # nothing was emitted


def test_valid_layout_round_trips_through_disk(tmp_path, signal_world):
    w = signal_world
    runs = _good_runs(w)
    run_dir = _write(tmp_path, w, runs)
    from_disk = run_protocol(
        run_dir, w.record.content_hash, n_boot=N_BOOT, bootstrap_seed=7, ranking_k=K
    )
    in_memory = evaluate(w, runs)
    assert report_to_json(from_disk) == report_to_json(in_memory)


def test_refuses_when_split_hash_differs(tmp_path, signal_world, capsys):
    _refuse(tmp_path, signal_world, _good_runs(signal_world), "split hash",
            hash_override="0" * 64, capsys=capsys)


def test_refuses_when_split_file_does_not_match_its_declared_hash(
    tmp_path, signal_world, capsys
):
    w = signal_world
    run_dir = _write(tmp_path, w, _good_runs(w))
    path = run_dir / "split.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc["test_idx"] = doc["test_idx"][:-1]               # tampered, hash unchanged
    path.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(ProtocolRefusal, match="does not match"):
        run_protocol(run_dir, w.record.content_hash, n_boot=N_BOOT, ranking_k=K)
    assert capsys.readouterr().out == ""


# ---- D16(e): records digest, tamper, peptide overlap ----------------------------
@pytest.mark.parametrize("column", ["label", "allele", "peptide"])
def test_refuses_tampered_records_under_an_unchanged_split_hash(
    tmp_path, signal_world, column, capsys
):
    """The tamper case the verifier found: flipped test labels with an unchanged
    split hash. The FrozenSplit hash does not cover labels; the records digest does."""
    w = signal_world
    run_dir = _write(tmp_path, w, _good_runs(w))
    path = run_dir / "split.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    i = int(doc["test_idx"][0])
    col = doc["records"][column]
    col[i] = (1 - col[i]) if column == "label" else col[i] + "X"
    path.write_text(json.dumps(doc), encoding="utf-8")
    assert doc["content_hash"] == w.record.content_hash        # hash really unchanged
    with pytest.raises(ProtocolRefusal, match="records digest"):
        run_protocol(run_dir, w.record.content_hash, n_boot=N_BOOT, ranking_k=K)
    out = capsys.readouterr()
    assert out.out == "" and out.err == ""


def test_a_coordinated_edit_is_caught_only_when_the_digest_is_pinned(
    tmp_path, signal_world
):
    w = signal_world
    labels = list(w.record.labels)
    labels[w.record.test_idx[0]] = 1 - labels[w.record.test_idx[0]]
    forged = w.record.with_records(labels=labels)      # digest recomputed: consistent
    assert forged.content_hash == w.record.content_hash
    assert forged.records_sha256 != w.record.records_sha256
    run_dir = tmp_path / "forged"
    write_run_dir(run_dir, forged, _good_runs(w))
    # nothing inside the directory is inconsistent, so without a pinned digest it loads
    run_protocol(run_dir, w.record.content_hash, n_boot=N_BOOT, bootstrap_seed=7,
                 ranking_k=K)
    with pytest.raises(ProtocolRefusal, match="pinned"):
        run_protocol(run_dir, w.record.content_hash, n_boot=N_BOOT, bootstrap_seed=7,
                     ranking_k=K, expected_records_digest=w.record.records_sha256)


def test_report_records_the_digest_and_whether_it_was_pinned(signal_world):
    w = signal_world
    runs = _good_runs(w)
    r = evaluate(w, runs)
    assert r["records"]["sha256"] == w.record.records_sha256
    assert r["records"]["verified_against_split_json"] is True
    assert r["records"]["pinned_by_caller"] is False
    pinned = evaluate(w, runs, expected_records_digest=w.record.records_sha256)
    assert pinned["records"]["pinned_by_caller"] is True


def test_refuses_when_a_peptide_occurs_in_more_than_one_partition(
    tmp_path, signal_world, capsys
):
    w = signal_world
    peptides = list(w.record.peptides)
    peptides[w.record.test_idx[0]] = peptides[w.record.train_idx[0]]   # leak
    leaky = w.record.with_records(peptides=peptides)
    assert leaky.content_hash == w.record.content_hash
    _refuse(tmp_path, dataclasses.replace(w, record=leaky), _good_runs(w),
            "peptide", capsys=capsys)


def test_split_json_carries_per_record_peptide_allele_and_label(tmp_path, signal_world):
    w = signal_world
    run_dir = _write(tmp_path, w, _good_runs(w))
    doc = json.loads((run_dir / "split.json").read_text(encoding="utf-8"))
    n = doc["n_records"]
    assert {len(doc["records"][k]) for k in ("allele", "label", "peptide")} == {n}
    assert len(doc["records_sha256"]) == 64


# ---- D16(a): seeds pinned -----------------------------------------------------
def test_refuses_with_fewer_than_the_expected_seeds(tmp_path, signal_world, capsys):
    runs = [r for r in _good_runs(signal_world) if r.seed != 5]
    _refuse(tmp_path, signal_world, runs, "seeds", capsys=capsys)


def test_refuses_an_extra_seed_on_every_variant(tmp_path, signal_world, capsys):
    """Training six seeds and submitting five is a forking path; six is refused."""
    w = signal_world
    runs = list(_good_runs(w))
    runs += [dataclasses.replace(r, seed=6) for r in runs if r.seed == 1]
    _refuse(tmp_path, w, runs, "seeds", capsys=capsys)


def test_refuses_five_seeds_that_are_not_the_preregistered_five(
    tmp_path, signal_world, capsys
):
    w = signal_world
    runs = [dataclasses.replace(r, seed=r.seed + 1) for r in _good_runs(w)]
    _refuse(tmp_path, w, runs, "seeds", capsys=capsys)


# ---- D16(d): structural control validation ----------------------------------------
def test_refuses_an_hla_only_control_that_varies_within_an_allele(
    tmp_path, signal_world, capsys
):
    """A noise model labelled hla_only would let a shortcut model pass."""
    w = signal_world
    model, _, pep, base = _signal_fns(w, 0.0, 0.3)

    def noisy_hla(s, rng):
        return rate_of(w) + noise(rng, w.n, 0.02)

    runs = assemble(w, model, noisy_hla, pep, base)
    _refuse(tmp_path, w, runs, "structural", capsys=capsys)


def _record_with_shared_test_peptide(w: World) -> tuple[SplitRecord, int, int]:
    first = int(w.test_idx[0])
    other = next(int(i) for i in w.test_idx
                 if w.record.alleles[int(i)] != w.record.alleles[first])
    peptides = list(w.record.peptides)
    peptides[first] = peptides[other] = "SHAREDPEPTIDE"
    return w.record.with_records(peptides=peptides), first, other


def test_refuses_a_peptide_only_control_that_scores_a_shared_peptide_differently(
    tmp_path, signal_world, capsys
):
    w = signal_world
    record, _, _ = _record_with_shared_test_peptide(w)
    world = dataclasses.replace(w, record=record)
    _refuse(tmp_path, world, _good_runs(w), "structural", capsys=capsys)


def test_accepts_a_peptide_only_control_that_is_invariant_to_the_allele(signal_world):
    w = signal_world
    record, first, other = _record_with_shared_test_peptide(w)
    pos = {int(i): j for j, i in enumerate(w.test_idx)}
    model, hla, pep, base = _signal_fns(w, 0.0, 0.3)

    def invariant_pep(s, rng):
        scores = pep(s, rng)
        scores[pos[other]] = scores[pos[first]]
        return scores

    report = evaluate(w, assemble(w, model, hla, invariant_pep, base,
                                  base_pep_fn=invariant_pep), record=record)
    verified = report["matching"]["structurally_verified"]["peptide_only"]
    assert verified["n_shared_peptide_groups"] == 1
    assert verified["n_runs_checked"] == 2 * len(SEEDS)


def test_the_structural_peptide_check_states_when_it_is_vacuous(passed_report):
    verified = passed_report["matching"]["structurally_verified"]["peptide_only"]
    assert verified["n_shared_peptide_groups"] == 0
    assert "vacuous" in verified["note"].lower()


# ---- the other FR-62 / FR-63 refusals --------------------------------------------
@pytest.mark.parametrize("missing_role", ["hla_only", "peptide_only"])
def test_refuses_when_a_model_has_no_control_of_a_mandated_kind(
    tmp_path, signal_world, missing_role, capsys
):
    runs = [r for r in _good_runs(signal_world)
            if not (r.matched_to == "model" and r.role == missing_role)]
    _refuse(tmp_path, signal_world, runs, "control", capsys=capsys)


@pytest.mark.parametrize("missing_role", ["hla_only", "peptide_only"])
def test_refuses_when_a_control_lacks_one_seed_the_model_reports(
    tmp_path, signal_world, missing_role, capsys
):
    runs = [r for r in _good_runs(signal_world)
            if not (r.matched_to == "model" and r.role == missing_role and r.seed == 3)]
    _refuse(tmp_path, signal_world, runs, "seed", capsys=capsys)


def test_refuses_when_a_control_reports_a_seed_the_model_does_not(
    tmp_path, signal_world, capsys
):
    w = signal_world
    runs = _good_runs(w)
    extra = next(r for r in runs if r.variant == "model__hla" and r.seed == 1)
    runs.append(dataclasses.replace(extra, seed=6))
    _refuse(tmp_path, w, runs, "seed", capsys=capsys)


@pytest.mark.parametrize("field,value", [
    ("epochs", 11),
    ("n_parameters", 1001),
    ("early_stopping_metric", "val_auroc"),
    ("tuning_trials", 6),
])
def test_refuses_when_control_config_differs_from_the_model(
    tmp_path, signal_world, field, value, capsys
):
    runs = [
        dataclasses.replace(r, config=dataclasses.replace(r.config, **{field: value}))
        if r.variant == "model__hla" else r
        for r in _good_runs(signal_world)
    ]
    _refuse(tmp_path, signal_world, runs, field, capsys=capsys)


def test_parameter_count_tolerance_is_off_by_default_and_recorded_when_used(
    signal_world,
):
    w = signal_world
    runs = [
        dataclasses.replace(r, config=dataclasses.replace(r.config, n_parameters=1001))
        if r.variant == "model__hla" else r
        for r in _good_runs(w)
    ]
    with pytest.raises(ProtocolRefusal, match="n_parameters"):
        evaluate(w, runs)
    report = evaluate(w, runs, n_parameters_rel_tol=0.01)
    assert report["matching"]["n_parameters_rel_tol"] == "0.010000"
    assert evaluate(w, _good_runs(w))["matching"]["n_parameters_rel_tol"] == "0.000000"


def test_refuses_when_test_indices_differ_between_runs(
    tmp_path, signal_world, capsys
):
    w = signal_world
    train_i = int(w.record.train_idx[0])
    runs = []
    for r in _good_runs(w):
        if r.variant == "model__pep" and r.seed == 2:
            scores = dict(r.scores)
            dropped = next(iter(scores))
            scores.pop(dropped)
            scores[train_i] = 0.5                     # a record outside the test set
            r = dataclasses.replace(r, scores=scores)
        runs.append(r)
    _refuse(tmp_path, w, runs, "test indices differ", capsys=capsys)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_refuses_non_finite_predictions(tmp_path, signal_world, bad, capsys):
    runs = []
    for r in _good_runs(signal_world):
        if r.variant == "model" and r.seed == 4:
            scores = dict(r.scores)
            scores[next(iter(scores))] = bad
            r = dataclasses.replace(r, scores=scores)
        runs.append(r)
    _refuse(tmp_path, signal_world, runs, "finite|NaN|inf", capsys=capsys)


def test_refuses_when_a_score_is_missing_for_a_test_record(
    tmp_path, signal_world, capsys
):
    runs = []
    for r in _good_runs(signal_world):
        if r.variant == "model" and r.seed == 1:
            scores = dict(r.scores)
            scores.pop(next(iter(scores)))
            r = dataclasses.replace(r, scores=scores)
        runs.append(r)
    _refuse(tmp_path, signal_world, runs, "missing a score", capsys=capsys)


def test_refuses_a_model_whose_controls_are_attached_to_nothing(
    tmp_path, signal_world, capsys
):
    runs = [
        dataclasses.replace(r, matched_to="ghost") if r.variant == "model__hla" else r
        for r in _good_runs(signal_world)
    ]
    _refuse(tmp_path, signal_world, runs, "control", capsys=capsys)


def test_baseline_must_be_matched_too(tmp_path, signal_world, capsys):
    runs = [r for r in _good_runs(signal_world) if r.variant != "base__hla"]
    _refuse(tmp_path, signal_world, runs, "control", capsys=capsys)


# ---- the CLI --------------------------------------------------------------------
def test_cli_exits_nonzero_and_prints_no_metric_on_refusal(tmp_path, signal_world):
    w = signal_world
    runs = [r for r in _good_runs(w) if r.variant != "model__hla"]
    run_dir = _write(tmp_path, w, runs)
    result = CliRunner().invoke(
        app, ["eval", "--protocol", "bias-controlled", "--run", str(run_dir),
              "--split", w.record.content_hash],
    )
    assert result.exit_code == 1
    assert "refused" in result.output.lower()
    assert not re.search(NO_METRIC_TEXT, result.output, re.I)


def test_cli_has_no_way_to_lower_the_bootstrap_or_change_its_seed():
    params = set(inspect.signature(__import__("nbx.cli", fromlist=["evaluate"])
                                   .evaluate).parameters)
    assert params == {"run", "split", "protocol", "out", "records_digest"}


def _spy_protocol(monkeypatch):
    captured: dict = {}
    real = protocol.run_protocol

    def spy(run_dir, split_hash, **kw):
        captured.update(kw)
        # the only place the production bootstrap is lowered: this test harness
        return real(run_dir, split_hash, **{**kw, "n_boot": N_BOOT, "ranking_k": K})

    monkeypatch.setattr(protocol, "run_protocol", spy)
    return captured


def test_cli_always_requests_the_production_bootstrap_and_fixed_seed(
    tmp_path, signal_world, monkeypatch
):
    captured = _spy_protocol(monkeypatch)
    w = signal_world
    run_dir = _write(tmp_path, w, _good_runs(w))
    out = tmp_path / "report.json"
    result = CliRunner().invoke(
        app, ["eval", "--protocol", "bias-controlled", "--run", str(run_dir),
              "--split", w.record.content_hash, "--out", str(out)],
    )
    assert result.exit_code == 0, result.output
    assert captured["n_boot"] == PRODUCTION_N_BOOT == 10_000
    assert captured["bootstrap_seed"] == DEFAULT_BOOTSTRAP_SEED
    written = json.loads(out.read_text(encoding="utf-8"))
    assert written["split_hash"] == w.record.content_hash
    assert written["preregistered_settings_used"] is False      # the spy lowered B


def test_cli_prints_the_presentation_the_researcher_needs(
    tmp_path, signal_world, monkeypatch
):
    _spy_protocol(monkeypatch)
    w = signal_world
    run_dir = _write(tmp_path, w, _good_runs(w))
    result = CliRunner().invoke(
        app, ["eval", "--protocol", "bias-controlled", "--run", str(run_dir),
              "--split", w.record.content_hash,
              "--records-digest", w.record.records_sha256],
    )
    text = " ".join(result.output.split())
    assert result.exit_code == 0, result.output
    assert "passed" in text and "control" in text
    assert "reference_only" in text                            # the baseline row
    assert "baseline kappa" in text
    assert "kappa change upper95" in text
    assert "peptide-only margin" in text and "90% interval" in text
    assert "test alleles" in text and "test records" in text
    assert "declared" in text and "structurally verified" in text
    assert "records digest" in text and "verified" in text and "pinned" in text
    assert "Holm" in text
    assert "settings NOT pre-registered" in text               # B was lowered here


def test_cli_rejects_an_unknown_protocol(tmp_path):
    result = CliRunner().invoke(
        app, ["eval", "--protocol", "other", "--run", str(tmp_path), "--split", "x"]
    )
    assert result.exit_code != 0


# ==========================================================================
# (vi) There is no way to switch a control off
# ==========================================================================
FORBIDDEN = {"force", "override", "skip_controls", "skip_control", "ignore_kappa",
             "bypass_gate", "disable_controls", "allow_unranked"}


def test_evaluate_package_has_no_bypass_parameters():
    root = Path(__file__).resolve().parent.parent / "src" / "nbx" / "evaluate"
    files = sorted(root.glob("*.py"))
    assert {f.name for f in files} >= {"gate.py", "metrics.py", "protocol.py", "stats.py"}
    offenders = []
    for path in files:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda)):
                a = node.args
                names = [x.arg for x in (*a.posonlyargs, *a.args, *a.kwonlyargs)]
                names += [x.arg for x in (a.vararg, a.kwarg) if x is not None]
                offenders += [(path.name, n) for n in names if n in FORBIDDEN]
            if isinstance(node, ast.Call):
                offenders += [(path.name, k.arg) for k in node.keywords
                              if k.arg in FORBIDDEN]
    assert offenders == []


def test_protocol_signature_has_no_way_to_name_a_control_out():
    params = set(inspect.signature(run_protocol).parameters)
    assert not (params & FORBIDDEN)
    assert {"run_dir", "split_hash"} <= params
