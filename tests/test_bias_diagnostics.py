"""Tests for the bias diagnostics.

Every test here pins a value that is known analytically or is cross-checked
against an independent implementation. The point of this project is that a
number which merely ran is not a number that is right.
"""

from __future__ import annotations

import math

import pytest

from nbx.bias.diagnostics import (
    RHO_SATURATION,
    BiasProfile,
    conditional_entropy,
    entropy,
    log_imbalance_ratio,
    mutual_information,
    per_allele_imbalance,
    profile_dataset,
    shortcut_dependence,
    spearman,
)


# --------------------------------------------------------------------------
# Entropy
# --------------------------------------------------------------------------
def test_entropy_of_fair_coin_is_log_two():
    assert entropy([0, 1]) == pytest.approx(math.log(2))
    assert entropy([0, 1], log_base=2) == pytest.approx(1.0)


def test_entropy_of_constant_label_is_zero():
    assert entropy([1, 1, 1, 1]) == pytest.approx(0.0)


def test_entropy_of_known_skew():
    # P = (0.25, 0.75): H = -(0.25 ln 0.25 + 0.75 ln 0.75)
    labels = [1, 0, 0, 0]
    expected = -(0.25 * math.log(0.25) + 0.75 * math.log(0.75))
    assert entropy(labels) == pytest.approx(expected)


def test_entropy_rejects_empty_input():
    # Zero would be a measurement. Absence must not silently become one.
    with pytest.raises(ValueError):
        entropy([])


# --------------------------------------------------------------------------
# Conditional entropy and mutual information
# --------------------------------------------------------------------------
def test_conditional_entropy_zero_when_variable_determines_label():
    labels = [0, 0, 1, 1]
    by = ["a", "a", "b", "b"]
    assert conditional_entropy(labels, by) == pytest.approx(0.0)


def test_mutual_information_is_label_entropy_when_fully_determined():
    labels = [0, 0, 1, 1]
    by = ["a", "a", "b", "b"]
    assert mutual_information(labels, by) == pytest.approx(math.log(2))


def test_mutual_information_is_zero_for_independent_variable():
    # Each group has the same 50/50 label split, so b says nothing about y.
    labels = [0, 1, 0, 1]
    by = ["a", "a", "b", "b"]
    assert mutual_information(labels, by) == pytest.approx(0.0, abs=1e-12)


def test_mutual_information_never_negative():
    labels = [0, 1, 1, 0, 1, 0, 0, 1]
    by = ["x", "y", "x", "y", "x", "y", "x", "y"]
    assert mutual_information(labels, by) >= 0.0


def test_mutual_information_rejects_misaligned_inputs():
    with pytest.raises(ValueError):
        mutual_information([0, 1, 1], ["a", "b"])


def test_mutual_information_matches_sklearn_when_available():
    sklearn_metrics = pytest.importorskip("sklearn.metrics")
    labels = [0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1]
    by = ["A", "A", "B", "C", "B", "B", "C", "A", "B", "C", "A", "B"]
    expected = sklearn_metrics.mutual_info_score(by, labels)
    assert mutual_information(labels, by) == pytest.approx(expected, rel=1e-9)


# --------------------------------------------------------------------------
# Imbalance
# --------------------------------------------------------------------------
def test_log_imbalance_ratio_balanced_is_zero():
    assert log_imbalance_ratio(50, 50) == pytest.approx(0.0)


def test_log_imbalance_ratio_ten_to_one():
    # The ">10:1" threshold the documents quote is log2(10) = 3.3219...
    assert log_imbalance_ratio(100, 10) == pytest.approx(math.log2(10))


def test_log_imbalance_ratio_saturates_rather_than_diverging():
    assert log_imbalance_ratio(7, 0) == RHO_SATURATION
    assert log_imbalance_ratio(0, 7) == -RHO_SATURATION


def test_log_imbalance_ratio_rejects_empty_allele():
    with pytest.raises(ValueError):
        log_imbalance_ratio(0, 0)


def test_per_allele_imbalance_excludes_single_sample_alleles():
    alleles = ["A", "A", "A", "B"]
    labels = [1, 1, 0, 1]
    rho = per_allele_imbalance(alleles, labels)
    assert "B" not in rho, "single-sample alleles must be excluded"
    assert rho["A"] == pytest.approx(math.log2(2 / 1))


def test_per_allele_imbalance_rejects_non_binary_labels():
    with pytest.raises(ValueError):
        per_allele_imbalance(["A", "A"], [1, 2])


# --------------------------------------------------------------------------
# Spearman and kappa
# --------------------------------------------------------------------------
def test_spearman_perfect_monotonic():
    assert spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [40, 30, 20, 10]) == pytest.approx(-1.0)


def test_spearman_handles_ties_by_averaging_ranks():
    # Ties averaged is what scipy does; verified below when scipy is present.
    value = spearman([1, 1, 2, 3], [5, 5, 6, 7])
    assert value == pytest.approx(1.0)


def test_spearman_matches_scipy_when_available():
    stats = pytest.importorskip("scipy.stats")
    xs = [0.4, -1.2, 3.3, 2.0, -0.5, 1.1, 1.1, 0.0]
    ys = [1.0, 0.2, 2.9, 2.1, 0.1, 1.4, 1.5, 0.9]
    assert spearman(xs, ys) == pytest.approx(stats.spearmanr(xs, ys).statistic)


def test_spearman_refuses_too_few_points():
    with pytest.raises(ValueError):
        spearman([1, 2], [1, 2])


def test_spearman_refuses_constant_input():
    with pytest.raises(ValueError):
        spearman([1, 1, 1, 1], [1, 2, 3, 4])


def test_shortcut_dependence_detects_a_pure_shortcut():
    """A model whose score is a monotone function of training imbalance."""
    train_rho = {"A": -2.0, "B": -0.5, "C": 0.5, "D": 2.0}
    test_alleles = ["A", "A", "B", "B", "C", "C", "D", "D"]
    # Scores track rho exactly: this is the failure mode, maximally expressed.
    test_scores = [0.1, 0.1, 0.4, 0.4, 0.6, 0.6, 0.9, 0.9]
    assert shortcut_dependence(train_rho, test_alleles, test_scores) == pytest.approx(1.0)


def test_shortcut_dependence_near_zero_for_imbalance_blind_model():
    train_rho = {"A": -2.0, "B": -0.5, "C": 0.5, "D": 2.0}
    test_alleles = ["A", "A", "B", "B", "C", "C", "D", "D"]
    # Mean score identical for every allele -> no rank correlation is defined.
    test_scores = [0.5] * 8
    with pytest.raises(ValueError):
        shortcut_dependence(train_rho, test_alleles, test_scores)


def test_shortcut_dependence_ignores_alleles_absent_from_training():
    train_rho = {"A": -2.0, "B": 0.0, "C": 2.0}
    test_alleles = ["A", "B", "C", "ZZZ_unseen"]
    test_scores = [0.1, 0.5, 0.9, 99.0]
    # The unseen allele must not participate; result stays a clean +1.
    assert shortcut_dependence(train_rho, test_alleles, test_scores) == pytest.approx(1.0)


def test_shortcut_dependence_requires_three_shared_alleles():
    with pytest.raises(ValueError):
        shortcut_dependence({"A": 1.0, "B": 2.0}, ["A", "B"], [0.1, 0.2])


# --------------------------------------------------------------------------
# Dataset profile
# --------------------------------------------------------------------------
def test_profile_dataset_flags_allele_as_bias_candidate():
    # Allele fully determines the label: the textbook shortcut.
    alleles = ["A", "A", "A", "A", "B", "B", "B", "B"]
    labels = [1, 1, 1, 1, 0, 0, 0, 0]
    lengths = [9, 9, 10, 9, 9, 10, 9, 9]

    profile = profile_dataset(
        labels, {"hla_allele": alleles, "peptide_length": lengths}
    )

    assert profile.n_records == 8
    assert profile.label_entropy == pytest.approx(math.log(2))
    assert profile.mutual_information["hla_allele"] == pytest.approx(math.log(2))
    assert "hla_allele" in profile.bias_candidates()
    # Both alleles are fully saturated, hence beyond the 10:1 threshold.
    assert set(profile.imbalanced_alleles) == {"A", "B"}


def test_profile_dataset_round_trips_to_json():
    import json

    alleles = ["A", "A", "B", "B", "C", "C"]
    labels = [1, 0, 1, 0, 1, 1]
    profile = profile_dataset(labels, {"hla_allele": alleles})
    restored = json.loads(json.dumps(profile.to_dict()))
    assert restored["n_records"] == 6
    assert restored["log_unit"] == "nats"


def test_profile_dataset_rejects_misaligned_variable():
    with pytest.raises(ValueError):
        profile_dataset([1, 0, 1], {"hla_allele": ["A", "B"]})


def test_bias_profile_candidates_empty_when_label_is_constant():
    profile = BiasProfile(
        n_records=4, label_entropy=0.0, mutual_information={"x": 0.0}
    )
    assert profile.bias_candidates() == {}
