"""Tests for evaluation metrics and the adaptive masking policy."""

from __future__ import annotations

import pytest

from nbx.bias.policy import (
    adaptive_policy,
    per_position_mi,
    select_threshold,
    uniform_policy,
)
from nbx.evaluate.metrics import (
    allele_stratified,
    auprc,
    auroc,
    fraction_ranked,
    matthews_corrcoef,
    precision_at_k,
    recall_at_k,
    top_k_immunogenic_fraction,
)


# ==========================================================================
# Metrics
# ==========================================================================
def test_auprc_perfect_ranking_is_one():
    assert auprc([0.9, 0.8, 0.2, 0.1], [1, 1, 0, 0]) == pytest.approx(1.0)


def test_auprc_matches_sklearn_average_precision():
    sk = pytest.importorskip("sklearn.metrics")
    scores = [0.11, 0.92, 0.37, 0.55, 0.78, 0.04, 0.61, 0.29]
    labels = [0, 1, 0, 1, 1, 0, 0, 1]
    assert auprc(scores, labels) == pytest.approx(
        sk.average_precision_score(labels, scores), rel=1e-9
    )


def test_auroc_matches_sklearn_including_ties():
    sk = pytest.importorskip("sklearn.metrics")
    scores = [0.5, 0.5, 0.5, 0.9, 0.1, 0.7, 0.7, 0.2]
    labels = [1, 0, 1, 1, 0, 0, 1, 0]
    assert auroc(scores, labels) == pytest.approx(
        sk.roc_auc_score(labels, scores), rel=1e-9
    )


def test_auprc_undefined_without_positives():
    with pytest.raises(ValueError):
        auprc([0.4, 0.6], [0, 0])


def test_auroc_undefined_when_one_class_absent():
    with pytest.raises(ValueError):
        auroc([0.4, 0.6], [1, 1])


def test_ties_are_broken_pessimistically():
    """Equal scores must not be ordered so as to flatter the model."""
    # All scores identical: the single positive is placed last.
    assert precision_at_k([0.5] * 4, [1, 0, 0, 0], k=1) == 0.0


def test_precision_and_recall_at_k():
    scores = [0.9, 0.8, 0.7, 0.6, 0.5]
    labels = [1, 0, 1, 0, 0]
    assert precision_at_k(scores, labels, 2) == pytest.approx(0.5)
    assert recall_at_k(scores, labels, 3) == pytest.approx(1.0)


def test_precision_at_k_rejects_k_beyond_candidate_list():
    with pytest.raises(ValueError):
        precision_at_k([0.5, 0.4], [1, 0], k=5)


def test_tesla_metrics_are_the_documented_quantities():
    scores = [1.0 - i / 200 for i in range(200)]
    labels = [1 if i < 50 else 0 for i in range(200)]
    # All 50 positives sit in the top 100, so fraction ranked is 1.0.
    assert fraction_ranked(scores, labels) == pytest.approx(1.0)
    # And the top 20 are all positive.
    assert top_k_immunogenic_fraction(scores, labels) == pytest.approx(1.0)


def test_matthews_matches_sklearn():
    sk = pytest.importorskip("sklearn.metrics")
    preds = [1, 0, 1, 1, 0, 0, 1, 0]
    labels = [1, 0, 0, 1, 0, 1, 1, 0]
    assert matthews_corrcoef(preds, labels) == pytest.approx(
        sk.matthews_corrcoef(labels, preds), rel=1e-9
    )


def test_matthews_is_zero_for_a_constant_predictor():
    assert matthews_corrcoef([1, 1, 1, 1], [1, 0, 1, 0]) == 0.0


# --------------------------------------------------------------------------
# Allele stratification (R3)
# --------------------------------------------------------------------------
def test_allele_stratified_reports_per_allele_and_weighted_aggregate():
    scores = [0.9, 0.1, 0.8, 0.2, 0.7, 0.3]
    labels = [1, 0, 1, 0, 1, 0]
    alleles = ["A", "A", "B", "B", "C", "C"]
    per_allele, aggregate = allele_stratified(auprc, scores, labels, alleles)
    assert set(per_allele) == {"A", "B", "C"}
    assert aggregate == pytest.approx(1.0)


def test_allele_stratified_omits_undefined_alleles_rather_than_zeroing_them():
    # Allele "C" has no positives, so AUPRC is undefined for it.
    scores = [0.9, 0.1, 0.8, 0.2, 0.7, 0.3]
    labels = [1, 0, 1, 0, 0, 0]
    alleles = ["A", "A", "B", "B", "C", "C"]
    per_allele, _ = allele_stratified(auprc, scores, labels, alleles)
    assert "C" not in per_allele, "an undefined metric must be omitted, not zeroed"
    assert set(per_allele) == {"A", "B"}


def test_allele_stratified_raises_when_nothing_is_computable():
    with pytest.raises(ValueError):
        allele_stratified(auprc, [0.5, 0.4], [0, 0], ["A", "A"])


def test_unweighted_aggregate_differs_from_weighted_when_support_is_uneven():
    scores = [0.9, 0.1] + [0.1, 0.9] * 10
    labels = [1, 0] + [1, 0] * 10
    alleles = ["SMALL", "SMALL"] + ["BIG"] * 20
    _, weighted = allele_stratified(auprc, scores, labels, alleles, True)
    _, unweighted = allele_stratified(auprc, scores, labels, alleles, False)
    assert weighted != pytest.approx(unweighted)


# ==========================================================================
# Masking policy
# ==========================================================================
def test_uniform_policy_masks_position_two_and_the_c_terminus():
    policy = uniform_policy(peptide_length=9, alleles=["A", "B"])
    assert policy.mask_for("A") == frozenset({1, 8})
    assert policy.mask_for("B") == frozenset({1, 8})
    assert policy.is_uniform()


def test_uniform_policy_covers_an_unseen_allele():
    policy = uniform_policy(9, alleles=["A"])
    assert policy.mask_for("NEVER_SEEN") == frozenset({1, 8})


def test_adaptive_policy_differs_per_allele_when_anchors_differ():
    mi = {
        "A": [0.01, 0.90, 0.02, 0.01, 0.02, 0.03, 0.01, 0.02, 0.88],
        "B": [0.01, 0.85, 0.80, 0.01, 0.02, 0.03, 0.01, 0.02, 0.90],
    }
    policy = adaptive_policy(mi, tau=0.5, peptide_length=9)
    assert policy.mask_for("A") == frozenset({1, 8})
    assert policy.mask_for("B") == frozenset({1, 2, 8})
    assert not policy.is_uniform()


def test_adaptive_policy_reduces_to_uniform_when_anchor_sets_match():
    """The clean-ablation property stated in Section 7.6.2 and FR-52."""
    shared = [0.01, 0.90, 0.02, 0.01, 0.02, 0.03, 0.01, 0.02, 0.88]
    mi = {"A": list(shared), "B": list(shared), "C": list(shared)}
    adaptive = adaptive_policy(mi, tau=0.5, peptide_length=9)
    uniform = uniform_policy(9, alleles=["A", "B", "C"])

    assert adaptive.is_uniform()
    for allele in ("A", "B", "C"):
        assert adaptive.mask_for(allele) == uniform.mask_for(allele)


def test_adaptive_policy_can_be_restricted_to_candidate_anchor_positions():
    mi = {"A": [0.9] * 9}
    policy = adaptive_policy(
        mi, tau=0.5, peptide_length=9, restrict_to=frozenset({1, 8})
    )
    assert policy.mask_for("A") == frozenset({1, 8}), (
        "positions outside the candidate anchor set must never be masked"
    )


def test_adaptive_policy_rejects_wrong_length_mi_vector():
    with pytest.raises(ValueError):
        adaptive_policy({"A": [0.1, 0.2]}, tau=0.05, peptide_length=9)


def test_per_position_mi_finds_the_informative_position():
    # Position 1 determines the label; every other position is noise.
    peptides, labels, alleles = [], [], []
    for i in range(40):
        y = i % 2
        peptides.append(["X", "K" if y else "A", "X"])
        labels.append(y)
        alleles.append("A*02:01")

    mi = per_position_mi(peptides, labels, alleles, peptide_length=3)
    values = mi["A*02:01"]
    assert values[1] > values[0]
    assert values[1] > values[2]


def test_per_position_mi_is_zero_for_a_constant_label_allele():
    peptides = [["A", "B", "C"]] * 5
    mi = per_position_mi(peptides, [1] * 5, ["X"] * 5, peptide_length=3)
    assert mi["X"] == [0.0, 0.0, 0.0]


def test_select_threshold_is_a_quantile_of_pooled_values():
    mi = {"A": [0.0, 0.1, 0.2], "B": [0.3, 0.4, 0.5]}
    assert select_threshold(mi, quantile=0.5) == pytest.approx(0.2)


def test_select_threshold_rejects_degenerate_quantile():
    with pytest.raises(ValueError):
        select_threshold({"A": [0.1]}, quantile=1.0)


def test_policy_serialises_to_json():
    import json

    policy = uniform_policy(9, ["A", "B"])
    restored = json.loads(json.dumps(policy.to_dict()))
    assert restored["is_uniform"] is True
    assert restored["default"] == [1, 8]
