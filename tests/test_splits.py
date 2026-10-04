"""Tests for allele-balanced split construction and freezing (R1)."""

from __future__ import annotations

import random

import pytest

from nbx.bias.diagnostics import log_imbalance_ratio
from nbx.splits.freeze import FrozenSplit, SplitCriterion, build_split


def _corpus(seed: int = 0, n_alleles: int = 12, per_allele: int = 60):
    """A synthetic corpus with deliberately varied per-allele imbalance."""
    rng = random.Random(seed)
    alleles, labels = [], []
    for i in range(n_alleles):
        allele = f"HLA-A*{i:02d}:01"
        p_pos = 0.05 + 0.9 * (i / max(1, n_alleles - 1))
        for _ in range(per_allele):
            alleles.append(allele)
            labels.append(1 if rng.random() < p_pos else 0)
    return alleles, labels


# --------------------------------------------------------------------------
# The balance criterion
# --------------------------------------------------------------------------
def test_every_admitted_allele_meets_the_balance_criterion():
    alleles, labels = _corpus()
    crit = SplitCriterion(max_abs_rho=1.0)
    split = build_split(alleles, labels, crit)

    counts: dict[str, list[int]] = {}
    for i in split.test_idx:
        c = counts.setdefault(alleles[i], [0, 0])
        c[labels[i]] += 1

    assert counts, "the split produced no test records at all"
    for allele, (n_neg, n_pos) in counts.items():
        rho = log_imbalance_ratio(n_pos, n_neg)
        assert abs(rho) <= crit.max_abs_rho + 1e-9, (
            f"{allele} violates the balance criterion: rho={rho:.3f}"
        )


def test_alleles_failing_the_criterion_are_reported_not_pooled():
    # One allele is entirely positive and cannot be balanced at any ratio.
    alleles = ["BAL"] * 40 + ["SKEW"] * 40
    labels = [i % 2 for i in range(40)] + [1] * 40

    split = build_split(alleles, labels, SplitCriterion(min_samples_per_allele=10))

    assert "SKEW" in split.reported_separately
    assert "SKEW" not in split.admitted_alleles
    assert all(alleles[i] != "SKEW" for i in split.test_idx), (
        "an allele that cannot be balanced must contribute no test records"
    )


def test_small_alleles_are_reported_separately():
    alleles = ["BIG"] * 40 + ["TINY"] * 3
    labels = [i % 2 for i in range(40)] + [1, 0, 1]
    split = build_split(alleles, labels, SplitCriterion(min_samples_per_allele=20))
    assert "TINY" in split.reported_separately


# --------------------------------------------------------------------------
# Partition integrity
# --------------------------------------------------------------------------
def test_partitions_are_disjoint_and_cover_the_corpus():
    alleles, labels = _corpus()
    split = build_split(alleles, labels)

    all_idx = list(split.train_idx) + list(split.val_idx) + list(split.test_idx)
    assert len(all_idx) == len(set(all_idx)), "a record appears in two partitions"
    assert set(all_idx) == set(range(len(alleles))), "records were lost or invented"


def test_overlapping_partitions_are_rejected_at_construction():
    with pytest.raises(ValueError, match="overlap"):
        FrozenSplit(
            criterion=SplitCriterion(),
            train_idx=(0, 1, 2),
            val_idx=(2, 3),           # 2 is in both
            test_idx=(4,),
            admitted_alleles=(),
            reported_separately=(),
            n_records=5,
        )


# --------------------------------------------------------------------------
# Freezing and identity
# --------------------------------------------------------------------------
def test_split_is_deterministic_for_the_same_seed():
    alleles, labels = _corpus()
    a = build_split(alleles, labels, SplitCriterion(seed=7))
    b = build_split(alleles, labels, SplitCriterion(seed=7))
    assert a.content_hash == b.content_hash
    assert a.test_idx == b.test_idx


def test_changing_the_seed_changes_the_identity():
    alleles, labels = _corpus()
    a = build_split(alleles, labels, SplitCriterion(seed=7))
    b = build_split(alleles, labels, SplitCriterion(seed=8))
    assert a.content_hash != b.content_hash


def test_changing_the_criterion_changes_the_identity():
    """A different rule is a different experiment, not an edit of this one."""
    alleles, labels = _corpus()
    a = build_split(alleles, labels, SplitCriterion(max_abs_rho=1.0))
    b = build_split(alleles, labels, SplitCriterion(max_abs_rho=0.5))
    assert a.content_hash != b.content_hash


def test_frozen_split_is_immutable():
    alleles, labels = _corpus()
    split = build_split(alleles, labels)
    with pytest.raises(AttributeError):
        split.test_idx = (1, 2, 3)  # type: ignore[misc]


def test_short_hash_is_stable_and_prefixed():
    alleles, labels = _corpus()
    split = build_split(alleles, labels)
    assert split.short_hash.startswith("sha256:")
    assert split.short_hash == split.short_hash


# --------------------------------------------------------------------------
# Criterion validation
# --------------------------------------------------------------------------
def test_criterion_rejects_fractions_that_do_not_sum_to_one():
    with pytest.raises(ValueError):
        SplitCriterion(train_fraction=0.7, val_fraction=0.2, test_fraction=0.2)


def test_criterion_rejects_single_sample_minimum():
    with pytest.raises(ValueError):
        SplitCriterion(min_samples_per_allele=1)


def test_build_split_rejects_empty_corpus():
    with pytest.raises(ValueError):
        build_split([], [])


def test_build_split_rejects_non_binary_labels():
    with pytest.raises(ValueError):
        build_split(["A", "A"], [0, 2])
