"""Tests for group-aware (peptide-grouped) split construction.

Why this file exists: `build_split` shuffled records inside each allele and
never saw the peptide, so one peptide could land in train and in test at the
same time. A peptide-only control, and any peptide-aware model, can memorise
such a peptide. These tests pin the grouped behaviour and, just as important,
pin that the ungrouped behaviour did not change.

`tests/test_splits.py` stays unchanged and still covers the ungrouped path.
"""

from __future__ import annotations

import json
import random

import pytest

from nbx.bias.diagnostics import log_imbalance_ratio
from nbx.splits.freeze import FrozenSplit, SplitCriterion, build_split

# --------------------------------------------------------------------------
# Corpora
# --------------------------------------------------------------------------


def _ungrouped_corpus(seed: int = 0, n_alleles: int = 12, per_allele: int = 60):
    """Byte-identical to `_corpus` in tests/test_splits.py.

    Copied rather than imported so that the pinned hash below is reproducible
    from this file alone.
    """
    rng = random.Random(seed)
    alleles, labels = [], []
    for i in range(n_alleles):
        allele = f"HLA-A*{i:02d}:01"
        p_pos = 0.05 + 0.9 * (i / max(1, n_alleles - 1))
        for _ in range(per_allele):
            alleles.append(allele)
            labels.append(1 if rng.random() < p_pos else 0)
    return alleles, labels


def _grouped_corpus(
    seed: int = 11,
    n_alleles: int = 8,
    per_allele: int = 120,
    n_peptides: int = 150,
):
    """A corpus where many peptides deliberately occur under several alleles.

    Each allele draws its records from a shared peptide pool, so a peptide is
    very likely to appear under more than one allele -- the situation that
    produced the leakage this module has to prevent.
    """
    rng = random.Random(seed)
    peptides = [f"PEP{i:04d}" for i in range(n_peptides)]
    alleles, labels, groups = [], [], []
    for i in range(n_alleles):
        allele = f"HLA-B*{i:02d}:01"
        p_pos = 0.25 + 0.5 * (i / max(1, n_alleles - 1))
        for _ in range(per_allele):
            alleles.append(allele)
            groups.append(rng.choice(peptides))
            labels.append(1 if rng.random() < p_pos else 0)
    return alleles, labels, groups


def _group_sets(split: FrozenSplit, groups):
    return (
        {groups[i] for i in split.train_idx},
        {groups[i] for i in split.val_idx},
        {groups[i] for i in split.test_idx},
    )


# --------------------------------------------------------------------------
# 1. The ungrouped path must not move
# --------------------------------------------------------------------------
def test_ungrouped_hash_is_unchanged_by_the_grouping_feature():
    """Pinned regression hash for the pre-grouping behaviour.

    Obtained by running, against the unmodified `nbx.splits.freeze` at commit
    487e5b2 (before `groups` existed):

        build_split(_ungrouped_corpus()[0], _ungrouped_corpus()[1]).content_hash

    If this value changes, every hash recorded by an earlier run of the
    ungrouped splitter has been invalidated.
    """
    alleles, labels = _ungrouped_corpus()
    split = build_split(alleles, labels)
    assert split.content_hash == (
        "defe124c86b844c6b3427de48f9c5b24ddf4f01e6c38098167a264916e563802"
    )


def test_ungrouped_seeded_hash_is_unchanged():
    """Second pinned value, criterion seed 7, obtained the same way."""
    alleles, labels = _ungrouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(seed=7))
    assert split.content_hash == (
        "0371f0a1afd4b9b7b056de073d73a3c8da8d3bf6ddcce358bdda44a6e2e8adb4"
    )


def test_groups_none_equals_the_positional_three_argument_call():
    alleles, labels = _ungrouped_corpus()
    crit = SplitCriterion(seed=3)
    a = build_split(alleles, labels, crit)
    b = build_split(alleles, labels, crit, groups=None)
    assert a.content_hash == b.content_hash
    assert a.quarantined_idx == ()


def test_default_criterion_records_no_grouping():
    crit = SplitCriterion()
    assert crit.group_by == "none"
    # "none" is omitted from the canonical form so that hashes recorded before
    # grouping existed stay valid.
    assert "group_by" not in json.loads(crit.canonical_json())


# --------------------------------------------------------------------------
# 2. The grouping invariant
# --------------------------------------------------------------------------
def test_no_group_appears_in_more_than_one_partition():
    alleles, labels, groups = _grouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(), groups=groups)
    train_g, val_g, test_g = _group_sets(split, groups)

    assert not (train_g & test_g), "a peptide is in train and in test"
    assert not (val_g & test_g), "a peptide is in val and in test"
    assert not (train_g & val_g), "a peptide is in train and in val"
    assert test_g, "the grouped split produced no test records at all"


def test_a_record_whose_group_is_in_test_is_never_trained_on():
    """The whole point: no path moves such a record into train or val."""
    alleles, labels, groups = _grouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(), groups=groups)
    test_g = {groups[i] for i in split.test_idx}
    allowed = set(split.test_idx) | set(split.quarantined_idx)
    for i, g in enumerate(groups):
        if g in test_g:
            assert i in allowed, (
                f"record {i} shares peptide {g} with the test set but sits in "
                "train or val"
            )


def test_multi_allele_peptide_is_quarantined_rather_than_leaked():
    """A hand-built corpus: TINY shares every peptide with BAL.

    BAL is large and balanced, so it is admitted and supplies the test block.
    TINY is below the minimum allele size, so it is reported separately -- and
    its records, sharing peptides with the test block, must not become
    training data.
    """
    peptides = [f"P{i:02d}" for i in range(40)]
    alleles = ["BAL"] * 40 + ["TINY"] * 6
    groups = peptides + peptides[:6]
    labels = [i % 2 for i in range(40)] + [1, 0, 1, 0, 1, 0]

    split = build_split(
        alleles, labels, SplitCriterion(min_samples_per_allele=20), groups=groups
    )
    assert "TINY" in split.reported_separately
    test_g = {groups[i] for i in split.test_idx}
    for i, allele in enumerate(alleles):
        if allele == "TINY" and groups[i] in test_g:
            assert i in set(split.quarantined_idx)


# --------------------------------------------------------------------------
# 3. The allele-balance criterion survives grouping
# --------------------------------------------------------------------------
def test_every_admitted_allele_still_meets_the_balance_criterion():
    alleles, labels, groups = _grouped_corpus()
    crit = SplitCriterion(max_abs_rho=1.0)
    split = build_split(alleles, labels, crit, groups=groups)

    counts: dict[str, list[int]] = {}
    for i in split.test_idx:
        c = counts.setdefault(alleles[i], [0, 0])
        c[labels[i]] += 1

    assert counts
    for allele, (n_neg, n_pos) in counts.items():
        rho = log_imbalance_ratio(n_pos, n_neg)
        assert abs(rho) <= crit.max_abs_rho + 1e-9, (
            f"{allele} violates the balance criterion: rho={rho:.3f}"
        )


def test_alleles_contributing_test_records_are_exactly_the_admitted_ones():
    alleles, labels, groups = _grouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(), groups=groups)
    contributing = {alleles[i] for i in split.test_idx}
    assert contributing <= set(split.admitted_alleles)
    assert not (contributing & set(split.reported_separately))


# --------------------------------------------------------------------------
# 4. Nothing is lost, nothing is invented
# --------------------------------------------------------------------------
def test_every_record_is_in_exactly_one_partition_or_quarantined():
    alleles, labels, groups = _grouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(), groups=groups)

    all_idx = (
        list(split.train_idx)
        + list(split.val_idx)
        + list(split.test_idx)
        + list(split.quarantined_idx)
    )
    assert len(all_idx) == len(set(all_idx)), "a record appears twice"
    assert set(all_idx) == set(range(len(alleles))), "records lost or invented"
    assert len(all_idx) == split.n_records


def test_quarantine_reasons_account_for_every_quarantined_record():
    alleles, labels, groups = _grouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(), groups=groups)
    reasons = dict(split.quarantine_reasons)
    assert sum(reasons.values()) == len(split.quarantined_idx)
    assert all(count > 0 for count in reasons.values())


def test_quarantined_records_overlapping_a_partition_are_rejected():
    with pytest.raises(ValueError, match="overlap"):
        FrozenSplit(
            criterion=SplitCriterion(group_by="peptide"),
            train_idx=(0, 1),
            val_idx=(2,),
            test_idx=(3,),
            admitted_alleles=(),
            reported_separately=(),
            n_records=4,
            quarantined_idx=(1,),  # already in train
        )


def test_to_dict_reports_quarantine_counts():
    alleles, labels, groups = _grouped_corpus()
    split = build_split(alleles, labels, SplitCriterion(), groups=groups)
    d = json.loads(json.dumps(split.to_dict()))
    assert d["n_quarantined"] == len(split.quarantined_idx)
    assert d["quarantine_reasons"] == dict(split.quarantine_reasons)
    assert d["criterion"]["group_by"] == "peptide"


# --------------------------------------------------------------------------
# 5. Identity: determinism and hash separation
# --------------------------------------------------------------------------
def test_grouped_split_is_deterministic_for_the_same_seed():
    alleles, labels, groups = _grouped_corpus()
    a = build_split(alleles, labels, SplitCriterion(seed=7), groups=groups)
    b = build_split(alleles, labels, SplitCriterion(seed=7), groups=groups)
    assert a.content_hash == b.content_hash
    assert a.test_idx == b.test_idx
    assert a.quarantined_idx == b.quarantined_idx


def test_grouped_split_does_not_depend_on_input_group_object_type():
    """Determinism must come from the seed, not from set or dict ordering."""
    alleles, labels, groups = _grouped_corpus()
    a = build_split(alleles, labels, SplitCriterion(seed=5), groups=groups)
    b = build_split(alleles, labels, SplitCriterion(seed=5), groups=tuple(groups))
    assert a.content_hash == b.content_hash


def test_changing_the_seed_changes_the_grouped_identity():
    alleles, labels, groups = _grouped_corpus()
    a = build_split(alleles, labels, SplitCriterion(seed=7), groups=groups)
    b = build_split(alleles, labels, SplitCriterion(seed=8), groups=groups)
    assert a.content_hash != b.content_hash


def test_grouped_and_ungrouped_splits_of_the_same_data_differ():
    alleles, labels, groups = _grouped_corpus()
    crit = SplitCriterion(seed=7)
    grouped = build_split(alleles, labels, crit, groups=groups)
    ungrouped = build_split(alleles, labels, crit)
    assert grouped.content_hash != ungrouped.content_hash
    assert grouped.criterion.group_by == "peptide"
    assert ungrouped.criterion.group_by == "none"
    assert json.loads(grouped.criterion.canonical_json())["group_by"] == "peptide"


def test_group_by_label_is_part_of_the_hash_on_its_own():
    """Two splits with identical indices but different grouping labels differ."""
    base = FrozenSplit(
        criterion=SplitCriterion(),
        train_idx=(0, 1),
        val_idx=(2,),
        test_idx=(3,),
        admitted_alleles=(),
        reported_separately=(),
        n_records=4,
    )
    tagged = FrozenSplit(
        criterion=SplitCriterion(group_by="peptide"),
        train_idx=(0, 1),
        val_idx=(2,),
        test_idx=(3,),
        admitted_alleles=(),
        reported_separately=(),
        n_records=4,
    )
    assert base.content_hash != tagged.content_hash


def test_quarantined_records_are_part_of_the_hash():
    common = dict(
        criterion=SplitCriterion(group_by="peptide"),
        train_idx=(0, 1),
        val_idx=(2,),
        test_idx=(3,),
        admitted_alleles=(),
        reported_separately=(),
    )
    a = FrozenSplit(**common, n_records=5, quarantined_idx=(4,))
    b = FrozenSplit(**common, n_records=5, quarantined_idx=())
    assert a.content_hash != b.content_hash


# --------------------------------------------------------------------------
# 6. Input validation
# --------------------------------------------------------------------------
def test_mismatched_group_length_raises():
    with pytest.raises(ValueError, match="parallel"):
        build_split(["A", "A"], [0, 1], groups=["p"])


def test_empty_corpus_with_groups_raises():
    with pytest.raises(ValueError, match="empty"):
        build_split([], [], groups=[])


def test_criterion_claiming_a_grouping_without_groups_raises():
    """A hash may not record a grouping that was not applied."""
    with pytest.raises(ValueError, match="group"):
        build_split(["A"] * 4, [0, 1, 0, 1], SplitCriterion(group_by="peptide"))


def test_criterion_rejects_an_empty_group_by_label():
    with pytest.raises(ValueError):
        SplitCriterion(group_by="")
