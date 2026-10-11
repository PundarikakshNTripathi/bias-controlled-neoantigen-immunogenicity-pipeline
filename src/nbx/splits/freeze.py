"""Allele-balanced split construction and freezing (requirement R1).

Two properties matter here and are both tested:

1.  Every allele admitted to the test set carries a bounded intra-allelic
    positive-to-negative ratio, so a model cannot succeed on it by
    memorising allele-level label frequencies.
2.  A split is identified by the content hash of its definition. Change the
    criterion and you get a different hash, and therefore a different
    experiment. There is no way to mutate a frozen split in place.

Alleles that cannot meet the balance criterion are recorded in
`reported_separately` and are never silently pooled into the headline
metrics (requirement FR-61).

Grouped splits (analysis plan section 3)
----------------------------------------

Records are rows of (peptide, allele). One peptide can occur under several
alleles, so a splitter that shuffles records inside each allele -- which is
what the ungrouped path below does -- can put the same peptide in train and
in test. A peptide-only control and any peptide-aware model can memorise such
a peptide, which is leakage.

Passing `groups` (the peptide of each record) switches on the grouped
algorithm. Its invariant: **no group value appears in more than one of train,
val and test.** No record is silently moved to satisfy that invariant; a
record that cannot be placed without breaking it is quarantined, that is,
excluded from every partition and counted under a reason.

The algorithm, in order, with `rng = random.Random(criterion.seed)` as the
only source of randomness:

1.  Group the records by allele, and inside each allele by group value. The
    unit of selection is a whole (allele, group) cell: all records of one
    group under one allele move together.
2.  For each allele in `sorted(..., key=str)` order:
    a.  Skip the allele (report it separately) if it has fewer than
        `min_samples_per_allele` records, or only one class.
    b.  Compute the same target test block as the ungrouped path: a block of
        about `test_fraction` of the allele's records, trimmed one record at a
        time until `|rho| <= max_abs_rho`.
    c.  Sort the allele's units by `str(group)` for a deterministic starting
        order, shuffle them with `rng`, then first-fit: accept a unit while it
        keeps the running positive and negative counts within the target. Unit
        order is independent of group multiplicity on purpose -- preferring
        single-allele peptides would bias the test set towards peptides that
        were only ever assayed against one allele.
    d.  If the realised block still violates `|rho| <= max_abs_rho`, drop
        whole units, each time the one contributing most to the majority side
        (ties broken by position in the shuffled order), until the block
        satisfies the criterion or is empty. An allele whose block ends up
        empty or single-class is reported separately and contributes no test
        records.
3.  Let `test_groups` be the group values of the selected records. Every
    record that is not itself selected but whose group is in `test_groups` is
    quarantined, under one of two reasons: its allele was admitted but the
    record was not in that allele's balanced block, or its allele was reported
    separately. Those records are the multi-allele peptides; putting them in
    train or val is exactly the leakage being prevented.
4.  Train and validation are assigned from the remaining records at group
    level: the remaining group values are sorted by `str`, shuffled with the
    same `rng`, and taken into train until train reaches
    `round(n_remaining * train_fraction / (train_fraction + val_fraction))`
    records; the rest go to validation. Because whole groups move, the
    realised fractions deviate from the nominal ones by up to the size of one
    group.

Hash compatibility: `SplitCriterion.canonical_json` omits `group_by` when it
is `"none"`, and `FrozenSplit.compute_hash` omits the quarantine list when it
is empty. So an ungrouped split keeps the hash it had before grouping existed
("none" and absent mean the same thing), while any grouped split -- which
always records `group_by` -- hashes differently from the ungrouped split of
the same data. The quarantine reason counts are a derived report and are not
part of the identity; the quarantined indices are.
"""

from __future__ import annotations

import hashlib
import json
import random
from collections.abc import Hashable, Sequence
from dataclasses import asdict, dataclass, field, replace

from nbx.bias.diagnostics import log_imbalance_ratio

__all__ = [
    "SplitCriterion",
    "FrozenSplit",
    "build_split",
    "NO_GROUPING",
    "DEFAULT_GROUPING",
    "QUARANTINE_ADMITTED_NOT_SELECTED",
    "QUARANTINE_ALLELE_REPORTED_SEPARATELY",
]

NO_GROUPING = "none"
"""`SplitCriterion.group_by` value meaning the records were not grouped."""

DEFAULT_GROUPING = "peptide"
"""Recorded when `groups` is supplied and the criterion names no grouping."""

QUARANTINE_ADMITTED_NOT_SELECTED = "group_in_test_not_selected_for_admitted_allele"
QUARANTINE_ALLELE_REPORTED_SEPARATELY = "group_in_test_allele_reported_separately"


@dataclass(frozen=True)
class SplitCriterion:
    """The rule set that defines a split. Hashed into the split identity."""

    max_abs_rho: float = 1.0
    """Maximum |log2(N+/N-)| an allele may have in the test set. 1.0 == 2:1."""

    min_samples_per_allele: int = 20
    """Alleles with fewer samples than this are reported separately."""

    train_fraction: float = 0.70
    val_fraction: float = 0.15
    test_fraction: float = 0.15

    seed: int = 20260928

    group_by: str = NO_GROUPING
    """What a partition-exclusive unit is: "none", or e.g. "peptide".

    Part of the split identity whenever it is not "none" (see the module
    docstring on hash compatibility). It records the grouping that was applied,
    so a grouped and an ungrouped split of the same data cannot be confused.
    """

    def __post_init__(self) -> None:
        total = self.train_fraction + self.val_fraction + self.test_fraction
        if not (0.999 < total < 1.001):
            raise ValueError(f"split fractions must sum to 1.0; got {total}")
        if self.max_abs_rho < 0:
            raise ValueError("max_abs_rho must be non-negative")
        if self.min_samples_per_allele < 2:
            raise ValueError(
                "min_samples_per_allele must be at least 2; a single-sample "
                "allele has no measurable imbalance"
            )
        if not isinstance(self.group_by, str) or not self.group_by.strip():
            raise ValueError(
                "group_by must be a non-empty label; use 'none' for no grouping"
            )

    @property
    def is_grouped(self) -> bool:
        return self.group_by != NO_GROUPING

    def canonical_json(self) -> str:
        payload = asdict(self)
        if payload.get("group_by") == NO_GROUPING:
            # Omitted, not written as "none", so that hashes recorded before
            # grouping existed remain valid. Absent and "none" are the same.
            payload.pop("group_by")
        return json.dumps(payload, sort_keys=True, separators=(",", ":"))


@dataclass(frozen=True)
class FrozenSplit:
    """An immutable split definition, identified by its content hash."""

    criterion: SplitCriterion
    train_idx: tuple[int, ...]
    val_idx: tuple[int, ...]
    test_idx: tuple[int, ...]
    admitted_alleles: tuple[str, ...]
    reported_separately: tuple[str, ...]
    n_records: int

    quarantined_idx: tuple[int, ...] = ()
    """Records excluded from every partition to keep groups partition-exclusive.

    Empty for an ungrouped split. These records are not training data, not
    validation data and not test data; they are reported, never used.
    """

    quarantine_reasons: tuple[tuple[str, int], ...] = ()
    """(reason, count) pairs summing to `len(quarantined_idx)`.

    A derived report: not part of the content hash, unlike the indices.
    """

    content_hash: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        held = set(self.quarantined_idx)
        overlap = (
            (set(self.train_idx) & set(self.val_idx))
            | (set(self.train_idx) & set(self.test_idx))
            | (set(self.val_idx) & set(self.test_idx))
            | (held & set(self.train_idx))
            | (held & set(self.val_idx))
            | (held & set(self.test_idx))
        )
        if overlap:
            raise ValueError(
                f"split partitions overlap on {len(overlap)} record(s); "
                "leakage between train and test invalidates every result"
            )
        if not self.content_hash:
            object.__setattr__(self, "content_hash", self.compute_hash())

    def compute_hash(self) -> str:
        definition: dict[str, object] = {
            "criterion": json.loads(self.criterion.canonical_json()),
            "train": list(self.train_idx),
            "val": list(self.val_idx),
            "test": list(self.test_idx),
            "admitted": list(self.admitted_alleles),
            "separate": list(self.reported_separately),
            "n_records": self.n_records,
        }
        if self.quarantined_idx:
            # Omitted when empty so that an ungrouped split keeps the hash it
            # had before quarantine existed.
            definition["quarantined"] = list(self.quarantined_idx)
        payload = json.dumps(definition, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @property
    def short_hash(self) -> str:
        return f"sha256:{self.content_hash[:6]}..{self.content_hash[-3:]}"

    def to_dict(self) -> dict:
        return {
            "content_hash": self.content_hash,
            "criterion": json.loads(self.criterion.canonical_json()),
            "n_records": self.n_records,
            "n_train": len(self.train_idx),
            "n_val": len(self.val_idx),
            "n_test": len(self.test_idx),
            "n_quarantined": len(self.quarantined_idx),
            "quarantine_reasons": dict(self.quarantine_reasons),
            "n_admitted_alleles": len(self.admitted_alleles),
            "n_reported_separately": len(self.reported_separately),
            "admitted_alleles": list(self.admitted_alleles),
            "reported_separately": list(self.reported_separately),
        }


def build_split(
    alleles: Sequence[Hashable],
    labels: Sequence[int],
    criterion: SplitCriterion | None = None,
    groups: Sequence[Hashable] | None = None,
) -> FrozenSplit:
    """Construct an allele-balanced, frozen split.

    The test partition is built allele by allele. For each allele that has
    enough samples, positives and negatives are drawn in a ratio that keeps
    |rho| within the criterion. Alleles that cannot meet the criterion
    contribute no test records and are listed in `reported_separately`.

    `groups` is an optional parallel sequence naming the unit that must not
    straddle partitions -- in this project the peptide. When it is given, no
    group value appears in more than one of train, val and test, and records
    that would break that are quarantined rather than moved; see the module
    docstring for the algorithm. When it is `None` the behaviour, and the
    content hash, are those of the ungrouped splitter.

    Determinism: the only source of randomness is `criterion.seed`, so the
    same inputs and the same criterion always produce the same hash.
    """
    criterion = criterion or SplitCriterion()

    als = list(alleles)
    ys = [int(y) for y in labels]
    if len(als) != len(ys):
        raise ValueError("alleles and labels must be parallel sequences")
    if not als:
        raise ValueError("cannot build a split from an empty corpus")
    if any(y not in (0, 1) for y in ys):
        raise ValueError("labels must be binary 0/1")

    if groups is None:
        if criterion.is_grouped:
            raise ValueError(
                f"criterion records group_by={criterion.group_by!r} but no "
                "groups were supplied; the hash would claim a grouping that "
                "was never applied"
            )
    else:
        gs = list(groups)
        if len(gs) != len(als):
            raise ValueError("groups must be a parallel sequence of alleles")
        if not criterion.is_grouped:
            criterion = replace(criterion, group_by=DEFAULT_GROUPING)
        return _build_grouped(als, ys, gs, criterion)

    rng = random.Random(criterion.seed)

    by_allele: dict[Hashable, dict[int, list[int]]] = {}
    for i, (a, y) in enumerate(zip(als, ys)):
        by_allele.setdefault(a, {0: [], 1: []})[y].append(i)

    admitted: list[str] = []
    separate: list[str] = []
    test_idx: list[int] = []
    remaining: list[int] = []

    for allele in sorted(by_allele, key=str):
        pos = list(by_allele[allele][1])
        neg = list(by_allele[allele][0])
        n_total = len(pos) + len(neg)

        if n_total < criterion.min_samples_per_allele or not pos or not neg:
            separate.append(str(allele))
            remaining.extend(pos + neg)
            continue

        rng.shuffle(pos)
        rng.shuffle(neg)

        # Largest balanced test block this allele can supply.
        take_pos, take_neg = _balanced_target(
            len(pos), len(neg), n_total, criterion
        )

        if take_pos == 0 or take_neg == 0:
            separate.append(str(allele))
            remaining.extend(pos + neg)
            continue

        admitted.append(str(allele))
        test_idx.extend(pos[:take_pos] + neg[:take_neg])
        remaining.extend(pos[take_pos:] + neg[take_neg:])

    # Train / validation from everything the test set did not claim.
    rng.shuffle(remaining)
    denom = criterion.train_fraction + criterion.val_fraction
    n_train = int(round(len(remaining) * (criterion.train_fraction / denom)))
    train_idx = sorted(remaining[:n_train])
    val_idx = sorted(remaining[n_train:])

    return FrozenSplit(
        criterion=criterion,
        train_idx=tuple(train_idx),
        val_idx=tuple(val_idx),
        test_idx=tuple(sorted(test_idx)),
        admitted_alleles=tuple(admitted),
        reported_separately=tuple(separate),
        n_records=len(als),
    )


def _balanced_target(
    n_pos: int, n_neg: int, n_total: int, criterion: SplitCriterion
) -> tuple[int, int]:
    """How many positives and negatives one allele may put in the test block.

    About `test_fraction` of the allele's records, split as evenly as the
    available classes allow, then trimmed one record at a time off the
    majority side until |rho| is inside `max_abs_rho`. Returns (0, 0) when no
    such block exists, which means the allele is reported separately.
    """
    want = max(1, int(round(n_total * criterion.test_fraction)))
    take_pos = min(n_pos, want // 2 + want % 2)
    take_neg = min(n_neg, want - take_pos)
    take_pos = min(take_pos, n_pos)

    while take_pos > 0 and take_neg > 0 and abs(
        log_imbalance_ratio(take_pos, take_neg)
    ) > criterion.max_abs_rho:
        if take_pos > take_neg:
            take_pos -= 1
        else:
            take_neg -= 1

    if take_pos == 0 or take_neg == 0:
        return 0, 0
    return take_pos, take_neg


def _build_grouped(
    als: list[Hashable],
    ys: list[int],
    gs: list[Hashable],
    criterion: SplitCriterion,
) -> FrozenSplit:
    """Group-exclusive split construction. See the module docstring."""
    rng = random.Random(criterion.seed)

    # allele -> group -> record indices. The unit of selection is one cell.
    cells: dict[Hashable, dict[Hashable, list[int]]] = {}
    for i, (a, g) in enumerate(zip(als, gs)):
        cells.setdefault(a, {}).setdefault(g, []).append(i)

    admitted: list[str] = []
    separate: list[str] = []
    test_idx: list[int] = []

    for allele in sorted(cells, key=str):
        by_group = cells[allele]
        n_pos = sum(1 for idxs in by_group.values() for i in idxs if ys[i] == 1)
        n_records = sum(len(idxs) for idxs in by_group.values())
        n_neg = n_records - n_pos

        if n_records < criterion.min_samples_per_allele or not n_pos or not n_neg:
            separate.append(str(allele))
            continue

        target_pos, target_neg = _balanced_target(n_pos, n_neg, n_records, criterion)
        if target_pos == 0 or target_neg == 0:
            separate.append(str(allele))
            continue

        # Deterministic starting order, then one shuffle from the seed.
        units: list[tuple[Hashable, list[int], int, int]] = []
        for group in sorted(by_group, key=str):
            idxs = by_group[group]
            up = sum(1 for i in idxs if ys[i] == 1)
            units.append((group, idxs, up, len(idxs) - up))
        rng.shuffle(units)

        # First fit: accept a whole unit while it stays inside the target.
        selected: list[tuple[Hashable, list[int], int, int]] = []
        got_pos = got_neg = 0
        for unit in units:
            if got_pos + unit[2] <= target_pos and got_neg + unit[3] <= target_neg:
                selected.append(unit)
                got_pos += unit[2]
                got_neg += unit[3]

        # Whole units leave the realised block coarser than the target, so the
        # criterion is enforced again here, by dropping whole units.
        while selected and (
            got_pos == 0
            or got_neg == 0
            or abs(log_imbalance_ratio(got_pos, got_neg)) > criterion.max_abs_rho
        ):
            if got_pos >= got_neg:
                excess = [u[2] - u[3] for u in selected]
            else:
                excess = [u[3] - u[2] for u in selected]
            worst = max(range(len(selected)), key=lambda j: (excess[j], -j))
            dropped = selected.pop(worst)
            got_pos -= dropped[2]
            got_neg -= dropped[3]

        if not selected:
            separate.append(str(allele))
            continue

        admitted.append(str(allele))
        for unit in selected:
            test_idx.extend(unit[1])

    # Any record sharing a group with the test set is withheld from training.
    in_test = set(test_idx)
    test_groups = {gs[i] for i in in_test}
    admitted_set = set(admitted)

    quarantined: list[int] = []
    reasons: dict[str, int] = {}
    remaining: list[int] = []
    for i in range(len(als)):
        if i in in_test:
            continue
        if gs[i] in test_groups:
            reason = (
                QUARANTINE_ADMITTED_NOT_SELECTED
                if str(als[i]) in admitted_set
                else QUARANTINE_ALLELE_REPORTED_SEPARATELY
            )
            reasons[reason] = reasons.get(reason, 0) + 1
            quarantined.append(i)
        else:
            remaining.append(i)

    # Train and validation are assigned whole groups at a time.
    members: dict[Hashable, list[int]] = {}
    for i in remaining:
        members.setdefault(gs[i], []).append(i)
    order = sorted(members, key=str)
    rng.shuffle(order)

    denom = criterion.train_fraction + criterion.val_fraction
    n_train_target = int(round(len(remaining) * (criterion.train_fraction / denom)))
    train_idx: list[int] = []
    val_idx: list[int] = []
    for group in order:
        if len(train_idx) < n_train_target:
            train_idx.extend(members[group])
        else:
            val_idx.extend(members[group])

    return FrozenSplit(
        criterion=criterion,
        train_idx=tuple(sorted(train_idx)),
        val_idx=tuple(sorted(val_idx)),
        test_idx=tuple(sorted(test_idx)),
        admitted_alleles=tuple(admitted),
        reported_separately=tuple(separate),
        n_records=len(als),
        quarantined_idx=tuple(sorted(quarantined)),
        quarantine_reasons=tuple(sorted(reasons.items())),
    )
