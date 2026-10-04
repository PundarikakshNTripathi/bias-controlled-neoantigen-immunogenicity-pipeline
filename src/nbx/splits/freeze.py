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
"""

from __future__ import annotations

import hashlib
import json
import random
from collections.abc import Hashable, Sequence
from dataclasses import asdict, dataclass, field

from nbx.bias.diagnostics import log_imbalance_ratio

__all__ = ["SplitCriterion", "FrozenSplit", "build_split"]


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

    def canonical_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))


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
    content_hash: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        overlap = (
            (set(self.train_idx) & set(self.val_idx))
            | (set(self.train_idx) & set(self.test_idx))
            | (set(self.val_idx) & set(self.test_idx))
        )
        if overlap:
            raise ValueError(
                f"split partitions overlap on {len(overlap)} record(s); "
                "leakage between train and test invalidates every result"
            )
        if not self.content_hash:
            object.__setattr__(self, "content_hash", self.compute_hash())

    def compute_hash(self) -> str:
        payload = json.dumps(
            {
                "criterion": json.loads(self.criterion.canonical_json()),
                "train": list(self.train_idx),
                "val": list(self.val_idx),
                "test": list(self.test_idx),
                "admitted": list(self.admitted_alleles),
                "separate": list(self.reported_separately),
                "n_records": self.n_records,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
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
            "n_admitted_alleles": len(self.admitted_alleles),
            "n_reported_separately": len(self.reported_separately),
            "admitted_alleles": list(self.admitted_alleles),
            "reported_separately": list(self.reported_separately),
        }


def build_split(
    alleles: Sequence[Hashable],
    labels: Sequence[int],
    criterion: SplitCriterion | None = None,
) -> FrozenSplit:
    """Construct an allele-balanced, frozen split.

    The test partition is built allele by allele. For each allele that has
    enough samples, positives and negatives are drawn in a ratio that keeps
    |rho| within the criterion. Alleles that cannot meet the criterion
    contribute no test records and are listed in `reported_separately`.

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
        want = max(1, int(round(n_total * criterion.test_fraction)))
        take_pos = min(len(pos), want // 2 + want % 2)
        take_neg = min(len(neg), want - take_pos)
        take_pos = min(take_pos, len(pos))

        # Trim the majority side until the ratio is inside the criterion.
        while take_pos > 0 and take_neg > 0 and abs(
            log_imbalance_ratio(take_pos, take_neg)
        ) > criterion.max_abs_rho:
            if take_pos > take_neg:
                take_pos -= 1
            else:
                take_neg -= 1

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
