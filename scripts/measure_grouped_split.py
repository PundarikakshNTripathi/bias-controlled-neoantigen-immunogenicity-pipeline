#!/usr/bin/env python
"""Measure the peptide-grouped split against the ungrouped one. Freezes nothing.

Builds both splits of the same corpus in memory and prints the counts that
decide whether grouping is worth its cost: how many peptides straddle
partitions without grouping (the leakage), and how many records grouping has
to quarantine to remove it (the price).

Nothing is written. Freezing is a separate, approved step, and
`artefacts/splits/` is immutable.

    python scripts/measure_grouped_split.py
    python scripts/measure_grouped_split.py --dataset artefacts/datasets/....parquet
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from pathlib import Path

import pandas as pd

from nbx.splits.freeze import FrozenSplit, SplitCriterion, build_split

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATASET = ROOT / "artefacts" / "datasets" / "iedb_tcell_majority.parquet"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def partition_of(split: FrozenSplit, n_records: int) -> list[str]:
    """One label per record: train, val, test or quarantined."""
    where = ["unassigned"] * n_records
    for name, idx in (
        ("train", split.train_idx),
        ("val", split.val_idx),
        ("test", split.test_idx),
        ("quarantined", split.quarantined_idx),
    ):
        for i in idx:
            where[i] = name
    return where


def straddling_peptides(split: FrozenSplit, peptides: list[str]) -> tuple[int, int]:
    """(peptides in more than one of train/val/test, records they cover).

    Quarantined records are not in any partition, so they cannot straddle.
    """
    seen: dict[str, set[str]] = {}
    for name, idx in (
        ("train", split.train_idx),
        ("val", split.val_idx),
        ("test", split.test_idx),
    ):
        for i in idx:
            seen.setdefault(peptides[i], set()).add(name)
    multi = {p for p, names in seen.items() if len(names) > 1}
    n_records = sum(1 for p in peptides if p in multi)
    return len(multi), n_records


def label_counts(split_idx: tuple[int, ...], labels: list[int]) -> tuple[int, int]:
    pos = sum(1 for i in split_idx if labels[i] == 1)
    return pos, len(split_idx) - pos


def report(
    title: str,
    split: FrozenSplit,
    alleles: list[str],
    labels: list[int],
    peptides: list[str],
) -> None:
    print(f"\n=== {title} ===")
    print(f"content_hash            {split.content_hash}")
    print(f"short_hash              {split.short_hash}")
    print(f"criterion.group_by      {split.criterion.group_by}")
    print(f"n_records               {split.n_records}")

    n_train, n_val = len(split.train_idx), len(split.val_idx)
    n_test, n_quar = len(split.test_idx), len(split.quarantined_idx)
    total = n_train + n_val + n_test + n_quar
    for name, n in (
        ("train", n_train),
        ("val", n_val),
        ("test", n_test),
        ("quarantined", n_quar),
    ):
        share = 100.0 * n / split.n_records if split.n_records else float("nan")
        print(f"  {name:<22}{n:>8}  ({share:.2f}% of records)")
    print(
        f"  {'accounted total':<22}{total:>8}  "
        f"(== n_records: {total == split.n_records})"
    )

    if split.quarantine_reasons:
        print("quarantine by reason")
        for reason, count in split.quarantine_reasons:
            print(f"  {reason:<50}{count:>8}")
        print(f"  {'sum':<50}{sum(c for _, c in split.quarantine_reasons):>8}")

    print(f"n_admitted_alleles      {len(split.admitted_alleles)}")
    print(f"n_reported_separately   {len(split.reported_separately)}")
    print(f"n_alleles_in_corpus     {len(set(alleles))}")

    pos, neg = label_counts(split.test_idx, labels)
    print(f"test positives          {pos}")
    print(f"test negatives          {neg}")
    for name, idx in (("train", split.train_idx), ("val", split.val_idx)):
        p, n = label_counts(idx, labels)
        print(f"{name} positives / negatives  {p} / {n}")
    if split.quarantined_idx:
        p, n = label_counts(split.quarantined_idx, labels)
        print(f"quarantined positives / negatives  {p} / {n}")
        print(
            "quarantined distinct alleles       "
            f"{len({alleles[i] for i in split.quarantined_idx})}"
        )
        print(
            "quarantined distinct peptides      "
            f"{len({peptides[i] for i in split.quarantined_idx})}"
        )

    n_multi, n_multi_records = straddling_peptides(split, peptides)
    train_p = {peptides[i] for i in split.train_idx}
    val_p = {peptides[i] for i in split.val_idx}
    test_p = {peptides[i] for i in split.test_idx}
    print(f"peptides in >1 partition           {n_multi}")
    print(f"records carrying such peptides     {n_multi_records}")
    print(f"test peptides also in train/val    {len(test_p & (train_p | val_p))}")
    print(
        "distinct peptides per partition    "
        f"train={len(train_p)} val={len(val_p)} test={len(test_p)}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    args = parser.parse_args()

    path = args.dataset if args.dataset.is_absolute() else ROOT / args.dataset
    df = pd.read_parquet(path)

    print(f"dataset                 {path}")
    print(f"file sha256             {file_sha256(path)}")
    print(f"rows x columns          {df.shape[0]} x {df.shape[1]}")

    alleles = df["allele"].astype(str).tolist()
    labels = df["label"].astype(int).tolist()
    peptides = df["peptide"].astype(str).tolist()

    print(f"distinct peptides       {len(set(peptides))}")
    print(f"distinct alleles        {len(set(alleles))}")
    print(f"positives / negatives   {sum(labels)} / {len(labels) - sum(labels)}")
    alleles_of_peptide: dict[str, set[str]] = {}
    for pep, allele in zip(peptides, alleles):
        alleles_of_peptide.setdefault(pep, set()).add(allele)
    multi_allele = {p for p, al in alleles_of_peptide.items() if len(al) > 1}
    print(f"peptides under >1 allele  {len(multi_allele)}")
    print(f"records carrying them     {sum(1 for p in peptides if p in multi_allele)}")
    print(
        "records per peptide       "
        f"max={max(Counter(peptides).values())}"
    )

    criterion = SplitCriterion()  # defaults; nothing tuned here
    print(f"criterion (defaults)    {criterion.canonical_json()}")

    ungrouped = build_split(alleles, labels, criterion)
    grouped = build_split(alleles, labels, criterion, groups=peptides)

    report("UNGROUPED (current behaviour)", ungrouped, alleles, labels, peptides)
    report("PEPTIDE-GROUPED", grouped, alleles, labels, peptides)

    print("\n=== comparison ===")
    u_multi, u_records = straddling_peptides(ungrouped, peptides)
    g_multi, g_records = straddling_peptides(grouped, peptides)
    print(f"peptides straddling partitions   ungrouped={u_multi}  grouped={g_multi}")
    print(f"records involved                 ungrouped={u_records}  grouped={g_records}")
    print(
        "admitted alleles                 "
        f"ungrouped={len(ungrouped.admitted_alleles)}  "
        f"grouped={len(grouped.admitted_alleles)}"
    )
    print(
        "test records                     "
        f"ungrouped={len(ungrouped.test_idx)}  grouped={len(grouped.test_idx)}"
    )
    print(
        "records usable for training      "
        f"ungrouped={len(ungrouped.train_idx)}  grouped={len(grouped.train_idx)}"
    )
    print(
        "records excluded entirely        "
        f"ungrouped={len(ungrouped.quarantined_idx)}  "
        f"grouped={len(grouped.quarantined_idx)}"
    )
    lost = set(ungrouped.admitted_alleles) - set(grouped.admitted_alleles)
    gained = set(grouped.admitted_alleles) - set(ungrouped.admitted_alleles)
    print(f"alleles admitted only ungrouped  {len(lost)}")
    print(f"alleles admitted only grouped    {len(gained)}")
    print("\nNothing was written. No split was frozen.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
