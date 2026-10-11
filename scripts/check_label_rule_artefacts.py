#!/usr/bin/env python
"""Compare the three label-rule artefacts, and guard the `strict` regression.

Reads the parquet files under `artefacts/datasets/` and recomputes every
number from the frames. Nothing here is a constant: if a count is printed, it
was just read off a frame.

Two things are checked rather than merely printed:

1.  `strict` must assert exactly the same (peptide, allele, label) set as the
    pre-existing `iedb_tcell.parquet`, which was built before the label rules
    existed with the drop-every-conflicting-pair policy. If that set changed,
    the pair refactor changed a label somewhere and the exit code is 1.
2.  Each artefact must still hash to the value its own report records.

Usage:
    python scripts/check_label_rule_artefacts.py
    python scripts/check_label_rule_artefacts.py --datasets artefacts/datasets
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from nbx.data.harmonise import (  # noqa: E402
    LABEL_RULES,
    content_hash,
    default_artefact_name,
)

BASELINE_NAME = "iedb_tcell"
"""The artefact built before the label rules existed. The regression target."""


def pair_label_set(frame: pd.DataFrame) -> set[tuple[str, str, int]]:
    return set(
        zip(
            frame["peptide"].astype(str),
            frame["allele"].astype(str),
            frame["label"].astype(int),
        )
    )


def load(datasets: Path, name: str) -> tuple[pd.DataFrame, dict]:
    frame = pd.read_parquet(datasets / f"{name}.parquet")
    report = json.loads((datasets / f"{name}.report.json").read_text(encoding="utf-8"))
    return frame, report


def describe(name: str, frame: pd.DataFrame, report: dict) -> None:
    n_pos_assays = frame.get("n_positive_assay_records")
    n_distinct_pairs = len(frame[["peptide", "allele"]].drop_duplicates())
    rule = report.get("label_rule", "(none: built before the rules existed)")
    print(f"--- {name}")
    print(f"  label rule                 {rule}")
    print(f"  rows (= pairs)             {len(frame):,}")
    print(f"  distinct (peptide, allele) {n_distinct_pairs:,}")
    print(f"  distinct peptides          {frame['peptide'].nunique():,}")
    print(f"  distinct alleles           {frame['allele'].nunique():,}")
    print(f"  label 1 / label 0          {int((frame['label'] == 1).sum()):,} / "
          f"{int((frame['label'] == 0).sum()):,}")
    if n_pos_assays is not None:
        n_assays = frame["n_assay_records"]
        n_neg_assays = n_assays - n_pos_assays
        print(f"  assay records behind pairs {int(n_assays.sum()):,}")
        print("  pairs kept, by pattern     "
              f"unanimous_positive {int((n_neg_assays == 0).sum()):,}, "
              f"unanimous_negative {int((n_pos_assays == 0).sum()):,}, "
              f"mixed {int(((n_pos_assays > 0) & (n_neg_assays > 0)).sum()):,}, "
              f"tied {int((n_pos_assays == n_neg_assays).sum()):,}")
        summary = report["pair_pattern_summary"]
        print("  pairs reaching the rule    "
              f"unanimous_positive {summary['unanimous_positive']:,}, "
              f"unanimous_negative {summary['unanimous_negative']:,}, "
              f"mixed {summary['mixed']:,}, tied {summary['tied']:,}, "
              f"total {summary['pairs_total']:,}")
        print(f"  pairs dropped by the rule  {report['pairs_dropped_by_label_rule']:,} "
              f"{report['pairs_dropped_by_reason']}")
    print(f"  recorded content hash      {report['content_hash']}")
    recomputed = content_hash(frame, report["criterion"])
    agree = recomputed == report["content_hash"]
    verdict = "match" if agree else "MISMATCH"
    print(f"  recomputed from parquet    {recomputed}  [{verdict}]")
    if not agree:
        raise SystemExit(f"{name}: the parquet does not hash to its recorded value")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", default=str(ROOT / "artefacts" / "datasets"))
    args = parser.parse_args()
    datasets = Path(args.datasets)

    frames: dict[str, pd.DataFrame] = {}
    for rule in LABEL_RULES:
        name = default_artefact_name(rule)
        frame, report = load(datasets, name)
        describe(name, frame, report)
        frames[rule] = frame

    baseline_path = datasets / f"{BASELINE_NAME}.parquet"
    print()
    if not baseline_path.exists():
        print(f"baseline {baseline_path} is absent; the regression guard cannot run")
        return 1
    baseline, baseline_report = load(datasets, BASELINE_NAME)
    describe(BASELINE_NAME, baseline, baseline_report)

    print()
    print("--- regression guard: strict vs the pre-rule artefact")
    strict_set = pair_label_set(frames["strict"])
    baseline_set = pair_label_set(baseline)
    print(f"  baseline rows              {len(baseline):,}")
    print(f"  baseline distinct triples  {len(baseline_set):,}")
    print(f"  strict rows                {len(frames['strict']):,}")
    print(f"  strict distinct triples    {len(strict_set):,}")
    print(f"  in strict only             {len(strict_set - baseline_set):,}")
    print(f"  in baseline only           {len(baseline_set - strict_set):,}")
    equal = strict_set == baseline_set
    print(f"  identical (peptide, allele, label) set: {equal}")
    if not equal:
        for triple in sorted(strict_set ^ baseline_set)[:20]:
            print(f"    differs: {triple}")
        return 1

    print()
    print("--- pairwise differences between the rules (from the frames)")
    for left in LABEL_RULES:
        for right in LABEL_RULES:
            if left >= right:
                continue
            a, b = pair_label_set(frames[left]), pair_label_set(frames[right])
            label_a = {(p, al): lab for p, al, lab in a}
            label_b = {(p, al): lab for p, al, lab in b}
            keys_a, keys_b = set(label_a), set(label_b)
            shared = keys_a & keys_b
            flipped = sum(1 for key in shared if label_a[key] != label_b[key])
            print(
                f"  {left} vs {right}: pairs only in {left} {len(keys_a - keys_b):,}, "
                f"only in {right} {len(keys_b - keys_a):,}, shared {len(shared):,}, "
                f"shared but oppositely labelled {flipped:,}"
            )

    print()
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
