#!/usr/bin/env python
"""Calibrate the margin test at the real admitted-allele regime. Freezes nothing.

The acceptance gate reads a one-sided p-value from an allele-clustered paired
bootstrap (`nbx.evaluate.stats.margin_over_control`). This script measures how often
that test reports a significant positive margin when there is NO true margin, at the
allele sizes and label balance of the peptide-grouped split of the primary corpus
(built in memory with the default `SplitCriterion`; nothing is written under
`artefacts/splits/`).

Decision rule (analysis plan, decision R3), fixed before any result: the test is
acceptable if the upper end of the exact (Clopper-Pearson) two-sided 95% binomial
interval of the false-significance rate, at one-sided alpha 0.05, is at most 0.07.
It must hold in both scenarios below. Otherwise the test is replaced, in this order,
by a percentile-t clustered bootstrap and then an allele-level permutation test.

Scenarios (the test set, its labels and alleles are the real ones; only scores vary):

  S1  model and control scores are independent standard-normal noise per record.
  S2  the control is a per-allele lookup (every record of an allele scores the
      allele's training positive fraction: a block of ties, the shape of a real
      HLA-only control); the model is the same lookup plus N(0, 0.02^2) noise, so it
      has no peptide skill.

Each replicate draws fresh scores for seeds 1 to 5 and a fresh set of bootstrap
resamples. A replicate is a false significance when delta > 0 and p < 0.05, the same
condition `ControlMargin.significant` applies.

Resumable: one JSON line per replicate; replicates already in the output file are
skipped, and every replicate's randomness depends only on its number.

    python scripts/calibrate_margin_test.py run --scenario S2 --n-rep 20 --n-boot 10000 \
        --out docs/results/calibration_pilot_S2.jsonl
    python scripts/calibrate_margin_test.py summarise docs/results/calibration_S1.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DATASET = ROOT / "artefacts" / "datasets" / "iedb_tcell_majority.parquet"
SEEDS = (1, 2, 3, 4, 5)
ALPHA = 0.05
RULE_UPPER_BOUND = 0.07
LOOKUP_NOISE_SD = 0.02
SCENARIOS = {"S1": 1, "S2": 2}

_REGIME: dict[str, np.ndarray] = {}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_regime(dataset: Path) -> dict[str, Any]:
    """The real test set of the grouped split: labels, allele index, training rates."""
    import pandas as pd

    from nbx.splits.freeze import SplitCriterion, build_split

    df = pd.read_parquet(dataset)
    alleles = df["allele"].astype(str).tolist()
    labels = df["label"].astype(int).tolist()
    peptides = df["peptide"].astype(str).tolist()
    split = build_split(alleles, labels, SplitCriterion(), groups=peptides)

    test_alleles = [alleles[i] for i in split.test_idx]
    test_labels = np.array([labels[i] for i in split.test_idx], dtype=np.int64)
    names = sorted(set(test_alleles))
    pos = {a: j for j, a in enumerate(names)}
    allele_idx = np.array([pos[a] for a in test_alleles], dtype=np.int64)

    train_pos: Counter[str] = Counter()
    train_all: Counter[str] = Counter()
    for i in split.train_idx:
        train_all[alleles[i]] += 1
        train_pos[alleles[i]] += labels[i]
    overall = sum(train_pos.values()) / max(sum(train_all.values()), 1)
    rate = np.array(
        [train_pos[a] / train_all[a] if train_all[a] else overall for a in names]
    )
    sizes = np.bincount(allele_idx)
    return {
        "labels": test_labels,
        "allele_idx": allele_idx,
        "lookup": rate[allele_idx],
        "meta": {
            "dataset": str(dataset),
            "dataset_sha256": file_sha256(dataset),
            "split_content_hash": split.content_hash,
            "criterion": split.criterion.canonical_json(),
            "n_test_records": int(test_labels.size),
            "n_test_positives": int(test_labels.sum()),
            "n_admitted_alleles": len(split.admitted_alleles),
            "n_test_alleles": len(names),
            "allele_test_size_min": int(sizes.min()),
            "allele_test_size_median": float(np.median(sizes)),
            "allele_test_size_max": int(sizes.max()),
            "alleles_without_training_records": int(
                sum(1 for a in names if not train_all[a])
            ),
        },
    }


def _init_worker(labels: np.ndarray, allele_idx: np.ndarray, lookup: np.ndarray) -> None:
    _REGIME["labels"] = labels
    _REGIME["allele_idx"] = allele_idx
    _REGIME["lookup"] = lookup


def run_replicate(task: tuple[str, int, int]) -> dict[str, Any]:
    from nbx.evaluate import stats

    scenario, rep, n_boot = task
    labels, allele_idx, lookup = (
        _REGIME["labels"], _REGIME["allele_idx"], _REGIME["lookup"]
    )
    n = labels.size
    rng = np.random.default_rng([20261026, SCENARIOS[scenario], rep])
    if scenario == "S1":
        model = [rng.normal(size=n) for _ in SEEDS]
        control = [rng.normal(size=n) for _ in SEEDS]
    else:
        model = [lookup + rng.normal(scale=LOOKUP_NOISE_SD, size=n) for _ in SEEDS]
        control = [lookup.copy() for _ in SEEDS]
    counts = stats.cluster_counts(
        int(allele_idx.max()) + 1, n_boot, seed=900_000 + 1000 * SCENARIOS[scenario] + rep
    )
    t0 = time.perf_counter()
    est = stats.margin_over_control(model, control, labels, allele_idx, counts)
    return {
        "scenario": scenario,
        "rep": rep,
        "n_boot": n_boot,
        "delta": est.delta,
        "p_value": est.p_value,
        "false_significance": bool(
            est.delta > 0.0 and est.p_value is not None and est.p_value < ALPHA
        ),
        "n_valid_resamples": est.n_valid_resamples,
        "seconds": round(time.perf_counter() - t0, 3),
    }


def clopper_pearson(hits: int, n: int, level: float = 0.95) -> tuple[float, float]:
    from scipy.stats import beta

    a = (1.0 - level) / 2.0
    lo = 0.0 if hits == 0 else float(beta.ppf(a, hits, n - hits + 1))
    hi = 1.0 if hits == n else float(beta.ppf(1.0 - a, hits + 1, n - hits))
    return lo, hi


def read_lines(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8").splitlines()
    return [json.loads(x) for x in lines if x.strip()]


def cmd_run(args: argparse.Namespace) -> int:
    out = args.out if args.out.is_absolute() else ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    if "artefacts" in out.parts and "splits" in out.parts:
        raise SystemExit("refusing to write under artefacts/splits/")
    regime = build_regime(args.dataset)
    meta_path = out.with_suffix(".regime.json")
    meta_path.write_text(json.dumps(regime["meta"], indent=2), encoding="utf-8")
    print("regime:", json.dumps(regime["meta"]))

    done = {int(r["rep"]) for r in read_lines(out) if r["scenario"] == args.scenario}
    todo = [
        (args.scenario, rep, args.n_boot)
        for rep in range(args.n_rep)
        if rep not in done
    ]
    print(f"{args.scenario}: {len(done)} replicates already in {out.name}, "
          f"{len(todo)} to run")
    t0 = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_init_worker,
        initargs=(regime["labels"], regime["allele_idx"], regime["lookup"]),
    ) as pool, out.open("a", encoding="utf-8") as handle:
        for k, row in enumerate(pool.map(run_replicate, todo, chunksize=1), start=1):
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            if k % 25 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)} done, {time.perf_counter() - t0:.0f}s elapsed",
                      flush=True)
    return 0


def cmd_summarise(args: argparse.Namespace) -> int:
    rows: list[dict[str, Any]] = []
    for p in args.files:
        rows += read_lines(p if p.is_absolute() else ROOT / p)
    verdicts = {}
    for scenario in sorted({str(r["scenario"]) for r in rows}):
        of_scenario = [r for r in rows if r["scenario"] == scenario]
        for n_boot in sorted({int(r["n_boot"]) for r in of_scenario}):
            sel = [r for r in of_scenario if int(r["n_boot"]) == n_boot]
            reps = {int(r["rep"]) for r in sel}
            n, hits = len(reps), sum(1 for r in sel if r["false_significance"])
            lo, hi = clopper_pearson(hits, n)
            secs = [float(r["seconds"]) for r in sel]
            meets = hi <= RULE_UPPER_BOUND
            verdicts[(scenario, n_boot)] = meets
            print(f"{scenario} B={n_boot}: replicates={n} false significances={hits} "
                  f"rate={hits / n:.4f} CP95=[{lo:.4f}, {hi:.4f}] "
                  f"rule(upper<={RULE_UPPER_BOUND})={'MET' if meets else 'NOT MET'} "
                  f"seconds/replicate mean={np.mean(secs):.2f} max={np.max(secs):.2f}")
    return 0 if all(verdicts.values()) else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--scenario", choices=sorted(SCENARIOS), required=True)
    r.add_argument("--n-rep", type=int, required=True)
    r.add_argument("--n-boot", type=int, required=True)
    r.add_argument("--workers", type=int, default=6)
    r.add_argument("--dataset", type=Path, default=DEFAULT_DATASET)
    r.add_argument("--out", type=Path, required=True)
    r.set_defaults(fn=cmd_run)
    s = sub.add_parser("summarise")
    s.add_argument("files", type=Path, nargs="+")
    s.set_defaults(fn=cmd_summarise)
    args = parser.parse_args()
    return int(args.fn(args))


if __name__ == "__main__":
    sys.exit(main())
