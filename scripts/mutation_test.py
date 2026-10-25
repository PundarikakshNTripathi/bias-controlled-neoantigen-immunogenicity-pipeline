#!/usr/bin/env python
"""Mutation test of the evaluation module's tests. The repository source is never edited.

Each mutation is one deliberate defect in `src/nbx/evaluate/` (the gate, the protocol
checks or the resampling statistics). For each one the script copies `pyproject.toml`,
`src/` and `tests/` to a temporary directory, applies the edit there (the target text
must occur exactly once), and runs the evaluation tests against the copy. A mutant is

  killed     some test failed, so the suite noticed the defect;
  survived   every test passed, so the suite did not notice it;
  error      pytest could not run (collection error, timeout): not counted as killed.

Before any mutant runs, the unmutated copy must pass, and `nbx` must import from the
copy, not from the editable install. A survivor is either a gap in the tests, to be
closed with a test that fails on the mutant, or an equivalent mutant (the defect cannot
change behaviour), which is listed with the reason.

    python scripts/mutation_test.py --out docs/results/mutation.json
    python scripts/mutation_test.py --only M01 M17
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEST_FILES = [
    "tests/test_gate.py",
    "tests/test_protocol.py",
    "tests/test_metrics_and_policy.py",
]
GATE = "src/nbx/evaluate/gate.py"
PROTOCOL = "src/nbx/evaluate/protocol.py"
STATS = "src/nbx/evaluate/stats.py"


@dataclass(frozen=True)
class Mutation:
    mid: str
    what: str
    edits: tuple[tuple[str, str, str], ...] = field(default_factory=tuple)


MUTATIONS = [
    Mutation("M01", "flip the one-sided direction of the bootstrap p-value",
             ((STATS, "float(np.sum(d <= 0.0))", "float(np.sum(d >= 0.0))"),)),
    Mutation("M02", "ignore the HLA-only control: its margin is never required",
             ((GATE, "if not hla_margin.significant:", "if False:"),)),
    Mutation("M03", "pessimistic ties: no other positive of a tie group ranks ahead",
             ((STATS, "spread = (pos - 1) / np.maximum(size - 1, 1)",
               "spread = 0.0 * pos"),)),
    Mutation("M04", "resample records instead of alleles",
             ((PROTOCOL,
               "                labels_arr,\n                allele_idx,\n"
               "                counts,\n            )",
               "                labels_arr,\n"
               "                np.arange(len(labels_arr)),\n"
               "                stats.cluster_counts(len(labels_arr), n_boot, "
               "bootstrap_seed),\n            )"),)),
    Mutation("M05", "skip seed averaging in the bootstrap (first seed only)",
             ((STATS, "boot = acc / len(models)",
               "boot = (weighted_auprc(models[0], y, a, counts)"
               " - weighted_auprc(controls[0], y, a, counts))"),)),
    Mutation("M06", "a missing kappa interval is treated as satisfying condition (ii)",
             ((GATE,
               "    # --- R2: both mandated null controls must have been run -----",
               "    if kappa_change_upper_bound is None:\n"
               "        kappa_change_upper_bound = -1.0\n"
               "    # --- R2: both mandated null controls must have been run -----"),)),
    Mutation("M07", "accept fewer or other seeds than the pre-registered five",
             ((PROTOCOL, "if seeds_list != sorted(EXPECTED_SEEDS):", "if False:"),)),
    Mutation("M08", "skip the matched-config check (epochs, stopping, tuning)",
             ((PROTOCOL, "if getattr(mc, fld) != getattr(cc, fld):", "if False:"),)),
    Mutation("M09", "skip the requested-split-hash check",
             ((PROTOCOL, "if record.content_hash != split_hash:", "if False:"),)),
    Mutation("M10", "skip the records-digest check",
             ((PROTOCOL,
               "if record.recomputed_records_digest() != record.records_sha256:",
               "if False:"),)),
    Mutation("M11", "skip the finite-score check",
             ((PROTOCOL, "if not np.all(np.isfinite(values)):", "if False:"),)),
    Mutation("M12", "accept runs that are missing test records",
             ((PROTOCOL, "if missing and not extra:", "if False:"),
              (PROTOCOL, "if missing or extra:", "if extra:"))),
    Mutation("M13", "accept a seed file that repeats a record index",
             ((PROTOCOL, "if len(set(idx.tolist())) != len(idx):", "if False:"),)),
    Mutation("M14", "skip the HLA-only constancy check",
             ((PROTOCOL, "if varies.any():", "if False:"),)),
    Mutation("M15", "skip the peptide-only identity check",
             ((PROTOCOL, "if np.ptp(values[js]) != 0.0:", "if False:"),)),
    Mutation("M16", "skip the peptide-overlap-between-partitions check",
             ((PROTOCOL, "if peptide in seen and seen[peptide] != part:", "if False:"),)),
    Mutation("M17", "drop the peptide-only condition (iv)",
             ((GATE, "if peptide_only_significantly_worse:", "if False:"),)),
    Mutation("M18", "drop the HLA-only-anchored shortcut condition (iii)",
             ((GATE, "if hla_kappa_change_upper_bound >= 0.0:", "if False:"),)),
    Mutation("M19", "skip the recomputed-split-hash check",
             ((PROTOCOL, "if recomputed != record.content_hash:", "if False:"),)),
]


def make_copy(dest: Path) -> None:
    dest.mkdir(parents=True)
    shutil.copy2(ROOT / "pyproject.toml", dest / "pyproject.toml")
    shutil.copytree(ROOT / "src", dest / "src",
                    ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
    shutil.copytree(ROOT / "tests", dest / "tests",
                    ignore=shutil.ignore_patterns("__pycache__"))


def env_for(copy: Path) -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(copy / "src")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def imports_from_copy(copy: Path) -> bool:
    out = subprocess.run(
        [sys.executable, "-c", "import nbx, sys; print(nbx.__file__)"],
        cwd=copy, env=env_for(copy), capture_output=True, text=True, timeout=120,
    )
    path = Path(out.stdout.strip()).resolve()
    return out.returncode == 0 and copy.resolve() in path.parents


def run_tests(copy: Path, timeout: int) -> tuple[int, str, float]:
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "pytest", *TEST_FILES, "-x", "-q",
             "-p", "no:cacheprovider"],
            cwd=copy, env=env_for(copy), capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return -1, "timeout", time.perf_counter() - t0
    return proc.returncode, proc.stdout + proc.stderr, time.perf_counter() - t0


def first_failure(output: str) -> str:
    for line in output.splitlines():
        if line.startswith("FAILED") or line.startswith("ERROR"):
            return line.strip()[:200]
    return ""


def apply(copy: Path, edits: tuple[tuple[str, str, str], ...]) -> None:
    for rel, old, new in edits:
        path = copy / rel
        text = path.read_text(encoding="utf-8")
        if text.count(old) != 1:
            raise SystemExit(
                f"{rel}: target text occurs {text.count(old)} times: {old!r}"
            )
        path.write_text(text.replace(old, new), encoding="utf-8")


def run_mutant(args: tuple[Mutation, Path, int]) -> dict[str, object]:
    mut, work, timeout = args
    copy = work / mut.mid
    make_copy(copy)
    apply(copy, mut.edits)
    if not imports_from_copy(copy):
        return {"id": mut.mid, "what": mut.what, "status": "error",
                "detail": "nbx did not import from the mutated copy"}
    code, output, secs = run_tests(copy, timeout)
    if code == 0:
        status = "survived"
    elif code == 1:
        status = "killed"
    else:
        status = "error"
    return {"id": mut.mid, "what": mut.what, "status": status,
            "first_failure": first_failure(output), "seconds": round(secs, 1),
            "pytest_exit_code": code}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=ROOT / "docs/results/mutation.json")
    parser.add_argument("--only", nargs="*", default=None)
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1500)
    args = parser.parse_args()

    chosen = [m for m in MUTATIONS if not args.only or m.mid in args.only]
    work = Path(tempfile.mkdtemp(prefix="nbx_mutation_"))
    try:
        base = work / "base"
        make_copy(base)
        if not imports_from_copy(base):
            raise SystemExit("nbx does not import from the unmutated copy; aborting")
        code, output, secs = run_tests(base, args.timeout)
        print(f"unmutated copy: pytest exit {code} in {secs:.0f}s")
        if code != 0:
            print(output[-2000:])
            raise SystemExit("the unmutated copy does not pass; aborting")

        with ThreadPoolExecutor(max_workers=args.parallel) as pool:
            results = list(pool.map(run_mutant,
                                    [(m, work, args.timeout) for m in chosen]))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    for r in results:
        print(f"{r['id']}  {r['status']:<9} {r['what']}  {r.get('first_failure', '')}")
    by = {s: sum(1 for r in results if r["status"] == s)
          for s in ("killed", "survived", "error")}
    print(by)
    control = next((r for r in results if r["id"] == "M01"), None)
    if control is not None and control["status"] != "killed":
        raise SystemExit("control mutant M01 was not killed: the harness is not "
                         "testing the mutated copy; results are void")
    out = args.out if args.out.is_absolute() else ROOT / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(
        {"test_files": TEST_FILES, "summary": by, "results": results}, indent=2),
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
