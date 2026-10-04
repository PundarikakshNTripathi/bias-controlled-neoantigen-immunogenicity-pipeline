#!/usr/bin/env python
"""The verification loop. One entry point, used by humans, hooks and CI.

Checks, in increasing cost:

  1. integrity    the bias math, the gate, the splits, the metrics
  2. lint         ruff
  3. tests        the full suite
  4. honesty      no fabricated metric values in result-bearing documents
  5. types        mypy (advisory; does not fail the run)

Exit code 0 means every blocking check passed. Anything else means stop.

Run `python scripts/verify.py --fast` for checks 1, 2 and 4 only.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

INTEGRITY_TESTS = [
    "tests/test_bias_diagnostics.py",
    "tests/test_gate.py",
    "tests/test_splits.py",
    "tests/test_metrics_and_policy.py",
]

# Documents that report results, and must not carry unsourced numbers.
RESULT_DOCS = ["README.md", "docs", "reports"]

METRIC_CLAIM = re.compile(
    r"(?i)\b(auprc|auroc|kappa|κ|mcc|matthews|ppv@?\d*|recall@?\d*|"
    r"fraction[\s_-]?ranked|accuracy|f1[\s_-]?score|p[\s_-]?value)\b"
    r"\s*(?:=|:|\bof\b|\bwas\b|\bis\b|\breached\b)\s*"
    r"(0?\.\d+|\d{1,3}(?:\.\d+)?\s*%|\d{1,3}\.\d+)"
)
PROVENANCE = (
    "illustrative", "example only", "not measured", "schematic", "placeholder",
    "artefact:", "artefacts/", "doi:", "arxiv:", "cited", "reported in",
    "published", "source:", "sha256:",
)

GREEN, RED, YELLOW, DIM, RESET = (
    "\033[32m", "\033[31m", "\033[33m", "\033[2m", "\033[0m"
)
if sys.platform == "win32":  # keep output readable in plain consoles
    try:
        import colorama  # type: ignore

        colorama.just_fix_windows_console()
    except Exception:
        GREEN = RED = YELLOW = DIM = RESET = ""


def run(name: str, cmd: list[str], blocking: bool = True) -> bool:
    print(f"{DIM}--- {name}{RESET}")
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    except FileNotFoundError:
        print(f"{YELLOW}SKIP{RESET} {name}: {cmd[0]} not available")
        return True
    output = (proc.stdout + proc.stderr).strip()
    if proc.returncode == 0:
        tail = output.splitlines()[-1] if output else "ok"
        print(f"{GREEN}PASS{RESET} {name}  {DIM}{tail[:100]}{RESET}")
        return True
    print(f"{RED}FAIL{RESET} {name}")
    print("\n".join(output.splitlines()[-40:]))
    return not blocking


def check_honesty() -> bool:
    """No metric value in a result-bearing document without a source."""
    print(f"{DIM}--- honesty{RESET}")
    offenders: list[str] = []

    targets: list[Path] = []
    for entry in RESULT_DOCS:
        p = ROOT / entry
        if p.is_file():
            targets.append(p)
        elif p.is_dir():
            targets.extend(p.rglob("*.md"))

    for path in targets:
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        for lineno, line in enumerate(text.splitlines(), start=1):
            match = METRIC_CLAIM.search(line)
            if not match:
                continue
            window = " ".join(
                text.splitlines()[max(0, lineno - 3): lineno + 2]
            ).lower()
            if any(marker in window for marker in PROVENANCE):
                continue
            offenders.append(
                f"{path.relative_to(ROOT)}:{lineno}: {match.group(0).strip()}"
            )

    if offenders:
        print(f"{RED}FAIL{RESET} honesty: unsourced metric value(s)")
        for o in offenders[:20]:
            print(f"  {o}")
        print(
            "\n  A number in a results document needs a source: cite the "
            "artefact or\n  publication, mark it illustrative, or render it as "
            "an em dash."
        )
        return False

    print(f"{GREEN}PASS{RESET} honesty  {DIM}no unsourced metric claims{RESET}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast", action="store_true", help="skip the full suite")
    args = parser.parse_args()

    py = sys.executable
    results: list[bool] = []

    existing = [t for t in INTEGRITY_TESTS if (ROOT / t).exists()]
    if existing:
        results.append(
            run("integrity core", [py, "-m", "pytest", "-q", "--no-header", *existing])
        )
    else:
        print(f"{RED}FAIL{RESET} integrity core: test files are missing")
        results.append(False)

    results.append(run("lint", [py, "-m", "ruff", "check", "src", "tests", "scripts"]))
    results.append(check_honesty())

    if not args.fast:
        results.append(run("full suite", [py, "-m", "pytest", "-q", "--no-header"]))
        run("types (advisory)", [py, "-m", "mypy", "src/nbx"], blocking=False)

    print()
    if all(results):
        print(f"{GREEN}VERIFICATION PASSED{RESET}")
        return 0
    print(f"{RED}VERIFICATION FAILED{RESET} — do not commit or report results")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
