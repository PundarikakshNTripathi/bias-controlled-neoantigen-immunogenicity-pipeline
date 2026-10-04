#!/usr/bin/env python
"""PostToolUse: run the fast checks after a Python edit and feed back failures.

Keeps the loop tight. A lint error or a broken integrity test is surfaced
immediately as additionalContext, while the edit is still in working memory,
instead of at commit time when the reasoning has to be reconstructed.

Only the fast lane runs here: ruff on the changed file, plus the integrity test
modules when something they cover has changed. The full suite belongs in
`python scripts/verify.py` and in CI.

Never blocks. A failing check is information, not a veto — the point is that
the next turn sees it.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

INTEGRITY_SOURCES = ("nbx/bias/", "nbx/evaluate/", "nbx/splits/")
INTEGRITY_TESTS = [
    "tests/test_bias_diagnostics.py",
    "tests/test_gate.py",
    "tests/test_splits.py",
    "tests/test_metrics_and_policy.py",
]
TIMEOUT = 120


def run(cmd: list[str], cwd: Path) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True, timeout=TIMEOUT
        )
        return proc.returncode, (proc.stdout + proc.stderr).strip()
    except FileNotFoundError:
        return 0, ""           # tool absent; not this hook's problem
    except subprocess.TimeoutExpired:
        return 1, f"timed out after {TIMEOUT}s"


def emit(context: str) -> None:
    json.dump(
        {
            "hookSpecificOutput": {
                "hookEventName": "PostToolUse",
                "additionalContext": context,
            }
        },
        sys.stdout,
    )
    sys.exit(0)


def main() -> None:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        sys.exit(0)

    project = Path(payload.get("cwd") or ".").resolve()
    tool_input = payload.get("tool_input", {}) or {}
    file_path = str(tool_input.get("file_path", ""))

    if not file_path.endswith(".py"):
        sys.exit(0)

    rel = file_path.replace("\\", "/")
    problems: list[str] = []

    code, out = run([sys.executable, "-m", "ruff", "check", file_path], project)
    if code != 0 and out:
        problems.append(f"ruff:\n{out[:1500]}")

    if any(marker in rel for marker in INTEGRITY_SOURCES):
        existing = [t for t in INTEGRITY_TESTS if (project / t).exists()]
        if existing:
            code, out = run(
                [sys.executable, "-m", "pytest", "-q", "--no-header", *existing],
                project,
            )
            if code != 0:
                tail = "\n".join(out.splitlines()[-25:])
                problems.append(
                    "INTEGRITY TESTS FAILING — these pin the project's core "
                    f"claims and must not stay red:\n{tail}"
                )

    if problems:
        emit(
            "Fast checks after your edit to "
            + rel
            + ":\n\n"
            + "\n\n".join(problems)
            + "\n\nFix these before moving on."
        )
    sys.exit(0)


if __name__ == "__main__":
    main()
