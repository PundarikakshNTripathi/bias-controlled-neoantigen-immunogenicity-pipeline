#!/usr/bin/env python
"""SessionStart: put the current real state of the project into context.

A session that starts knowing which modules are built, what the git state is,
and whether the integrity core still passes will not spend its first three
turns rediscovering those things — or, worse, assuming them.

Everything here is measured at session start. Nothing is asserted from a
static table.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

MODULES = [
    ("M1", "Data harmonisation", "src/nbx/data/harmonise.py"),
    ("M2", "Bias diagnostics", "src/nbx/bias/diagnostics.py"),
    ("M3", "Representation service", "src/nbx/represent/cache.py"),
    ("M4", "Model core", "src/nbx/model/backbone.py"),
    ("M5", "Masking policy", "src/nbx/bias/policy.py"),
    ("M6", "Evaluation harness", "src/nbx/evaluate/protocol.py"),
    ("M7", "Release packaging", "src/nbx/release/package.py"),
]


def sh(cmd: list[str], cwd: Path) -> str:
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=20)
        return (p.stdout or "").strip()
    except Exception:
        return ""


def main() -> None:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        payload = {}

    project = Path(payload.get("cwd") or ".").resolve()
    lines: list[str] = []

    branch = sh(["git", "rev-parse", "--abbrev-ref", "HEAD"], project)
    dirty = sh(["git", "status", "--porcelain"], project)
    if branch:
        state = f"{len(dirty.splitlines())} uncommitted file(s)" if dirty else "clean"
        lines.append(f"Git: branch `{branch}`, {state}.")

    built, pending = [], []
    for tag, name, rel in MODULES:
        path = project / rel
        if not path.exists():
            pending.append(f"{tag} {name}")
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pending.append(f"{tag} {name}")
            continue
        (built if "NotImplementedError" not in text else pending).append(
            f"{tag} {name}"
        )

    lines.append("Built: " + (", ".join(built) if built else "nothing yet"))
    lines.append("Contract only: " + (", ".join(pending) if pending else "none"))

    # Does the integrity core still hold? Measure, do not assume.
    integrity = [
        t
        for t in ("tests/test_bias_diagnostics.py", "tests/test_gate.py")
        if (project / t).exists()
    ]
    if integrity:
        try:
            proc = subprocess.run(
                [sys.executable, "-m", "pytest", "-q", "--no-header", *integrity],
                cwd=project,
                capture_output=True,
                text=True,
                timeout=120,
            )
            lines.append(
                "Integrity core: PASSING"
                if proc.returncode == 0
                else "Integrity core: **FAILING** — fix before any other work."
            )
        except Exception:
            lines.append("Integrity core: could not be checked this session.")

    lines.append(
        "Reminder: no invented numbers, unmeasured renders as an em dash, "
        "the acceptance gate has no override."
    )

    json.dump(
        {
            "hookSpecificOutput": {
                "hookEventName": "SessionStart",
                "additionalContext": "Project state at session start:\n"
                + "\n".join(f"- {line}" for line in lines),
            }
        },
        sys.stdout,
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
