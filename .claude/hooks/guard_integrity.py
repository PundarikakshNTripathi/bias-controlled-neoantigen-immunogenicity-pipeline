#!/usr/bin/env python
"""PreToolUse guard: refuse edits that would violate the project's integrity rules.

Reads the hook payload on stdin, writes a PreToolUse decision on stdout.

Three things are blocked:

1.  Any write under `artefacts/splits/`. A frozen split is immutable; a
    different criterion is a different experiment with a different hash.

2.  Adding a bypass parameter to the acceptance gate (`force=`, `override=`,
    `skip_controls=`, `ignore_kappa=`). The gate exists precisely so that it
    cannot be argued with at 2am.

3.  Introducing what looks like a fabricated result: a metric name followed by
    a concrete number, in a file whose job is to report results. An
    illustrative value is fine when it says so; a measured value is fine when
    it cites the artefact it came from.

Exit code is always 0. The decision travels in the JSON body, which is how
PreToolUse hooks are expected to answer.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import PurePosixPath

# --------------------------------------------------------------------------
# What counts as a result-bearing file. Source and tests are exempt: a test
# asserting `== pytest.approx(0.87)` is pinning behaviour, not claiming a
# finding.
# --------------------------------------------------------------------------
RESULT_BEARING_SUFFIXES = (".md", ".rst", ".txt", ".ipynb")
EXEMPT_PATH_PARTS = ("tests", "test_", "fixtures", "conftest")

METRIC_CLAIM = re.compile(
    r"""(?ix)
    \b(
        auprc | auroc | au-?roc | au-?prc
      | kappa | κ
      | mcc | matthews
      | ppv@?\d* | precision@?\d* | recall@?\d*
      | fraction[\s_-]?ranked
      | accuracy | f1[\s_-]?score
      | p[\s_-]?value
    )\b
    \s*(?:=|:|\bof\b|\bwas\b|\bis\b|\breached\b)\s*
    (?:
        0?\.\d+            # 0.873 or .873
      | \d{1,3}(?:\.\d+)?\s*%   # 87.3 %
      | \d{1,3}\.\d+
    )
    """
)

# An edit is allowed to carry a number when it says where the number came from.
PROVENANCE_MARKERS = (
    "illustrative",
    "example only",
    "not measured",
    "schematic",
    "placeholder",
    "artefact:",
    "artefacts/",
    "doi:",
    "arxiv:",
    "cited",
    "reported in",
    "published",
    "source:",
    "sha256:",
)

BYPASS_PATTERNS = re.compile(
    r"""(?ix)
    \b(
        force
      | override
      | skip[_-]?controls?
      | ignore[_-]?kappa
      | bypass[_-]?gate
      | disable[_-]?controls?
      | allow[_-]?unranked
    )\s*(?::\s*bool)?\s*=
    """
)

GATE_FILES = ("evaluate/gate.py", "evaluate\\gate.py")


def _deny(reason: str) -> None:
    json.dump(
        {
            "hookSpecificOutput": {
                "hookEventName": "PreToolUse",
                "permissionDecision": "deny",
                "permissionDecisionReason": reason,
            }
        },
        sys.stdout,
    )
    sys.exit(0)


def _allow() -> None:
    sys.exit(0)


def main() -> None:
    try:
        payload = json.load(sys.stdin)
    except Exception:
        _allow()  # never break the session over a malformed payload

    tool = payload.get("tool_name", "")
    if tool not in ("Write", "Edit", "NotebookEdit", "MultiEdit"):
        _allow()

    tool_input = payload.get("tool_input", {}) or {}
    raw_path = str(tool_input.get("file_path") or tool_input.get("notebook_path") or "")
    if not raw_path:
        _allow()

    norm = raw_path.replace("\\", "/")
    path = PurePosixPath(norm)
    lower = norm.lower()

    # --- rule 1: frozen splits are immutable -----------------------------
    if "/artefacts/splits/" in lower or lower.endswith("artefacts/splits"):
        _deny(
            "Blocked by project rule 4: files under artefacts/splits/ are frozen. "
            "A split is identified by the content hash of its definition, so a "
            "change to the criterion must produce a NEW split via "
            "`nbx split freeze`, not an edit to an existing one. "
            "If the criterion genuinely needs to change, create a new split and "
            "record why in the commit message."
        )

    # --- what text is this edit introducing? -----------------------------
    added = " ".join(
        str(tool_input.get(k, ""))
        for k in ("content", "new_string", "new_source", "replace_all")
    )
    if not added.strip():
        _allow()

    # --- rule 3: no bypass parameters on the gate ------------------------
    if any(g in lower for g in GATE_FILES) or "evaluate/gate" in lower:
        match = BYPASS_PATTERNS.search(added)
        if match:
            _deny(
                f"Blocked by project rule 3: this edit introduces "
                f"`{match.group(0).strip()}` into the acceptance gate. The gate has "
                "no override by design — that is the single property the whole "
                "project exists to demonstrate. If the gate is wrong, change the "
                "RULE in gate.py and its truth table in tests/test_gate.py "
                "deliberately, and say so in the commit."
            )

    # --- rule 1: no fabricated metric values in result-bearing prose -----
    if path.suffix.lower() in RESULT_BEARING_SUFFIXES and not any(
        part in lower for part in EXEMPT_PATH_PARTS
    ):
        claim = METRIC_CLAIM.search(added)
        if claim:
            context = added.lower()
            if not any(marker in context for marker in PROVENANCE_MARKERS):
                _deny(
                    f"Blocked by project rule 1: this edit writes "
                    f"`{claim.group(0).strip()}` into {path.name}, which reports "
                    "results, without saying where the number came from.\n\n"
                    "No experiment has produced a measured value yet. Either:\n"
                    "  - render it as an em dash (the honest empty state), or\n"
                    "  - mark it illustrative if it is an example, or\n"
                    "  - cite the artefact or publication it came from "
                    "(e.g. 'artefacts/runs/<hash>' or 'reported in Zhang et al. 2026').\n\n"
                    "This guard is the mechanical form of the project's own thesis."
                )

    _allow()


if __name__ == "__main__":
    main()
