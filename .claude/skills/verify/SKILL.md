---
name: verify
description: Run the full verification loop and report the real output. Use before declaring any module finished, before any commit that touches the integrity core, and whenever a previous session claimed something works.
allowed-tools: Bash(python scripts/verify.py:*) Bash(python -m pytest:*) Bash(nbx:*) Read Grep
disable-model-invocation: false
---

## Verify

Establish ground truth by execution. Do not summarise; show.

1. Run `python scripts/verify.py`.
2. Run `nbx status`.
3. If anything fails, paste the failure verbatim and stop. Report what broke
   and what it means — do not begin fixing until the failure is stated.
4. If everything passes, paste the summary lines and say explicitly which
   modules are built and which remain contracts only.

Never write "tests pass" without the line of output that says so.

If `$ARGUMENTS` names a module or a file, narrow the run to it and say that you
narrowed it.
