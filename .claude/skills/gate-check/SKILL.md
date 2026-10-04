---
name: gate-check
description: Apply the acceptance rule to a claimed result before it is reported anywhere. Use whenever a model result, comparison or leaderboard entry is about to be written into a document, figure or commit.
allowed-tools: Read Grep Glob Bash(python -m pytest:*) Bash(nbx:*)
---

## Gate check

Before this result goes anywhere, establish that it is allowed to.

Walk the acceptance rule in order and answer each with evidence from an
artefact, not from memory:

1. **Were both null controls run?** HLA-only and peptide-only, at matched
   capacity and training budget. If either is missing, the status is PENDING
   and **no headline metric may be displayed at all** (FR-63). Stop here.
2. **Is the margin over the HLA-only control positive and significant?** If
   not, metrics may be shown but **no rank is assigned**.
3. **Did |κ| rise against the baseline?** If it did, the result is reported as
   *not surviving the controls*. No rank.
4. **Only if 1–3 all pass** may a rank be assigned, and the margin must be
   displayed alongside it. The margin can never be hidden.

Then state the resulting `GateStatus` and quote the exact numbers you used,
with the artefact path each came from.

If the honest answer is PENDING — which it is until an experiment has actually
run — say PENDING. That is the correct state, not a gap to paper over.

Consider delegating to the `bias-auditor` subagent for anything that will
appear in the report or the thesis.
