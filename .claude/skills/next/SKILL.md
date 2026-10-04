---
name: next
description: Decide what to work on next, based on the real state of the repo and the semester deadlines. Use at the start of a working session.
allowed-tools: Read Grep Glob Bash(git log:*) Bash(git status:*) Bash(nbx:*) Bash(python scripts/verify.py:*)
---

## What next

Work out the real state, then propose one thing.

1. Run `nbx status` and `git log --oneline -10`.
2. Check `python scripts/verify.py --fast` is green. If it is not, that is the
   next thing, and nothing else is.
3. Read the priority order below and find the first item not yet done.
4. Propose **one** concrete next task, small enough to finish in a sitting, and
   say what "done" means for it in terms that can be checked by running a
   command.

## Priority order

The ordering follows the project's own logic: the protocol ships before the
model, and bias measurement precedes training.

1. **M1 data harmonisation** — nothing downstream is interpretable without a
   canonical corpus carrying provenance.
2. **M2 applied to real data** — the bias diagnostics are built and tested, but
   have not been run on an actual corpus. Computing per-allele ρ on a public
   dataset and showing the distribution is the first genuine result, and it
   evaluates H1 and H2.
3. **M6 evaluation harness orchestration** — the gate, the metrics and the
   split construction exist; the runner that executes a variant with its
   matched controls across seeds does not.
4. **M4 model core with null-control variants** — including reproducing the
   published baselines before any immunogenicity claim (FR-44).
5. **M5 adaptive masking** — the policy is implemented; the ablation against
   uniform masking is not.
6. **M3 representation service** — the geometric channel is the expensive part;
   plan compute before starting.
7. **M7 release packaging.**

## Reminder

Deliverable D1 is the benchmark protocol, not a model. If a choice is between
making the protocol solid and making a model look better, the protocol wins.
