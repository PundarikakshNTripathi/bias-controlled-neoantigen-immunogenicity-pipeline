---
name: next
description: Run the next roadmap stage end to end. Plans it in plan mode, waits for approval, then builds, verifies, documents, commits and pushes it, and updates docs/ROADMAP.md. Also resumes a stage interrupted by a usage limit. Use at the start of any working session.
allowed-tools: Read Grep Glob Edit Write Bash(git status:*) Bash(git diff:*) Bash(git log:*) Bash(git add:*) Bash(git commit:*) Bash(git push origin main) Bash(python:*) Bash(.venv/Scripts/python.exe:*) Bash(nbx:*)
---

# /next: run one roadmap stage

`docs/ROADMAP.md` lists eight stages and holds the state. This skill moves exactly one stage forward
per run. The researcher types only `/next`.

## 0. Orient (cheap)

1. Read `docs/ROADMAP.md`: the Status block, the "In progress" block, and the section of the
   current stage only. Do not read other stages.
2. Confirm the virtual environment is active (`python -c "import sys; print(sys.prefix)"` shows a
   `.venv` path). If not, say so in one line and stop.
3. `git status --short`. If `docs/ROADMAP.md`, `.claude/skills/next/SKILL.md`, the session hook or
   `CLAUDE.md` are uncommitted setup changes, they are committed first as `chore: add roadmap
   workflow` as soon as plan approval allows it (not before, if still in plan mode).

## 1. Resume, if a stage is IN PROGRESS

If `State: IN PROGRESS`: say in one line which stage and which step is next, run
`git diff --stat`, decide whether any uncommitted work is complete (finish it, or discard it and say
so), and go to step 3. Do not re-plan: the approved checklist in "In progress" is the plan.

## 2. Otherwise plan the current stage, in plan mode

1. Call `EnterPlanMode` (load its schema with `ToolSearch` if it is not available).
2. In plan mode, do the read-only investigation the stage needs. Where the stage lists "decided in
   the plan", compute the evidence read-only and put the concrete proposal in the plan, so approving
   the plan approves the decision. Use `AskUserQuestion` only for a real choice the evidence cannot
   settle.
3. Write the plan: context, a numbered checklist of at most about twelve steps (each with the
   command that proves it and a commit boundary), decisions with recommendations, risks, what will
   be documented, and the final verification. Honour every rule in `CLAUDE.md` and the approved
   decisions in the roadmap.
4. Call `ExitPlanMode`. After approval, copy the checklist into "In progress" in `docs/ROADMAP.md`,
   set `State: IN PROGRESS`, commit the setup files and this change, and go to step 3.

## 3. Execute without further prompting

- Work the checklist in order. Test first for anything load-bearing. After each step run its check,
  commit (conventional prefix, one concern, no attribution lines) and tick the checklist in the same
  commit. Use `python scripts/verify.py --fast` between steps and the full run only at the end.
- Stay inside the approved plan. Do not ask for permission for things the plan already covers. If a
  step cannot be done as written, or needs something outside the plan, stop, say exactly what is
  needed and why, record it under "In progress", and do not improvise a rule or a threshold.
- Agents only where the stage names them: one at a time, always with `model: "sonnet"`, never in the
  background. Check an agent's claims by running the command yourself before acting on them.
- Long computations (calibration, training) run as background commands. Check them at long intervals
  with a cheap listing, not by reading logs repeatedly.
- Never run anything that edits `artefacts/splits/`, adds `force`/`override`/`skip_controls`, or
  writes a number into documentation without a command and artefact that produced it.

## 4. Finish

1. `python scripts/verify.py` (full). Paste the actual output tail. If it fails, fix it; do not mark done.
2. Document: a dated entry in `docs/progress.md` (what was built, measured with commands, caveats,
   what is not done); README State and Commands tables and the `CLAUDE.md` architecture table if
   state changed; a deviation-log entry in `docs/analysis_plan.md` for any departure from the
   pre-registration. Results go in `docs/results/` with the command and report hash beside every table.
3. In `docs/ROADMAP.md`: set the stage to DONE with date and commit hash, empty "In progress", set
   `Current module` to the next TODO stage and `State: TODO`.
4. Commit the documentation, then `git push origin main`. If the push fails, report the error and stop.
5. Final message, at most ten lines: what was done, the key outcomes and limits, the verification
   result, and the name of the next stage. Then stop. Do not begin the next stage.

## Limits

If the session is cut by a usage limit, nothing is lost beyond the step in progress: every finished
step is committed and ticked. The next `/next` resumes from section 1.
