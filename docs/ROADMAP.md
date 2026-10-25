# Roadmap

Eight stages from the current state to the final submission. Run one with `/next`: it plans
the stage in plan mode, waits for approval, then builds, verifies, documents, commits and
pushes it, and updates this file. The same command resumes an interrupted stage.

## Status

Current module: Stage 1
State: IN PROGRESS

| Stage | Title | Tier | State | Finished | Commit |
|---|---|---|---|---|---|
| 1 | Harden the evaluation | A | IN PROGRESS | | |
| 2 | Data and frozen splits | A | TODO | | |
| 3 | The model and its controls | A | TODO | | |
| 4 | H1 on the published benchmarks | A | TODO | | |
| 5 | Full runs and verdicts | A | TODO | | |
| 6 | H2 and robustness | A | TODO | | |
| 7 | Release and documentation | B | TODO | | |
| 8 | Report drafts and final audit | A | TODO | | |

Tier A is needed for the academic pass condition. Tier B should be done. If time or allowance
runs short, cut in this order and say so in the report: H3 (inside Stage 6), the sensitivity
analyses, Stage 7. Never cut the null controls or the honesty rules.

## In progress

Stage 1: Harden the evaluation. Plan approved 2026-10-26. Tick each step in the commit that finishes it.

- [x] 1. Commit the setup files; this checklist in place, State IN PROGRESS.
- [x] 2. R10: hooks use `.venv\Scripts\python.exe`; run each hook by hand; `verify.py --fast`.
- [x] 3. R1 as an added condition (iii): truth table first, then `gate.py` (upper 95% bound of
  |kappa_model| - |kappa_HLA-only| strictly below zero; missing gives PENDING; no override).
- [x] 4. Wire (iii) into `protocol.py` (report block, `kappa_gate_rule` text) with protocol tests.
- [ ] 5. `scripts/calibrate_margin_test.py` (resumable, deterministic) and a 20-replicate pilot at
  B = 10,000; timing recorded; B and replicate count decided from it.
- [ ] 6. Run the calibration (S1 noise-vs-noise, S2 lookup control vs zero-skill model, empirical allele
  sizes from the grouped majority split, target 1,500 replicates each); write
  `docs/results/calibration.md` with the R3 verdict.
- [ ] 7. Only if R3 fails: percentile-t, then allele permutation, test-first, recalibrated.
- [ ] 8. `scripts/mutation_test.py`, 17 roadmap mutations plus removal of (iii); `docs/results/mutation.md`.
- [ ] 9. A test for every non-equivalent survivor (fails on the mutant, passes on the code).
- [ ] 10. Pre-registration addendum in the deviation log of `docs/analysis_plan.md`.
- [ ] 11. Full `verify.py`; `docs/progress.md`, README, CLAUDE.md; roadmap Stage 1 DONE; push.

Notes: `PROMPTS.md` and `starter.md` are not part of this stage and stay uncommitted.

## Decisions already approved (2026-10-11, plan v2 section 3)

R1 gate condition anchored to the HLA-only control. R2 ratified design choices (peptide-only
condition under FAILED_SHORTCUT, seed-free tie expectation, ranking_k 20 and 100, MCC at the
top-n_pos threshold). R3 calibration rule: acceptable if the upper bound of the 95% binomial
interval of the false-significance rate at one-sided alpha 0.05 is at most 0.07, else switch to
a percentile-t clustered bootstrap, then an allele-level permutation test. R4 assay scope:
functional T cell assays primary, all assays as sensitivity. R5 minimum test records per
admitted allele: 20, or 10 if fewer than 40 alleles remain. R6 group-label validation.
R7 primary corpus ImmuBPI, IEDB the replication, overlap measured, licence resolved before
use beyond download. R8 H2: per-allele |rho| against the shift of the peptide-aware model's
allele-mean score, Spearman, at least 1,000 permutations, replicate in the second corpus.
R9 commits and pushes are approved at the end of every stage once verification passes.
R10 hooks use the virtual-environment interpreter. R11 PyTorch in an optional `ml` extra,
newest CUDA build that passes a kernel smoke test, else CPU. R12 H4, language-model embeddings
and structure features dropped; H3 only if time allows. R13 no hyperparameter tuning, seeds 1
to 5, matched settings. R14 the FR-44 target is fixed by literature check before training.

Decisions that depend on data (the assay-type classification, the licence reading, the freeze
criteria, the FR-44 target and tolerance) are not approved in advance. The stage plan shows
them, computed read-only, so approving the plan approves them.

## Stages

### Stage 1: Harden the evaluation (Tier A)
Goal: close the open checks on the evaluation module before any real result depends on it.
Scope:
- Calibrate the margin test at about 69 alleles with unequal sizes, at least 400 replicates, a
  binomial interval, and the R3 rule; replace the test if the rule fails.
- Mutation test of the suite: break gate, protocol and stats in 18 ways in a temporary copy
  (flip the one-sided direction; ignore the HLA-only control; pessimistic ties; resample
  records not alleles; skip seed averaging; missing kappa passes; fewer than five seeds; skip
  matched-config, split-hash, records-digest, finite-score, missing-record, duplicate-index,
  HLA-only-constancy, peptide-only-identity and peptide-overlap checks; drop the peptide-only
  condition). Every survivor the suite should catch becomes a test.
- Pre-registration addendum (deviation log) recording R1 to R14, the H1 definition of
  "competitive" (not significantly worse under the same clustered bootstrap), and early
  stopping on validation AUPRC for every variant (a deviation from the synopsis wording).
- Hooks call `.venv\Scripts\python.exe` (R10).
- Gate condition from R1, test-first.
Done when: calibration and mutation results are logged, the gate truth table includes R1 and no
override exists, the addendum is committed, `python scripts/verify.py` passes.
Agents: none.

### Stage 2: Data and frozen splits (Tier A)
Goal: two harmonised corpora, the published test sets, and frozen splits.
Scope: group validation and the minimum test block (R5, R6); `nbx split freeze`; IEDB v3 with
the ratified assay scope (R4); ImmuBPI harmonisation with allele-name mapping, pair-level
rules, reconciliation to its published row count (FR-06), overlap with IEDB measured; the two
published test sets harmonised with training-overlap rows removed and tool scores kept; freeze
both training corpora and record the hashes; a test that the frozen splits share no peptide.
Decided in the plan, from read-only output shown there: the assay-type classification table,
the licence reading, the freeze criteria with dry-run counts.
Done when: artefacts are hashed and recorded, frozen hashes are recorded, the peptide-overlap
test passes on the real splits.
Agents: data-steward optional (`model: "sonnet"`).

### Stage 3: The model and its controls (Tier A)
Goal: the pre-registered model and its matched controls, trained once end to end.
Scope: PyTorch install and GPU kernel smoke test (R11); FR-44 target and tolerance fixed by one
literature-verifier agent plus the stage's own check of the locator, committed before training
(R14); M3 minimal (tokeniser, segment positions, content-keyed deterministic cache; no
language-model or structure channels); M4 backbone (4 layers, 4 heads, width 128, dropout 0.1,
Adam 5e-5) with HLA-only and peptide-only controls that withhold input structurally, with
invariance and parameter-matching tests; training loop with resumable checkpoints and early
stopping on validation AUPRC; prediction writer in the harness layout; one pilot run with
measured timing.
Done when: controls pass the invariance tests, one real seed has trained and its predictions
load in the harness, timing is recorded. This is tripwire 1 (18 October).
Agents: literature-verifier (`model: "sonnet"`), ml-engineer optional (`model: "sonnet"`).

### Stage 4: H1 on the published benchmarks (Tier A)
Goal: evaluate H1 without training any third-party tool.
Scope: an explicit, tested single-run path for deterministic third-party scores; the per-allele
lookup control built from a corpus's training split only; run each published tool on each
benchmark; pooled and allele-stratified metrics, margin over the control, gate status; outcome
stated per benchmark against the refutation condition; published figures compared only where a
locator exists.
Done when: `docs/results/h1.md` exists with every number tied to a command and a report hash.
Agents: literature-verifier optional (`model: "sonnet"`).

### Stage 5: Full runs and verdicts (Tier A)
Goal: the main results with the gate applied and independently audited.
Scope: FR-44 baseline reproduction against the committed target; the full matrix (model,
HLA-only, peptide-only; seeds 1 to 5; both corpora) as a resumable background run; `nbx eval`
with the split hash and records digest pinned; one bias-auditor pass (`model: "sonnet"`) and
fixes; `docs/results/model.md`.
Done when: results file committed with sources, audit findings resolved or recorded. Runs launched
by 20 October is tripwire 2.
Agents: bias-auditor (`model: "sonnet"`).

### Stage 6: H2 and robustness (Tier A, with Tier B and C parts)
Goal: the second hypothesis and whether conclusions hold up.
Scope: H2 per R8, tested on synthetic data first, each corpus separately, replication stated;
sensitivity analyses chosen in the plan from the measured cost (ungrouped vs grouped splits, label
rules, assay scope, minimum block size), the subset written down before running; H3 only if the
plan, which states the date and remaining time, says go, and then with the validation-split
threshold and within-allele permutation null.
Done when: `docs/results/h2.md` and `sensitivity.md` exist and each conclusion is marked robust, not
robust, or not tested; H3 is evaluated or recorded as not attempted with the reason.

### Stage 7: Release and documentation (Tier B)
Goal: reproducibility and an accurate repository.
Scope: `nbx.release.package`; a constraints file from the working environment; a Dockerfile (built
only if Docker is installed); a reproduction script that states plainly what it does not retrain;
tables and figures generated from artefacts by `nbx report`, each with a table equivalent and a
report hash; `nbx status` made to read an explicit registry; README, CLAUDE.md and progress log
checked against the code; licence and attribution.
Done when: `nbx report` regenerates the tables from artefacts, and no state claim in the repository
is untrue.

### Stage 8: Report drafts and final audit (Tier A)
Goal: evidence-complete chapter drafts and a final check.
Scope: drafts in `docs/report/`, one per chapter in the order of
`02_requirements_and_uiux/10_final_project_walkthrough.md` section 8, every number tagged with its
source, a limits paragraph per chapter; one literature-verifier pass (`model: "sonnet"`) over claims
about prior work; a final bias-auditor and repro-verifier pass (`model: "sonnet"`, one at a time)
comparing claims with artefacts and hashes; the claims table checked; clean tree and green CI.
Done when: audits recorded, hashes match, the release tag proposal is in the plan for approval.
