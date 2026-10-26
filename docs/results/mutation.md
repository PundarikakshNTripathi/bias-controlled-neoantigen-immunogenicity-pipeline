# Mutation test of the evaluation tests

Date: 2026-10-26. Stage 1 of `docs/ROADMAP.md`. This describes the test suite, not the science.

## What was done

`scripts/mutation_test.py` introduces one deliberate defect at a time into a temporary copy of
`src/nbx/evaluate/` (the gate, the protocol checks and the resampling statistics) and runs
`tests/test_gate.py`, `tests/test_protocol.py` and `tests/test_metrics_and_policy.py` against the
copy. The repository source is never edited. Before any mutant runs, the unmutated copy must
pass and `nbx` must import from the copy. A mutant is killed when a test fails, and survives when
every test passes. The first mutant, M01 (flip the direction of the one-sided p-value), is a
control: if it were not killed the results would be void.

Command: `python scripts/mutation_test.py --out docs/results/mutation.json`

The 19 mutants are the 17 named in the roadmap, the removal of the new HLA-only-anchored
condition (iii), and the removal of the recomputed-split-hash check (M19), which the roadmap
list folds into "split-hash". They are listed in `scripts/mutation_test.py`.

## First run (suite at commit 63e4caa)

Source: `docs/results/mutation.json`, sha256 beginning `cfbc0cb2c60aee9d`.

- Killed: 16. Survived: 3. Error: 0.
- Survivors: M04 (resample records instead of alleles in the margin bootstrap), M05 (bootstrap
  only the first seed instead of averaging over the five), M13 (accept a seed file that repeats a
  record index).

None of the three is an equivalent mutant: each changes behaviour the analysis plan specifies
(the cluster is the allele, section 5; the statistic is the seed average; the loader refuses a
malformed layout). They were gaps in the tests.

## Tests added for the survivors (commit 0fa26f7)

| Mutant | New test | Why it fails on the mutant |
|:---|:---|:---|
| M04 | `test_the_margin_bootstrap_resamples_alleles_not_records` | the clusters passed to the bootstrap must number the test alleles, with resample multiplicities summing to that number |
| M05 | `test_the_bootstrap_averages_over_every_seed_not_just_the_first` | seed 1 has a real margin and seeds 2 to 5 a reversed one; the seed-averaged margin is negative and must lie inside its own bootstrap interval |
| M13 | `test_refuses_a_seed_file_that_repeats_a_record_index` | a seed file with a repeated index must be refused with no output |

Each test was run on the unmutated code (passes) and against its mutant (fails) before the
full rerun.

## Rerun on the final suite (commit 0fa26f7)

Source: `docs/results/mutation_after_new_tests.json`, sha256 beginning `4312376ef8389f89`.
Log: `docs/results/mutation_after_new_tests.log`.

- Killed: 19. Survived: 0. Error: 0.

## Limits

- Nineteen hand-chosen defects are a small sample of what could go wrong. Killing all of them
  shows the suite notices these defects; it does not show that no defect would escape.
- Only the three evaluation test files were run. The bias, split and harmonisation tests were
  not part of this measurement.
- Some checks were killed by a single test. That is enough to show the check is covered, not
  that it is covered well.
- The mutants change source text, so they test the code paths as written; a defect in the
  data that a trainer would write (a model that cheats in a way no check anticipates) is not
  covered, and the report says that the matching of epochs, stopping metric and tuning effort
  between a model and its controls is declared metadata that the harness compares but cannot
  verify.
