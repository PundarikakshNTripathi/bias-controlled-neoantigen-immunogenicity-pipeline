# Calibration of the margin test at the real admitted-allele regime

Date: 2026-10-26. Stage 1 of `docs/ROADMAP.md`, decision R3 of `docs/analysis_plan.md`.
This describes the statistical test, not the science: every score here is synthetic noise
or a lookup, laid over the real test set. Nothing was frozen and nothing was written under
`artefacts/splits/`.

## Question and rule

The acceptance gate reads a one-sided p-value from an allele-clustered paired bootstrap
(`nbx.evaluate.stats.margin_over_control`, one-sided alpha 0.05, seeds 1 to 5). When there is
no true margin, how often does it report a significant positive one?

Rule, fixed before any result (R3): the test is acceptable if the upper end of the exact
(Clopper-Pearson) two-sided 95% binomial interval of the false-significance rate is at most
0.07, in both scenarios S1 and S2. Otherwise it is replaced by a percentile-t clustered
bootstrap, then an allele-level permutation test.

## Regime

Source: `docs/results/calibration_S1_B1000.regime.json` (identical in the other regime files),
written by `scripts/calibrate_margin_test.py`.

- Dataset `artefacts/datasets/iedb_tcell_majority.parquet`, sha256 beginning `ad0f63bacc5d27c9`.
- Peptide-grouped split built in memory with the default `SplitCriterion`; split content hash
  beginning `eeed58f4e4fa47a2`.
- 69 admitted alleles; 6,475 test records, 3,137 positive.
- Allele test sizes: minimum 3, median 21, maximum 1,716.

## Scenarios

- S1: model and control scores are independent standard-normal noise per record.
- S2: the control is a per-allele lookup (each allele's training positive fraction, a block of
  ties like a real HLA-only control); the model is the same lookup plus N(0, 0.02^2) noise.
- S3 (supplementary, added after S2 was seen to be conservative, see below): as S2 with noise
  of standard deviation 1e-6, so the allele order never changes. This is the exact null for a
  lookup control. S3 was added after S1 and S2 were declared; the verdict is read from S1 and S2
  and S3 is shown beside it.

A false significance is a positive margin with p below 0.05, the condition
`ControlMargin.significant` applies.

## Results

Command: `python scripts/calibrate_margin_test.py summarise <the five files below>`.
Replicate files (sha256 prefixes): `calibration_S1_B1000.jsonl` `73d4f8ab6d0c0b74`,
`calibration_S2_B1000.jsonl` `00184ea55c8a137e`, `calibration_S3_B1000.jsonl`
`4169bce2e38af04f`, `calibration_S1_B10000.jsonl` `f87fcda4a720f501`,
`calibration_S2_B10000.jsonl` `b4a02d8aaf6c3bb2`, all in `docs/results/`.

Main run, B = 1,000 resamples, 1,500 replicates per scenario:

| Scenario | False significances | Rate | Clopper-Pearson 95% | Upper at most 0.07 |
|:---|---:|---:|:---|:---|
| S1 | 77 | 0.0513 | 0.0407 to 0.0637 | met |
| S2 | 0 | 0.0000 | 0.0000 to 0.0025 | met |
| S3 (supplementary) | 85 | 0.0567 | 0.0455 to 0.0696 | met, narrowly |

Spot check, B = 10,000 resamples (the production setting), 100 replicates per scenario:

| Scenario | False significances | Rate | Clopper-Pearson 95% | Upper at most 0.07 |
|:---|---:|---:|:---|:---|
| S1 | 7 | 0.0700 | 0.0286 to 0.1389 | not met by construction |
| S2 | 0 | 0.0000 | 0.0000 to 0.0362 | met |

## Verdict

**R3 rule met on S1 and S2 (the main run). The test is unchanged; the replacement path was not
needed.**

Limits that belong beside the verdict:

- The spot checks cannot meet a 0.07 upper bound at 100 replicates: even a calibrated test
  (true rate 0.05) gives an upper bound near 0.11 at that size. The R3 rule applies to the
  main run. The S1 spot check at B = 10,000 (7 of 100) is consistent with the main run's rate
  and with 0.05, and it is not evidence of a problem or of its absence. It confirms only that
  the larger B was not obviously worse.
- S2 is conservative, not calibrated. Its noise can swap alleles whose training rates are
  close, which gives a small systematic negative margin (mean delta -0.000474 over 1,500
  replicates, a positive delta in 22% of them). A rate of zero there is a weak check against
  liberalness. That is why S3 was added.
- S3 is the sharper test of the realistic null and it is the closest of the three to the
  boundary: observed rate 0.0567, interval 0.0455 to 0.0696. The interval contains 0.05 and its
  upper end is below 0.07 by 0.0004. S3 meets the rule, but a slightly liberal test (a true
  rate near 0.06) is not excluded by these data. Gate results should be read with this in mind.
- The scenarios have no real peptide signal and no real allele-score structure beyond a lookup.
  They do not cover a model with partial skill; the margin test's power is not measured here.
- Sources for the per-scenario delta and p-value quantiles quoted in the log of this stage are
  the replicate files above; the p-values of S1 and S3 at B = 1,000 are close to uniform on
  (0, 1), as a calibrated test gives.

## Cost

Mean seconds per replicate, from the `seconds` field of the replicate files (several jobs
shared the machine, so these are upper figures): S1 about 57 at B = 1,000 and about 798 at
B = 10,000; S2 about 35 and about 277; S3 about 42 at B = 1,000.
