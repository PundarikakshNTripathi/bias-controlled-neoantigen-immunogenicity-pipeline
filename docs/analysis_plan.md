# Analysis plan (pre-registration)

Written and committed before any model has been trained or evaluated in this
repository. Its purpose is to fix, in advance and in version control, every analytic
choice that could otherwise be adjusted after seeing a result. Any later departure is
recorded in the deviation log at the end, with a reason and a date; the original text
is not edited.

Status of this document: the researcher approved the analysis plan of 2026-10-11 as a
whole. The one choice still marked **pending** (the shortcut-dependence tolerance in
section 5) changes the acceptance gate and needs separate explicit approval; it is not
binding until that marker is removed in a later commit.

## 1. Hypotheses and refutation conditions

Taken from the project synopsis (Table 4.1). Each is paired with the outcome that
refutes it, and a refutation is a reportable result.

| ID | Hypothesis | Refuted if |
|---|---|---|
| H1 | Current immunogenicity benchmarks are not bias-controlled, and an HLA-only null model is competitive with published peptide-HLA models on them. | The HLA-only control performs significantly worse than published models on every evaluated benchmark. |
| H2 | Intra-allelic label imbalance, measured by mutual information between allele and label, predicts the magnitude of a model's shortcut dependence. | No significant association between the mutual-information diagnostic and the measured allele sensitivity of predictions. |
| H3 | Adaptive per-allele masking improves mutation-level discrimination relative to uniform position-2 and C-terminal masking. | Adaptive masking performs no better than, or worse than, uniform masking on held-out mutation-variant evaluation. |
| H4 | Geometric interface features reduce a model's dependence on allele identity relative to a sequence-only model of matched capacity. | Adding geometric features leaves allele sensitivity and the margin over the HLA-only control statistically unchanged. |

H1 and H2 are the academic-deliverable hypotheses and require at least two public
corpora (SRS section 9.2). H3 and H4 are targets, not pass conditions.

## 2. Unit of analysis and labels

- The unit is the (peptide, allele) pair, not the assay row. Rows that share a pair are
  repeated measurements of one observation.
- **Primary label rule: majority vote over the pair's assays, ties
  dropped.** All `Positive*` qualitative levels count as positive. The pair's assay
  count and positive fraction are retained as columns.
- **Sensitivity analyses, declared now:** (a) drop every pair whose assays disagree;
  (c) label a pair positive if any assay is positive. Each is a separate artefact with
  its own content hash. Conclusions are stated for the primary rule and reported as
  robust or not robust against (a) and (c); the primary is not changed after results
  are seen.
- Scope: human HLA class I, linear unmodified peptides, 8 to 14 standard residues,
  alleles resolvable to a two-field name with an available pseudo-sequence. Every
  excluded record is counted under a named reason in the corpus report.

## 3. Splits

- Frozen before any model is trained; a different criterion is a different split with a
  different hash. Files under `artefacts/splits/` are never edited.
- **All records of a peptide go to a single partition.** The test set
  keeps the existing allele-balance criterion (|rho| bounded in the test block; default
  bound and minimum allele size are those in `SplitCriterion`, recorded in the split
  hash). Alleles that cannot meet the criterion are reported separately and never
  pooled into headline metrics.
- A stricter variant that also separates peptides by sequence similarity may be
  reported as a robustness check; it is not the primary split.

## 4. Metrics and the acceptance rule

- Headline discrimination metric: AUPRC. Reported alongside: AUROC, PPV@k, Recall@k,
  the two TESLA ranking metrics, MCC, and per-allele figures aggregated with explicit
  weights (the weighting is stated in every table).
- Shortcut dependence kappa(f) = Spearman(rho(a), E[f(x) | a]) over alleles in the test
  set, computed by `nbx.bias.diagnostics.shortcut_dependence`.
- A variant is rankable only if it exceeds the matched HLA-only control by a
  significant margin and shortcut dependence does not rise, as implemented in
  `nbx.evaluate.gate.evaluate_gate`. There is no override.

## 5. Significance of the margin over the HLA-only control

- Statistic: delta = AUPRC(model) minus AUPRC(HLA-only control) on the same test set and
  the same seed index.
- Test: paired bootstrap that resamples alleles (clusters) with replacement, B = 10,000,
  one-sided, alpha = 0.05, applied to the seed-averaged delta. Per-seed deltas are also
  reported as mean and standard deviation across seeds.
- Only the HLA-only comparison gates. The peptide-only comparison is diagnostic.
- Seeds: 1, 2, 3, 4, 5, fixed now. Results are reported as means with across-seed
  variance, never as the best seed.
- Shortcut-dependence tolerance (pending explicit approval, because it changes the
  gate rule): a rise in |kappa| is declared only if the allele-clustered bootstrap
  interval for the change excludes zero. Until approved, the gate's existing rule
  (point estimate, tolerance zero) applies.

## 6. Matched controls

Every reported variant is accompanied by an HLA-only and a peptide-only control with
the same backbone, capacity, number of epochs, early-stopping criterion, tuning budget
and seeds. Controls withhold their input structurally. They cannot be disabled, and
the evaluation harness refuses to print a headline metric when one is missing.

## 7. Masking policy estimation

- Per-position mutual information within each allele subset is compared with a null
  from 1,000 within-subset label permutations; a position is eligible for masking only
  if its observed value exceeds the 95th percentile of its own null.
- The threshold tau is selected on the validation split, never on the test split.
- Positions are restricted to candidate anchor positions when the comparison with the
  uniform scheme is run, so that T-cell-contact positions are not masked merely because
  they are informative in a biased corpus.

## 8. Baseline reproduction

Before any immunogenicity claim, binding and presentation configurations are reproduced
within a stated tolerance of published figures (FR-44). Each published figure is
recorded with its source and locator beside the reproduced figure. A failed
reproduction is documented and later comparisons use the reproduced baseline, labelled
as such.

## 9. What counts as a negative result

A variant that does not exceed the HLA-only control significantly, or whose kappa rises,
is reported in the results exactly as such. Hyperparameters, thresholds and labelling
rules are not adjusted after seeing test-set outcomes. If adaptive masking does not
beat uniform masking, that is the finding.

## 10. Known limits

Conclusions are bounded by the quality, coverage and heterogeneity of public assay
data, and no new ground truth is generated. A single-source corpus cannot test the
claim that the same allele carries contrasting ratios across datasets; that needs a
second corpus. Plug-in mutual information is biased upward for small groups, which
section 7 addresses. The set of alleles meeting the balance criterion may be small;
if so, that is reported as a finding about data availability.

## Deviation log

### 2026-10-11: harness and gate amendments after the design audit

No model had been trained or evaluated in this repository at the time; every check
described here used synthetic data. An independent adversarial audit of the M6 harness
found that it could rank a model that scores by the shortcut. The amendments below
were approved by the researcher in the master plan of the same date.

1. **Section 5, shortcut-dependence tolerance: withdrawn.** The clause "a rise in
   |kappa| is declared only if the allele-clustered bootstrap interval for the change
   excludes zero" is withdrawn. It was an error in the original drafting: it makes the
   gate more permissive, because a model needs less evidence to avoid a failure. It is
   replaced by a burden-of-proof rule. A variant is rankable only if (i) its margin over
   the matched HLA-only control is significant, (ii) the upper 95% bound of the change
   in |kappa| relative to the baseline, from the allele-clustered bootstrap, does not
   exceed a non-inferiority margin of zero, and (iv) it is not significantly worse than
   the matched peptide-only control. A failure of (ii) or (iv) is reported as
   `FAILED_SHORTCUT` with the failed condition named; filing (iv) under that status is a
   design choice recorded here for ratification.
2. **Condition (iii), an absolute ceiling on |kappa|, is not adopted.** Its value is a
   scientific judgement that has not been given. The report lists it as not done. Until
   it is set, a model exactly as shortcut-dependent as its baseline, with a zero-width
   interval, can satisfy (ii); the report states this limit.
3. **Section 4, ties: exact expectation, seed-free.** Reported AUPRC and the other
   ranking metrics use the exact expectation over uniformly random orderings within tied
   scores, so the tie handling has no seed. The planning note of the same date that
   called this "equivalent" to threshold-grouped average precision was wrong: they differ
   on tied blocks, and the grouped form was measured to be biased against an
   allele-constant control. `metrics.auprc` and its tests are unchanged.
4. **Section 5, reporting pair.** The one-sided p-value (alpha 0.05) is reported with the
   matching 90% two-sided percentile interval and the one-sided 95% lower bound. A 95%
   two-sided interval is not reported beside a one-sided p-value.
5. **Pinned settings.** Trained-model seeds are exactly 1, 2, 3, 4, 5. The production
   command always uses 10,000 bootstrap resamples and a fixed bootstrap seed; the Python
   API may lower the resample count for tests, and the report then records that the
   pre-registered settings were not used.
6. **Multiplicity.** The report states the number of model variants compared and a
   Holm-adjusted p-value for each model's margin. The gate reads the unadjusted
   one-sided p-value; using the adjusted one would be a rule change requiring approval.
7. **Integrity checks added to the harness.** An HLA-only control must give identical
   predictions within every allele, a peptide-only control must give identical scores to
   records that share a peptide, `split.json` carries per-record peptide, allele and
   label with a digest that the harness verifies, and a peptide may not occur in more
   than one of train, validation and test. Matched epochs, parameters, early-stopping
   metric and tuning trials remain declared metadata that the harness cannot verify.

Not done at this date: calibration of the bootstrap at the real admitted-allele count
and unequal allele sizes, and a mutation test of the suite. A first attempt at both,
run by an audit agent, was terminated by a rate limit before it reported, and nothing
from it is used. They remain to be run, and their results will be logged here. Until
then the bootstrap is known to be calibrated only on small synthetic regimes, where an
independent check measured it as mildly liberal, and no statement is made about its
behaviour at the real admitted-allele count.


### 2026-10-26: decisions R1 to R14 recorded; gate condition (iii)

No model has been trained or evaluated in this repository at this date. This entry records,
before any model run, the decisions the researcher approved on 2026-10-11 (master plan v2,
section 3) and the one change to the acceptance gate that follows from them. The original
text above is not edited. Source of the decisions: `docs/ROADMAP.md`.

1. **R1, gate condition (iii), implemented.** The earlier entry (2026-10-11, items 1 and 2)
   left condition (iii) unimplemented and noted that a model exactly as shortcut-dependent as
   its baseline, with a zero-width interval, could satisfy (ii). Condition (iii) is now
   implemented as an added condition: the upper 95% bound of
   `|kappa_model| - |kappa_HLA-only control|`, from the same allele-clustered bootstrap, must
   be strictly below zero. Conditions (ii) and (iv) are unchanged, so the rule only tightens.
   A missing or undefined interval gives `PENDING`. The absolute ceiling on `|kappa|` is
   considered and not adopted, because its value would be an arbitrary constant. Code and
   truth table: commit 74840e0 (`gate.py`, `tests/test_gate.py`, `tests/test_protocol.py`).
2. **R2, ratified design choices.** The peptide-only condition is filed under
   `FAILED_SHORTCUT`; ties use the seed-free expected average precision; ranking_k is 20 and
   100, both reported; MCC uses the top-n_pos threshold.
3. **R3, calibration rule.** The margin test is acceptable if the upper end of the 95%
   binomial interval of its false-significance rate, at one-sided alpha 0.05, is at most
   0.07, at about 69 admitted alleles of unequal size. Otherwise it is replaced by a
   percentile-t clustered bootstrap, then by an allele-level permutation test, whichever
   first meets the rule. The number of replicates is "at least 400". At 400 replicates a test
   whose true rate is 0.05 has an exact upper bound of about 0.076 and would fail the rule
   more often than not, so the run uses more replicates (1,500 per scenario); the 0.07 bound
   and the decision order are unchanged. The interval used is the exact Clopper-Pearson
   interval. Source: `scripts/calibrate_margin_test.py`, written before any result.
4. **R4, assay scope.** Functional T cell assays are primary; all assays are a declared
   sensitivity analysis. The classification of assay types is shown and ratified before any
   model run.
5. **R5, minimum test records per admitted allele.** 20, or 10 if fewer than 40 alleles
   remain at 20. Alleles below it are reported separately and never pooled into headline
   metrics.
6. **R6, group-label validation.** The splitter raises on null, NaN or mixed-type group
   labels.
7. **R7, corpora.** The ImmuBPI training corpus is primary and the IEDB T cell corpus is the
   replication. The overlap between them is measured and reported, and cross-corpus
   statements use de-duplicated pairs. The licence question (MIT in the repository,
   Apache-2.0 on Zenodo) is resolved before any use beyond the local download.
8. **R8, H2 operationalisation.** Per allele, `|rho(a)|` is tested against the size of the
   shortcut in the peptide-aware model's allele-mean score, by Spearman correlation with the
   allele as the unit, at least 1,000 permutations and allele-bootstrap intervals, in each
   corpus separately. H2 is supported only if sign and significance replicate in the second
   corpus. The HLA-only control is reported as the reference upper end, not as the test.
9. **R9, commits and pushes** are approved at the end of every stage once verification
   passes.
10. **R10, hooks** run with the virtual-environment interpreter (commit 9a994ea).
11. **R11, PyTorch** is installed into an optional `ml` extra, newest CUDA build that passes a
    kernel smoke test on the GPU, otherwise CPU.
12. **R12, scope cuts.** H4, language-model embeddings and structure features are dropped from
    the academic deliverable and named as future work. H3 is attempted only if time allows.
13. **R13, no tuning.** No hyperparameter search for any variant; seeds 1 to 5; matched
    settings. The model may underperform a tuned one; no claim of state of the art is made.
14. **R14, FR-44 target.** The published target figure, with source and locator, is fixed by
    a literature check before any training and committed first. If the reproduction misses
    its stated tolerance, later comparisons use the reproduced baseline, labelled as such.

Two further clarifications, both fixed before any evaluation:

- **H1, the meaning of "competitive".** The HLA-only control is competitive with a published
  model on a benchmark if its AUPRC is not significantly lower than that model's under the
  same paired, allele-clustered bootstrap (one-sided, alpha 0.05, opposite tail). The
  refutation condition in section 1 is unchanged. How the per-model outcomes on a benchmark
  are aggregated is written into the Stage 4 plan before any figure is computed.
- **Early stopping.** Every variant, model and controls alike, stops on validation AUPRC.
  The Synopsis wording (a validation AUROC) is older; `configs/ablation.yaml` already names
  `val_auprc`, in line with AUPRC as the headline metric.

Results of the two checks of the evaluation module, run after the entry above (they describe
the test and its test suite, not any model):

- **R3 calibration: rule met, test unchanged.** On the S1 and S2 scenarios at 1,500 replicates
  each, the upper end of the exact 95% interval of the false-significance rate is below 0.07
  (S1 rate 0.0513, upper 0.0637; S2 zero in 1,500, upper 0.0025). S2 is conservative, so a
  supplementary scenario S3 (added after S1 and S2 were declared) gave rate 0.0567 with upper
  end 0.0696: inside the rule by 0.0004, and a slightly liberal test is not excluded. The
  replacement tests (percentile-t, permutation) were not needed and were not built. Source:
  `docs/results/calibration.md`, commit 6bc7ae7.
- **Mutation test: 19 of 19 killed after three added tests.** The first run killed 16 and
  three survived (M04, M05, M13), none equivalent; tests were added in commit 0fa26f7 and the
  rerun killed all 19. This shows the suite notices those 19 defects, not that none could
  escape. Source: `docs/results/mutation.md`, commit 8d5838e.
