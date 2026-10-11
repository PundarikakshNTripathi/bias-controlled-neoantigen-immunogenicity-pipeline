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

(empty)
