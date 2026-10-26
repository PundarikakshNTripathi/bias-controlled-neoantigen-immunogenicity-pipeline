# Progress log

Append one entry per working session. Be specific about what is broken.

---

## 2026-10-04 — repository bootstrap

**Built.** Integrity core: bias diagnostics (I(b;y), rho(a), kappa(f)), adaptive
masking policy M(a), the acceptance gate, allele-balanced split construction and
freezing, and the evaluation metrics (AUPRC, AUROC, PPV@k, Recall@k, the two
TESLA ranking metrics, MCC, allele-stratified aggregation).

**Verified.** 85 tests passing. The mathematics is pinned against analytically
known cases and cross-checked against `scipy.stats.spearmanr`,
`sklearn.metrics.mutual_info_score`, `average_precision_score`, `roc_auc_score`
and `matthews_corrcoef`. Lint and types clean.

**Not built.** M1 harmonisation, M3 representation service, M4 model core, M6
protocol orchestration, M7 release packaging. All have written contracts.

**No experiment has been run.** No corpus has been ingested. Every metric
everywhere renders as an em dash, correctly.

**Next.** M1 on the IEDB corpus, then M2 applied to it — that pair produces the
first real result and evaluates H1 and H2.

---

## 2026-10-05 — M1 Step 0: raw IEDB and pseudo-sequence inspection

No code written. Environment fixed first: the "integrity core FAILING" session
report was a missing `.venv` (system Python lacked pytest/ruff); with the venv,
`python scripts/verify.py` passes.

**Downloaded (gitignored, hashes in `data/raw/MANIFEST.txt`).**
IEDB `tcell_full_v3.zip` and `mhc_ligand_full_single_file.zip`; ImmuBPI
`dataset/const/IEDB_MHC_pseudo.dat`.

**Measured from the real T cell file** (`tcell_full_v3.csv`, 161 columns, two
header rows: group then field; columns must be addressed as `Group|Field`):
- 578,203 assay rows. Epitope object type: 571,323 linear peptide, 6,808
  non-peptidic, 72 discontinuous.
- Qualitative measurement: 363,553 Negative; 185,577 Positive; 18,111
  Positive-Low; 5,827 Positive-High; 5,135 Positive-Intermediate.
  How to map the four Positive* levels to a binary label is a decision to record.
- MHC restriction class: 186,512 class I, 250,633 class II, 3,885 non classical,
  137,173 missing. 132,906 rows have no MHC restriction name.
- 173,133 rows carry an `HLA-` name, 458 distinct strings, but these include class
  II and noise: `HLA-A*02:01 T163A mutant`, `HLA-DR7`, `HLA-Cw12`, `HLA class I`.
  Allele reconciliation (FR-02) is a real task, not a formality.
- `Reference|PMID` missing on 73,347 rows.

**MHC ligand file:** 112 columns, same two-row header, uncompressed size
10,042,937,020 bytes. Must be streamed in chunks, not loaded whole. Not yet
profiled.

**Pseudo-sequences.** `IEDB_MHC_pseudo.dat` has 11,373 lines, every one name +
34-character sequence (includes non-human, e.g. BoLA). Names use the form
`HLA-A02:01` (no asterisk). Its origin is not stated in the ImmuBPI README
(searched for pseudo/netmhc/licen), and the format matches NetMHCpan's
`MHC_pseudo.dat`, so redistribution terms are unverified: it stays in
`data/raw/`, not vendored. Repo licence MIT (GitHub API); Zenodo record
lists Apache-2.0 — discrepancy confirmed, still to resolve.

**Not done / open.** FR-06 reference statistics not yet located. No harmonised
dataset exists. M1 remains contract only.

---

## 2026-10-05 — M1 built for the IEDB T cell export

**Built.** `nbx.data.harmonise` (IEDB `tcell_full_v3` only), `nbx.data.alleles`
(HLA nomenclature reconciliation with named drop reasons), a priority-chain
`PseudoSequenceRegistry`, a content-hashed parquet writer that refuses to
overwrite a differing artefact, and `nbx data harmonise`. `pyarrow` added as a
dependency. 23 tests in `tests/test_harmonise.py` cover FR-01 to FR-06, including
that a record without provenance is dropped rather than defaulted.

**Measured** (`artefacts/datasets/iedb_tcell.report.json`, content hash
`3cc846db…7821`): 578,203 raw rows; 46,608 canonical records; 146 alleles;
28,714 distinct peptides; 10,925 positive and 35,683 negative. Kept, dropped and
merged records sum to the raw row count. `repro-verifier` re-ran the command and
reproduced the hash and a byte-identical artefact, and independently reproduced the
label counts to within four records (see B*44:01 below).

**Findings that affect interpretation.**
- Dropping contradictory-label (peptide, allele) pairs removes 25,034 positive and
  8,611 negative records and moves the positive:negative ratio from about 1:1 to
  about 1:3. The policy is part of the hashed criterion. It is a choice that has not
  been compared against alternatives (keep-and-flag, majority vote).
- `HLA-B*44:01` is absent from the pseudo-sequence table and its four records are
  dropped; the drop reason is mislabelled `allele_two_digit_mapping_ambiguous`.
- The two-digit-representative mapping fired zero times on real data; it is
  exercised only by synthetic tests.
- Serological names (13,240 records) and unspecified-resolution alleles (26,479)
  are dropped.
- 22,220 records have no PMID; they are retained and marked `pmid_unreported`.

**Not done.** `mhc_ligand_full` (about 10 GB uncompressed) is not ingested;
TESLA and HiTIDE are not ingested. FR-06 has no published reference to compare
against. `write_corpus` does not re-hash an existing parquet before refreshing its
report. M2 has not been run on this corpus; no hypothesis has been evaluated.
`nbx bias profile` and `nbx split freeze` are not wired to the CLI.

**Next.** Decide the conflict policy, then run M2 on `iedb_tcell.parquet`.

---

## 2026-10-11 — M2 run on the harmonised IEDB T cell corpus

**Built.** `nbx.bias.profile` (assembles the report from `nbx.bias.diagnostics`,
`nbx.bias.policy` and `nbx.splits.freeze`; no mathematics reimplemented) and the
`nbx bias profile` command. 19 tests in `tests/test_bias_profile.py`, including
refusal of empty, single-class, non-binary and column-missing input, determinism,
and refusal to write under a `splits` directory. Artefact:
`artefacts/bias/iedb_tcell.bias.json`, sha256 `49044635…` after the schema
correction below; the file records the parquet's sha256 `38388c89…` as its input.

**Verified.** `repro-verifier` recomputed H(y), five mutual-information values and
every per-allele imbalance count independently (pandas, numpy, scikit-learn, not
`nbx.bias`) and agreed to about 1e-15 and exactly on counts. It also found one
misleading element in my first report, fixed before commit: "alleles whose mask
differs from the uniform mask" was 61 of 61, but 44 of those masks were empty and an
empty mask never equals the uniform one. That field is removed; the report now gives
alleles with a non-empty mask, with a mask equal to uniform, and with positions 2 and
the C-terminus masked, plus the scope as a fraction of the corpus.

**Measured** (46,608 records, 146 alleles; commands above):
- I(b;y) as a share of H(y): hla_allele 0.2006, assay_type 0.2241,
  c_terminal_residue 0.0565, position_2_residue 0.0498, peptide_length 0.0138,
  source_dataset 0 by construction (one source). No variable exceeds 25% of H(y), so
  the candidate list is empty.
- 124 alleles have a defined rho; 22 single-sample alleles are excluded; 48 of the
  124 lie beyond 10:1 and 30 are saturated.
- Masking policy, 9-mers only, alleles with at least 20 nine-mer records: 61 of 146
  alleles and 27,923 of 46,608 records. Only 17 of the 61 masks are non-empty,
  none equals the uniform {P2, C-terminus} mask, and two alleles mask position 2 and
  two mask the C-terminus. tau is the 90th percentile of pooled per-position MI,
  chosen on this same corpus, so by construction at most about a tenth of positions
  can exceed it.
- Split admissibility at the default criterion (`max_abs_rho` 1.0, minimum 20
  samples): 67 alleles admitted, 79 reported separately (77 below the sample
  minimum, 2 with a single class). Computed in memory; nothing was frozen.

**Does this corpus reproduce the published pattern?** Partly, and it cannot test the
rest. Zhang et al. report that more than half of alleles show a within-allele ratio
beyond tenfold (as summarised in the project synopsis, section 2.3). Here 48 of 124
alleles (about 39%) do, on a different and filtered corpus, so the published figure is
not reproduced at that threshold. Substantial imbalance is present. The claim that
the same allele has contrasting ratios across datasets cannot be tested on one
source. H1 and H2 need trained models and a second corpus; neither is evaluated.

**Cautions.**
- Plug-in MI is inflated for small groups. The 17 non-empty masks all come from
  alleles with 22 to 53 nine-mer records, so the positions that clear tau are as
  likely to be small-sample noise as anchor structure. This is an inference from
  group sizes, not a tested result. `assay_type` has 314 distinct values, 233 with
  fewer than five records, so its share of H(y) is inflated too.
- The whole profile inherits the undecided duplicate-label-conflict policy.

**Not done.** No split frozen. `nbx split freeze` still not wired. No bias
estimator correction chosen. No second corpus. `artefacts/` is gitignored, so the
profile is not version-controlled; its identity rests on the hashes above.

**Needs your decision before splits are frozen:** the conflict policy, because it
changes which records, and so which alleles, enter the split.

---

## 2026-10-11 — pre-registration and corpus v2 (pair-level labels)

**Decided.** The analysis plan is committed as `docs/analysis_plan.md` before any
model exists: hypotheses and refutation conditions, pair as the unit of analysis,
`majority` as the primary label rule with `strict` and `any_positive` as declared
sensitivity analyses, peptide-grouped splits, an allele-clustered paired bootstrap
for the margin over the HLA-only control, and a label-permutation null for the
masking policy. One item is still pending explicit approval: the shortcut-dependence
tolerance, because it changes the gate rule.

**Built.** `label_rule` on the harmonisation criterion (part of the content hash);
pair columns `n_assay_records`, `n_positive_assay_records`, `positive_fraction`;
record-level and pair-level reconciliation that raises if either fails; a distinct
drop reason for alleles that resolve to a two-field name but are absent from the
pseudo-sequence table (it was mislabelled before); a tamper check that re-hashes an
existing parquet before refreshing its report; `--label-rule` on `nbx data harmonise`;
`scripts/check_label_rule_artefacts.py`. 139 tests pass (35 in `test_harmonise.py`).

**Measured** (command: `nbx data harmonise ... --label-rule <rule>`):

| | majority | strict | any_positive |
|---|---|---|---|
| pairs | 49,395 | 46,608 | 50,782 |
| alleles | 147 | 146 | 152 |
| positive / negative | 13,107 / 36,288 | 10,925 / 35,683 | 15,099 / 35,683 |

50,782 pairs reached the rule: 10,925 unanimous positive, 35,683 unanimous negative,
4,174 mixed, of which 1,387 are ties. The `strict` artefact has exactly the
(peptide, allele, label) set of the earlier artefact (regression guard).

**Verified.** `repro-verifier` re-ran all three commands (hashes and parquet bytes
identical), recomputed every count from the raw export with independent code
(agreement exact), checked majority and any-positive semantics on random mixed pairs,
confirmed both reconciliation identities, and confirmed the tamper check raises.

**Caveats found.**
- Majority is over assay rows, not independent studies. 692 of 2,787 majority-decided
  pairs have all rows in one PMID group; a per-paper majority differs on 23 pairs and
  ties on 309. 22,223 pairs have no PMID.
- `records_without_reference_pmid` actually counts pairs; the name is misleading.
- Drop-reason names for two alleles reasons changed meaning, so drop tables from
  before and after this stage are not directly comparable.
- `artefacts/bias/iedb_tcell.bias.json` was computed on the pre-rule artefact and has
  not been re-run on the primary `majority` artefact.

**New findings that need your decision.**
1. *Assay-type scope.* The IEDB export contains 314 distinct assay types, including
   pMHC-multimer binding rows recorded as "qualitative binding", which are not T cell
   functional responses. `assay_type` carries the most information about the label of
   the profiled variables. Whether to restrict to functional immunogenicity assays is
   a scoping decision that changes the corpus; it is not made here.
2. *The ImmuBPI immunogenicity corpus is usable as a second source.* Downloaded from
   Zenodo 18610279 (md5 matches the record; download only, nothing vendored).
   `deepneo_tcr_10fold.csv` has 19,933 rows, 18,866 pairs, 136 alleles, 5,426
   positive and 14,507 negative, with a `source` column. Only 19 pairs carry
   conflicting labels (our raw export has 4,174), so it is curated differently. The
   positive fraction differs strongly by source (for example 0.038 for Bobisse, 0.167
   for IEDB, 0.514 for IMMA2, 0.946 for Calis) and 76 of the 136 alleles occur in more
   than one source, so the cross-source contrast of imbalance can be tested within it.
   Its allele names use another convention (`HLA-A-0201`). Licence still unresolved:
   GitHub says MIT, Zenodo says Apache-2.0.

**Housekeeping.** `PROMPTS.md` has a one-line change (a blank line) that this work did
not make; it is not included in these commits.

**Not done.** No split frozen; peptide-grouped splitting not yet built; the second
corpus is not harmonised; no model; M6 orchestration not started.

---

## 2026-10-11 (later) — peptide-grouped splitter and M2 on all three corpora

**Built.** `build_split(..., groups=peptide)` in `nbx.splits.freeze`: no peptide in more
than one of train, validation and test; records withheld to guarantee that are
quarantined and counted, never moved. `SplitCriterion.group_by` is part of the hash
when set; ungrouped hashes are byte-identical to before (pinned in tests, and
confirmed against the original code in 12 cases). 23 tests in
`tests/test_splits_grouped.py`; `scripts/measure_grouped_split.py` measures without
freezing. `nbx bias profile` now computes split admissibility with the grouped
splitter. 163 tests pass.

**Measured** on `iedb_tcell_majority` (49,395 pairs): ungrouped splitting puts 4,299
peptides (16,877 records) in more than one partition and 1,958 test peptides also in
train or validation; grouped splitting puts 0. Grouped: train 30,944, validation
6,631, test 6,475, quarantined 5,345 (959 positive, 4,386 negative); 69 alleles
admitted, 78 reported separately; the number admitted is unchanged by grouping.
`repro-verifier` recomputed these from the index tuples, found no counterexample in
400 random adversarial corpora, and confirmed determinism across `PYTHONHASHSEED`
values.

**Caveats.**
- The smallest admitted allele has 3 test records; nothing enforces a minimum
  test-block size (the 20-sample minimum applies to the allele's corpus total).
- Quarantine is label-skewed, so the train and validation positive rate rises
  relative to the corpus.
- Only exact-peptide leakage is removed; near-identical peptides are not.
- Group labels are not validated: NaN or None groups, or `1` versus `"1"`, are not
  rejected. Real peptides have none of these (30,921 unique raw and after strip and
  upper, no nulls), but the builder should raise on them.
- On the real, pair-level corpus every (allele, peptide) cell holds one record, so the
  whole-group selection path is exercised only by the synthetic tests.

**M2 on all three rules** (README table). The strict-rule run reproduces the earlier
figures exactly (48 of 124 alleles beyond 10:1). Primary `majority`: 45 of 124.

**Decisions still needed from you.** (1) Whether to restrict to functional T cell
assay types; (2) whether to enforce a minimum test-block size per admitted allele
(this changes the split criterion and its hash); (3) approval to freeze a split;
(4) approval of the shortcut-dependence tolerance (gate rule change); (5) approval
to push.

## 2026-10-11 (end of day) — M6 evaluation harness and gate redesign

**Built.** `nbx.evaluate.protocol` (run directory contract, refusal conditions,
content-hashed report), `nbx.evaluate.stats` (exact expectation of average precision
over random orderings within ties; weighted form for the bootstrap), a redesigned
`nbx.evaluate.gate.evaluate_gate`, and `nbx eval` (options `--run`, `--split`,
`--protocol`, `--out`, `--records-digest`; the bootstrap size and seed are not
options). Tests were written first for the gate truth table; `tests/test_protocol.py`
covers the harness. `python scripts/verify.py` passes (251 tests).

**Why the gate changed.** An independent adversarial audit of the first harness found
that it could rank a model that scores alleles by their training label rate, because
the kappa condition compared the model with the baseline's point estimate and the
baseline was equally biased. The rule now requires a significant margin over the
HLA-only control, an upper 95% bound on the change in |kappa| that does not exceed
zero, and no significant deficit against the peptide-only control. Details and the
withdrawn pre-registration clause are in the deviation log of `docs/analysis_plan.md`.
A first draft of the plan called the tie handling "equivalent" to threshold-grouped
average precision; that was wrong (they differ on tied blocks) and is corrected there.

**Verified** by `repro-verifier` on a fresh copy: an exhaustive grid of gate inputs
with no violation of the stated rule; the exact tie expectation against brute-force
enumeration of orderings; the bootstrap against an independent implementation; every
refusal case refuses; the three synthetic scenarios (allele-only model not
significant, allele shortcut plus peptide signal failed, peptide-driven low-kappa
model passed).

**Not done.**
- Bootstrap calibration at the real admitted-allele count with unequal allele sizes,
  and a mutation test of the suite. The audit agent that was running them was
  terminated by a rate limit; no result from it is used.
- An absolute ceiling on |kappa| (the value is a scientific decision).
- FR-65 and the published-model acceptance reproduction.
- Matched-configuration fields remain declared metadata. A coordinated forgery of
  `split.json` is refused only when the digest is pinned with `--records-digest`.
- No split has been frozen; no model has been trained.

**Decisions still needed from you.** The |kappa| ceiling; ratification of filing the
peptide-only condition under `FAILED_SHORTCUT`, of the neutral tie definition, of
`ranking_k` = 20 and 100, and of the MCC rule; the assay-scope list; the minimum test
records per admitted allele; which corpus is primary; the H2 operationalisation; and
approval to freeze splits and to install `torch`.

---

## 2026-10-26 — Stage 1: the evaluation hardened

**Built.**
- Gate condition (iii), decision R1: the upper 95% bound of |kappa_model| minus
  |kappa_HLA-only control| must be strictly below zero. Added to conditions (ii) and
  (iv), so the rule only tightens; a missing interval gives `PENDING`. Truth table in
  `tests/test_gate.py` first; wired into `nbx.evaluate.protocol` and `nbx eval`.
- Hooks run with the virtual-environment interpreter (R10).
- `scripts/calibrate_margin_test.py` and `scripts/mutation_test.py`.
- Decisions R1 to R14 recorded in `docs/analysis_plan.md`.

**Measured** (details and sources in `docs/results/calibration.md` and `mutation.md`).
- Calibration at the real regime (69 admitted alleles, 6,475 test records), 1,500
  replicates per scenario: false-significance rate 0.0513 (S1, upper 0.0637) and zero
  (S2, upper 0.0025). The R3 rule (upper at most 0.07) is met and the test is unchanged.
  Supplementary S3: 0.0567, upper 0.0696, inside the rule by 0.0004.
- Mutation test: 19 defects, first run 16 killed and 3 survived (M04, M05, M13); three
  tests added; rerun 19 killed.
- `python scripts/verify.py` passes (264 tests collected, lint, honesty and types clean).

**Not done / limits.** The 100-replicate B = 10,000 spot check of S1 (7 of 100) cannot meet
the rule at that size and is not evidence either way. S2 is conservative. Mutation testing
covers 19 hand-chosen defects in the evaluation module only. No data were frozen and no
model was trained. The roadmap tripwire dates for Stages 3 to 5 have passed and are
revisited when Stage 3 is planned. FR-65 and the published-model reproduction remain open.

