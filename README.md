# nbx

A bias-controlled evaluation protocol and computational pipeline for neoantigen
immunogenicity prediction.

Undergraduate (B.Tech., Computer Science and Engineering) major project, Babu
Banarasi Das University, Lucknow. Author: Pundarikaksh Narayan Tripathi.
Supervisor: Mr. Sachin Bhardwaj.

Status: research code under active development. The data, bias-diagnostic and
evaluation layers are built and tested; no model has been trained and no
immunogenicity result has been produced. See [State](#state).

---

## Motivation

For a mutated peptide to be a useful neoantigen it must be bound by the patient's
HLA class I molecule, presented at the cell surface, and recognised by a T cell
receptor. Predictors for the first two stages are reported to perform well
(presentation AUROC above 0.97 in the literature). Immunogenicity prediction does
not: in the TESLA consortium benchmark only about 6% of computationally
predicted epitopes were experimentally confirmed (Wells et al., Cell 2020).

Zhang et al. (Cell Genomics, April 2026) argue that this is not only a data
scarcity problem. Published immunogenicity models exploit intra-HLA label
imbalance: within a single allele the ratio of positive to negative peptides
frequently differs by more than tenfold, and the same allele shows different
ratios in different datasets, so the ratio reflects experimental design rather than
biology. They report that a control model given only the HLA allele, with the
peptide withheld, outperformed fifteen previously evaluated models on the IEDB
benchmark by 3.18 to 15.95 percent AUROC. A model that cannot represent T cell
recognition in principle ranks above models built to do so, which indicates that
the benchmark does not measure the intended task.

This repository addresses that evaluation problem and the pipeline needed to
study it.

## Contributions

The system is a six-layer pipeline from data ingestion to release. The
contributions are the parts that control for the shortcut:

| ID | Contribution | Module |
|---|---|---|
| O1 | Evaluation protocol: allele-balanced splits, allele-stratified metrics, mandatory matched HLA-only and peptide-only null controls, and an acceptance rule | M6, `nbx.evaluate`, `nbx.splits` |
| O2 | Information-theoretic bias diagnostics computed on the corpus before training | M2, `nbx.bias.diagnostics` |
| O3 | Adaptive per-allele masking derived from per-position mutual information, compared against no masking and uniform position-2/C-terminus masking | M5, `nbx.bias.policy` |
| O4 | Geometric interface features from an open structure predictor as a candidate bias-resistant input | M3, M4 |
| O5 | Reproducible open implementation | M7 |

The evaluation protocol (O1) is the primary claim, since it is the part that can be
defended within the project's time and compute budget. O3 and O4 are hypotheses
(H3, H4) that may be refuted; a negative result is reported as such. The
predictor is not claimed to be state of the art.

Bias diagnostics (L2) sit upstream of modelling (L4), so that the measured bias of a
corpus conditions training instead of annotating results afterwards.

### Quantities

Defined once in `src/nbx/bias/diagnostics.py` and tested against analytic cases and
against `scipy` and `scikit-learn`:

```
I(b;y)  = H(y) - H(y|b)               mutual information between a bias candidate b and the label y
rho(a)  = log2(N+(a) / N-(a))         per-allele log imbalance ratio, saturated at +/-10
kappa(f)= spearman(rho(a), E[f(x)|a]) shortcut dependence of a trained model f
M(a)    = { j : I(x_j; y | a) > tau } per-allele masking policy
```

A model that has learned peptide-level determinants should have kappa near zero. A
clearly positive kappa indicates that the model scores an allele by how that allele
was labelled in training.

### Acceptance rule

Implemented in `src/nbx/evaluate/gate.py` and specified by the truth table in
`tests/test_gate.py`. The burden of proof is on the model. A variant is rankable only
if all three conditions hold:

1. its margin over the matched HLA-only control, in the headline metric, is
   significant (one-sided, allele-clustered bootstrap);
2. the upper 95% bound of the change in |kappa| relative to the baseline does not
   exceed a non-inferiority margin of zero;
3. it is not significantly worse than the matched peptide-only control.

A variant is in one of five states:

| State | Meaning |
|---|---|
| `PENDING` | A mandated null control or a required input is missing. No metric is displayed. |
| `NOT_SIGNIFICANT` | Controls present; margin over the HLA-only control not significant. No rank. |
| `FAILED_SHORTCUT` | Margin significant but condition 2 or 3 failed; the failed condition is named. Reported, not ranked. |
| `PASSED` | All three conditions hold. Ranked, margin shown. |
| `CONTROL` | The row is itself a null control. Never ranked. |

An earlier version compared a model's kappa with the baseline's point estimate. An
adversarial audit showed that a model scoring alleles purely by their training label
rate passed whenever its baseline was equally biased. The rule was replaced; the
change and its reasons are in the deviation log of `docs/analysis_plan.md`. An
absolute ceiling on |kappa| is not implemented because its value is a scientific
judgement that has not been made; until it is set, a model exactly as shortcut
dependent as its baseline, with a zero-width interval, can satisfy condition 2.

The gate has no override parameter by design. Quantities that have not been
measured render as an em dash and never as zero.

---

## State

`nbx status` reports the current state; this table may lag it.

| Module | Package | State |
|---|---|---|
| M1 Data harmonisation | `nbx.data.harmonise` | IEDB T cell export only; MHC ligand, TESLA and HiTIDE not ingested |
| M2 Bias diagnostics | `nbx.bias.diagnostics`, `nbx.bias.profile` | built, tested; run once on the IEDB T cell corpus (see below) |
| M3 Representation service | `nbx.represent.cache` | contract only |
| M4 Model core | `nbx.model.backbone` | contract only |
| M5 Masking policy | `nbx.bias.policy` | built, tested |
| M6 Evaluation harness | `nbx.evaluate` | built and tested on synthetic predictions (gate, metrics, tie-exact AUPRC, clustered bootstrap, `nbx eval`); never run on real model output; FR-65 and the published-model reproduction not done |
| M7 Release packaging | `nbx.release` | contract only |

Hypotheses H1 to H4 have not been evaluated. Until they are, results are shown as
an em dash.

### M1: IEDB T cell corpus

Input: the IEDB `tcell_full_v3` full export (CSV, two header rows, downloaded
2026-10-05; checksums in `data/raw/MANIFEST.txt`, which is gitignored along with the
raw data). Allele pseudo-sequences (34 residues, NetMHCpan convention) are read
from `IEDB_MHC_pseudo.dat` in the ImmuBPI repository; the file's origin and licence
are not stated upstream, so it is downloaded locally and not vendored.

The unit of analysis is the (peptide, allele) pair; raw assay rows are repeated
measurements of a pair. Retained pairs satisfy: human HLA class I with a resolvable
two-field allele present in the pseudo-sequence table, linear unmodified peptide of 8
to 14 standard residues, qualitative outcome mapped to a binary label (all `Positive*`
levels to 1, `Negative` to 0), and non-empty source and assay method. A record
missing provenance is dropped, never defaulted. Pairs whose assays disagree are
labelled by one of three rules, fixed in advance in `docs/analysis_plan.md`:

| Rule | Mixed pairs | Role |
|---|---|---|
| `majority` | label by strict majority of assay records; ties dropped | primary |
| `strict` | every disagreeing pair dropped | sensitivity analysis |
| `any_positive` | positive if any assay record is positive | sensitivity analysis |

Output of the runs on 2026-10-11 (from `artefacts/datasets/iedb_tcell_<rule>.report.json`,
reproducible with the command below):

| Quantity | majority | strict | any_positive |
|---|---|---|---|
| Pairs | 49,395 | 46,608 | 50,782 |
| Alleles | 147 | 146 | 152 |
| Positive / negative | 13,107 / 36,288 | 10,925 / 35,683 | 15,099 / 35,683 |
| Pairs dropped by the rule | 1,387 (ties) | 4,174 | 0 |
| Assay records behind kept pairs | 94,491 | 64,476 | 98,121 |

The raw export has 578,203 assay rows. For every rule, rows kept, dropped (each under
one named reason) and collapsed into pairs sum to that count; the module raises if
they do not. Content hashes are recorded in each report. A pair counts as mixed in
4,174 cases, of which 1,387 are ties. The three class balances above are a consequence
of the labelling rule, not a property of T cell biology.

Caveats that bear on interpretation:

- A majority is taken over assay rows, not over independent studies. 692 of the 2,787
  pairs decided by majority have all their rows in a single PMID group, and 22,223
  pairs carry no PMID. A per-paper majority would agree with the row majority on all
  but 23 of those 2,787 pairs, with a further 309 per-paper ties.
- The export mixes assay types (it has 314 distinct values) whose meaning differs, for
  example pMHC-multimer binding and functional T cell assays. `assay_type` carries the
  most information about the label of the variables profiled, and whether all of these
  assays measure immunogenicity is a scoping question that is open.
- `HLA-B*44:01` is absent from the pseudo-sequence table and its 4 records are
  excluded.
- No published reference statistic for this subset has been located, so the
  reproduction check (FR-06) reports no reference, rendered as an em dash, which is
  not a match.

### M2: bias profile of the IEDB T cell corpus

```bash
nbx bias profile artefacts/datasets/iedb_tcell_majority.parquet --report artefacts/bias/iedb_tcell_majority.bias.json
```

The report records the input file's sha256 and is byte-for-byte reproducible. Run on
2026-10-11 for each label rule:

| Quantity | majority | strict | any_positive |
|---|---|---|---|
| Alleles with a defined imbalance ratio | 124 | 124 | 126 |
| Of those, beyond 10:1 in either direction | 45 | 48 | 48 |
| I(allele; y) as a share of H(y) | 0.197 | 0.201 | 0.199 |
| I(assay type; y) as a share of H(y) | 0.261 | 0.224 | 0.272 |
| Masking-policy scope (alleles / records) | 61 / 29,872 | 61 / 27,923 | 61 / 30,737 |
| Alleles with a non-empty adaptive mask | 18 | 17 | 20 |
| Split: admitted / reported separately | 69 / 78 | 67 / 79 | 70 / 82 |
| Split: test records / quarantined records | 6,475 / 5,345 | 5,992 / 5,064 | 6,718 / 5,610 |

- The published claim is that more than half of alleles exceed a tenfold ratio (Zhang
  et al. 2026, as summarised in the project synopsis). Here it is 45 of 124 for the
  primary rule, on a filtered single-source corpus, so that figure is not reproduced.
  Contrasting ratios across datasets cannot be tested with one source.
- Mutual information with the label is a bias candidate signal, not a verdict. No
  variable exceeds a quarter of the label entropy for `majority` except `assay_type`
  (0.261), whose value is inflated because most of its 314 values are rare, and which
  mixes assays of different meaning. `source_dataset` is zero by construction.
- The adaptive masking policy covers 9-mers for 61 of the 147 alleles. Its threshold is
  chosen on the corpus being profiled, and plug-in mutual information is inflated for
  small allele subsets, so the policy block is exploratory.
- The split block is computed in memory with peptide grouping (all records of a
  peptide in one partition). Records withheld to guarantee that are quarantined and
  counted, never moved. No split has been frozen.
- The smallest admitted allele has 3 test records, so per-allele test statistics can be
  very noisy; nothing yet enforces a minimum test-block size.
- No hypothesis (H1 to H4) has been evaluated.

#### Why splits are grouped by peptide

In the `majority` corpus 6,832 of 30,921 peptides occur under more than one allele
(25,306 of 49,395 records). A record-level split shuffles within allele and does not
see the peptide: measured on this corpus it places 4,299 peptides (16,877 records) in
more than one of train, validation and test, and 1,958 test peptides also occur in
train or validation. That lets any peptide-aware model, and the peptide-only control,
score test records by recognising training peptides. The grouped split has none
(`scripts/measure_grouped_split.py`). Its cost is 5,345 records (10.8%) withheld from
train and validation, of which 959 are positive and 4,386 negative, so the withheld
set is not label-neutral. Near-identical but non-identical peptides can still cross
partitions; only exact peptide leakage is removed.

### M6: evaluation harness

`nbx.evaluate.protocol` is model-agnostic. It consumes the predictions of runs, not
models, and produces the evaluation. Inputs are a frozen split (`split.json`, with the
per-record peptide, allele and label and a digest of them) and, for each variant, a
`config.json` and one `seed_<n>.parquet` of test-set scores per seed. Statistics:

- **Metrics.** AUPRC (headline), AUROC, PPV@k and Recall@k, MCC, pooled and per
  allele, with across-seed mean and standard deviation. Ties in the scores are handled
  by the exact expectation of average precision over uniformly random orderings within
  tied groups, so tie handling has no seed (`nbx.evaluate.stats`).
- **Margin test.** Paired, allele-clustered bootstrap of the seed-averaged AUPRC
  difference against the HLA-only control, one-sided, alpha 0.05, 10,000 resamples,
  fixed bootstrap seed. The one-sided p-value is reported with the matching 90%
  two-sided interval and the one-sided 95% lower bound. A Holm-adjusted p-value and the
  number of variants compared are reported; the gate reads the unadjusted p-value.
- **Pinned settings.** Trained-model seeds must be exactly 1 to 5. The command line
  cannot lower the resample count or change the bootstrap seed.

The harness refuses, with a non-zero exit and no metric printed, when: the split hash
or records digest does not recompute; a control is missing; seeds or matched
configuration (epochs, parameter count, early-stopping metric, tuning trials) differ;
a score is missing, duplicated or non-finite; the test set has a single class; a
peptide occurs in more than one partition; an HLA-only control is not constant within
each allele; or a peptide-only control scores two records of the same peptide
differently. The report labels control rows `CONTROL`, flags the baseline
`reference_only`, and shows the margin beside every model metric.

Tested on synthetic predictions only (251 tests pass under `python scripts/verify.py`;
the count is that of the run on 2026-10-11). Three synthetic scenarios pin the
behaviour that motivated the design: a model that scores only by allele is
`NOT_SIGNIFICANT`; a model that adds a peptide signal to an allele shortcut is
`FAILED_SHORTCUT`; a model driven by the peptide with low kappa is `PASSED`. These
show the machinery distinguishes the cases; they say nothing about any real model.

Known limits, also stated in the harness report:

- Matched epochs, parameters, early-stopping metric and tuning trials are declared
  metadata, not verified.
- A coordinated edit of `split.json` (alter a label and recompute the digest) is
  detected only when the digest is pinned with `--records-digest`.
- The baseline's identity is read from the run directory, not pinned in the analysis
  plan.
- The bootstrap's behaviour at the real admitted-allele count and unequal allele sizes
  has not been calibrated, and a mutation test of the suite has not been run. Both are
  recorded as not done in `docs/analysis_plan.md`.
- FR-65 (control sets with inverted intra-allelic imbalance) and the acceptance
  reproduction of a published model on a published benchmark are not done.

---

## Installation and verification

Python 3.10 or later.

```bash
git clone https://github.com/PundarikakshNTripathi/bias-controlled-neoantigen-immunogenicity-pipeline.git
cd bias-controlled-neoantigen-immunogenicity-pipeline

python -m venv .venv
# Windows:      .venv\Scripts\activate
# Linux/macOS:  source .venv/bin/activate

pip install -e ".[dev]"

python scripts/verify.py      # lint, honesty check, tests, types
nbx doctor                    # integrity core self-check
nbx status
```

The virtual environment must be active when the repository's Claude Code hooks run,
because they invoke `python`.

### Reproducing the harmonised corpus

Download the IEDB `tcell_full_v3.zip` export and `IEDB_MHC_pseudo.dat` into
`data/raw/` (not committed), then:

```bash
nbx data harmonise \
  --iedb-tcell data/raw/iedb/tcell_full_v3.zip \
  --mhc-pseudo data/raw/immubpi/IEDB_MHC_pseudo.dat \
  --out artefacts/datasets \
  --label-rule majority          # or strict, or any_positive
```

This writes `iedb_tcell_<rule>.parquet` and `iedb_tcell_<rule>.report.json`. The content hash covers
both the data and the filter criterion, so a changed criterion produces a different
hash, and the command refuses to overwrite an existing artefact whose hash differs.
IEDB data are released under CC BY 4.0; cite IEDB when using derived data
(Vita et al., Nucleic Acids Research 47(D1), 2019).

### Commands

| Command | Function | State |
|---|---|---|
| `nbx data harmonise` | Ingest a public corpus into one canonical dataset | built (IEDB T cell) |
| `nbx bias profile` | Compute I(b;y), rho(a), per-position MI and the masking policy | built |
| `nbx split freeze` | Build and freeze an allele-balanced split | CLI not wired |
| `nbx repr cache` | Compute and cache representations | not built |
| `nbx train --matrix` | Run the ablation matrix with mandatory null controls | not built |
| `nbx eval --split <hash>` | Evaluate prediction runs under the protocol | built; tested on synthetic runs only |
| `nbx report` | Package the benchmark report | not built |
| `nbx doctor`, `nbx status` | Integrity check, module state | built |

Commands for unbuilt modules exit with code 2 and name the module they depend on.

---

## Reproducibility rules

These are enforced by tests and by hooks in `.claude/`:

1. No metric value, count or p-value appears in code or documentation unless a
   command in this repository computed it or it is cited to a source with a
   locator.
2. Unmeasured quantities render as an em dash.
3. The acceptance gate cannot be bypassed.
4. Frozen splits under `artefacts/splits/` are immutable. A different criterion is a
   different experiment with a different content hash.
5. Every trained variant ships with matched HLA-only and peptide-only controls at
   equal capacity, training budget and tuning effort.
6. Negative results are reported. Hyperparameters are not tuned until a hypothesis
   passes.
7. Results are reported as means with across-seed variance, not as best runs.

`CLAUDE.md` holds the full project rules and `PROMPTS.md` the prompts used to drive
development with Claude Code. `SETUP.md` covers environment and tooling setup.

---

## Scope and limitations

Computational and machine-learning work only, on public de-identified data: IEDB,
TESLA-derived evaluation sets, HiTIDE/NCI harmonised corpora, PDB and UniProt.
The project scope is HLA class I presentation to CD8-positive T cells; class II is
excluded.

Conclusions are bounded by the quality, coverage and heterogeneity of existing
public assay data; no new experimental ground truth is generated. Out of scope:
patient-identifiable or controlled-access data, wet-laboratory work, and any
clinical or diagnostic claim. This is research software, not a medical device, and
candidate-ranking output carries a permanent research-use-only notice.

## References

- Zhang Z., Wei L., Wu W., Qi H., Wang X. Cross-task interpretability through
  unified modeling reveals a universal shortcut bias in neoantigen prediction.
  *Cell Genomics* 6(7), 101214 (2026). doi:10.1016/j.xgen.2026.101214. Code:
  [WangLabTHU/ImmuBPI](https://github.com/WangLabTHU/ImmuBPI) (repository states MIT);
  data: [Zenodo 18610279](https://zenodo.org/records/18610279) (record lists
  Apache-2.0). The two licence statements differ and have not yet been reconciled.
- Wells D.K. et al. Key parameters of tumor epitope immunogenicity revealed through
  a consortium approach improve neoantigen prediction. *Cell* 183(3), 818-834 (2020).
  doi:10.1016/j.cell.2020.09.015.
- Reynisson B. et al. NetMHCpan-4.1 and NetMHCIIpan-4.0. *Nucleic Acids Research*
  48(W1), W449-W454 (2020). doi:10.1093/nar/gkaa379.
- Albert B.A. et al. Deep neural networks predict class I major histocompatibility
  complex epitope presentation and transfer learn neoantigen immunogenicity. *Nature
  Machine Intelligence* 5, 861-872 (2023). doi:10.1038/s42256-023-00694-6.
- Vita R. et al. The Immune Epitope Database (IEDB): 2018 update. *Nucleic Acids
  Research* 47(D1), D339-D343 (2019). doi:10.1093/nar/gky1006.

## Licence

MIT; see [LICENSE](LICENSE). Third-party models, datasets and baselines remain under
their own terms.
