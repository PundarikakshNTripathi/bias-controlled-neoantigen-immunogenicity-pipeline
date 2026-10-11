# nbx

A bias-controlled evaluation protocol and computational pipeline for neoantigen
immunogenicity prediction.

Undergraduate (B.Tech., Computer Science and Engineering) major project, Babu
Banarasi Das University, Lucknow. Author: Pundarikaksh Narayan Tripathi.
Supervisor: Mr. Sachin Bhardwaj.

Status: research code under active development. No model has been trained and no
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
`tests/test_gate.py`. A result is positive only if discrimination improves over the
matched HLA-only control by a significant margin and shortcut dependence does not
increase. A variant is in one of five states:

| State | Meaning |
|---|---|
| `PENDING` | A mandated null control or a required metric is missing. No metric is displayed. |
| `NOT_SIGNIFICANT` | Controls present; margin over the HLA-only control not significant. No rank. |
| `FAILED_SHORTCUT` | Margin significant but kappa increased. Reported, not ranked. |
| `PASSED` | Discrimination improved and kappa did not increase. Ranked, margin shown. |
| `CONTROL` | The row is itself a null control. Never ranked. |

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
| M6 Evaluation harness | `nbx.evaluate` | gate and metrics built; orchestration (`protocol.py`) not built |
| M7 Release packaging | `nbx.release` | contract only |

Hypotheses H1 to H4 have not been evaluated. Until they are, results are shown as
an em dash.

### M1: IEDB T cell corpus

Input: the IEDB `tcell_full_v3` full export (CSV, two header rows, downloaded
2026-10-05; checksums in `data/raw/MANIFEST.txt`, which is gitignored along with the
raw data). Allele pseudo-sequences (34 residues, NetMHCpan convention) are read
from `IEDB_MHC_pseudo.dat` in the ImmuBPI repository; the file's origin and licence
are not stated upstream, so it is downloaded locally and not vendored.

Retained records satisfy: human HLA class I with a resolvable two-field allele,
linear unmodified peptide of 8 to 14 standard residues, qualitative outcome
mapped to a binary label (all `Positive*` levels to 1, `Negative` to 0), and
non-empty source, assay method and negative-sample origin. A record missing
provenance is dropped, never defaulted. Records duplicated on (peptide, allele) are
merged with the union of their provenance; pairs with contradictory labels are
dropped and counted.

Output of the run on 2026-10-05 (from `artefacts/datasets/iedb_tcell.report.json`,
reproducible with the command below):

| Quantity | Value |
|---|---|
| Raw assay rows | 578,203 |
| Canonical records | 46,608 |
| Alleles | 146 |
| Distinct peptides | 28,714 |
| Positive / negative | 10,925 / 35,683 |
| Content hash | `3cc846db315c0c2c9ddfc2af232638cb0dac91d29fdba59bb5cfed77029a7821` |

Kept, dropped and merged records sum to the raw row count; every drop is assigned
to exactly one named reason in the report.

Two properties of this corpus need to be kept in view. First, dropping
contradictory-label pairs removes 25,034 positive and 8,611 negative records, so the
positive-to-negative ratio moves from roughly 1:1 before to roughly 1:3 after;
this single policy choice changes the label balance and is part of the hashed
criterion. Second, `HLA-B*44:01` is absent from the pseudo-sequence table and is
therefore excluded. No published reference statistic for this subset has been
located, so the reproduction check (FR-06) reports no reference rather than a
match.

### M2: bias profile of the IEDB T cell corpus

```bash
nbx bias profile artefacts/datasets/iedb_tcell.parquet --report artefacts/bias/iedb_tcell.bias.json
```

The report records the input file's sha256 and is byte-for-byte reproducible. Read
from the run of 2026-10-11 on the corpus above:

- Of 146 alleles, 124 have a defined imbalance ratio (22 single-sample alleles are
  excluded), and 48 of those 124 lie beyond 10:1 in either direction. The published
  claim is that more than half of alleles do (Zhang et al. 2026, as summarised in the
  project synopsis); on this filtered, single-source corpus the fraction is lower, so
  that figure is not reproduced here. Contrasting ratios across datasets cannot be
  tested with one source.
- Mutual information with the label is reported as a bias candidate signal, not a
  verdict. No variable exceeds a quarter of the label entropy. `source_dataset` is
  zero by construction because the corpus has one source, and the `assay_type`
  figure is inflated by many rare values.
- The adaptive masking policy was computed for 9-mers over 61 of the 146 alleles.
  Its threshold is chosen on the same corpus (no frozen validation split exists) and
  plug-in mutual information is inflated for small allele subsets, so the policy
  block is exploratory. No hypothesis (H1 to H4) has been evaluated.

Split admissibility is reported in memory at the default criterion; no split has been
frozen.

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
  --out artefacts/datasets
```

This writes `iedb_tcell.parquet` and `iedb_tcell.report.json`. The content hash covers
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
| `nbx eval --protocol bias-controlled` | Evaluate under the protocol | not built |
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
