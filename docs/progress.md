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
