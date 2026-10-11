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
