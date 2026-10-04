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
