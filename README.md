# nbx

**A bias-controlled benchmark protocol for neoantigen immunogenicity prediction.**

Final-year B.Tech. major project · Babu Banarasi Das University, Lucknow
Author: Pundarikaksh Narayan Tripathi

---

## The problem

A mutated peptide must clear three stages to be a useful neoantigen: bind the
patient's HLA molecule, be presented at the cell surface, and be recognised by a
T cell. Binding and presentation are, in engineering terms, handled.
Immunogenicity is not — the TESLA consortium confirmed only about **6%** of
computationally predicted epitopes as immunogenic.

Work published in *Cell Genomics* in April 2026 established why this is worse
than a data-scarcity problem. Deep-learning immunogenicity predictors exploit
**intra-HLA label imbalance** as a statistical shortcut. Within an allele, the
positive-to-negative ratio often differs by more than tenfold, and the *same*
allele shows opposite ratios in different datasets — which proves the ratio is
an artefact of experimental design rather than a property of the biology.

Then the result that defines this project: three input variants were trained —
peptide-only, peptide+HLA, and **HLA-only**, with the peptide withheld
entirely. The HLA-only model, which cannot represent T cell recognition even in
principle, performed **best**, exceeding fifteen previously evaluated models on
the IEDB benchmark by 3.18–15.95% AUROC.

> If a model structurally incapable of the task outranks models designed for it,
> the benchmark is not measuring the task.

## What this repository is

A **complete computational pipeline** for neoantigen immunogenicity prediction:
ingestion and harmonisation, bias quantification, representation, model,
evaluation, release. Six layers, L1 to L6.

What distinguishes it from the pipelines that already exist is structural. Bias
measurement sits **upstream** of training rather than downstream of results, and
the evaluation layer can **refuse** — a model that has not beaten its own null
control cannot be ranked, and a quantity that has not been measured cannot be
displayed as though it had.

### Is this a benchmark or a pipeline?

Both, and the distinction matters.

The **system** is a pipeline. A trained immunogenicity predictor comes out of the
end of it, together with the diagnostics and the debiasing method.

The **headline claim** is the evaluation protocol, because that is what the field
is currently missing and what this work can defend within one semester. Two of
the four research objectives are nevertheless modelling contributions, not
benchmark contributions: the adaptive per-allele masking scheme (O3) and the
geometric structural channel as a bias-resistant signal (O4).

So: the model is real and is built. It is simply not claimed to be state of the
art, because claiming that would require exactly the kind of benchmark this
project argues does not currently exist.

### The four quantities

| | | |
|---|---|---|
| `I(b;y) = H(y) − H(y\|b)` | mutual information | detects bias candidates |
| `ρ(a) = log₂(N⁺(a)/N⁻(a))` | per-allele imbalance | saturates at ±10 |
| `κ(f) = spearman(ρ(a), E[f(x)\|a])` | **shortcut dependence** | the headline number |
| `M(a) = {j : I(x_j;y\|a) > τ}` | adaptive masking policy | generalises fixed anchor masking |

κ near zero is what a model that learned biology should show. A clearly positive
κ means the model scores an allele according to how that allele was labelled in
training.

### The acceptance rule

Implemented in `nbx/evaluate/gate.py`, with a truth table in `tests/test_gate.py`:

> A result counts as positive only if discrimination improves **and** shortcut
> dependence does not rise.

Four reachable states — `PENDING` (a mandated null control is missing; no
headline metric may be displayed at all), `NOT_SIGNIFICANT`, `FAILED_SHORTCUT`,
`PASSED`. There is deliberately no `force`, `override` or `skip_controls`
parameter, and a hook blocks adding one.

---

## State

Run `nbx status` for the live answer. At the time of writing:

| Module | Package | State |
|---|---|---|
| M1 Data harmonisation | `nbx.data.harmonise` | contract only |
| **M2 Bias diagnostics** | `nbx.bias.diagnostics` | **built, tested** |
| M3 Representation service | `nbx.represent.cache` | contract only |
| M4 Model core | `nbx.model.backbone` | contract only |
| **M5 Masking policy** | `nbx.bias.policy` | **built, tested** |
| M6 Evaluation harness | `nbx.evaluate` | gate + metrics built; orchestration pending |
| M7 Release packaging | `nbx.release` | contract only |

No experiment has been run. Every metric in every interface therefore renders as
an em dash (—). That is the designed empty state, not a gap waiting for sample
values.

The mathematics in M2 and M5, the acceptance gate, the split construction and
the metrics are covered by tests that pin them against analytically known cases
and cross-check them against `scipy` and `scikit-learn`.

---

## Quick start

```bash
git clone https://github.com/PundarikakshNTripathi/bias-controlled-neoantigen-immunogenicity-pipeline.git
cd bias-controlled-neoantigen-immunogenicity-pipeline

python -m venv .venv
# Windows:  .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate

pip install -e ".[dev]"

python scripts/verify.py     # should print VERIFICATION PASSED
nbx doctor
nbx status
```

### Commands

| Command | Does |
|---|---|
| `nbx init` | Scaffold the artefact store |
| `nbx data harmonise` | Ingest public corpora into one canonical dataset |
| `nbx bias profile` | Compute I(b;y), ρ(a) and derive the masking policy |
| `nbx split freeze` | Build and freeze an allele-balanced split (R1) |
| `nbx repr cache` | Compute and cache representations |
| `nbx train --matrix` | Run the ablation matrix with mandated null controls |
| `nbx eval --protocol bias-controlled` | Evaluate under the protocol |
| `nbx report` | Package the benchmark report (deliverable D1) |
| `nbx doctor` / `nbx status` | Integrity check / module state |

Commands for unbuilt modules exit with code 2 and point at the prompt that
builds them. They do not pretend to work.

---

## Working with Claude Code

This repository ships a harness that enforces its own research-integrity rules.

```bash
claude            # from the repo root
```

**Slash commands**

| Command | Does |
|---|---|
| `/next` | Decide what to work on, from the real state of the repo |
| `/module M1` | Build one module test-first against its contract |
| `/verify` | Run the verification loop and show the real output |
| `/gate-check` | Apply the acceptance rule before a result is reported |
| `/bias-report` | Compute and present a corpus bias profile |

**Subagents** — `bias-auditor` (adversarial result audit), `repro-verifier`
(verify by execution), `literature-verifier` (check every citation),
`ml-engineer`, `data-steward`.

**Hooks** — block writes to frozen splits, block adding a gate override, block
unsourced metric values in result-bearing documents, run fast checks after every
Python edit, and inject real project state at session start.

The guards are deliberately strict. If one blocks something legitimate, that is
a conversation about changing the rule, not a reason to disable it.

See **[SETUP.md](SETUP.md)** to get running, **[PROMPTS.md](PROMPTS.md)** for the
prompt library, and **[CLAUDE.md](CLAUDE.md)** for the full project rules.

---

## Scope

Computational and machine-learning work only, on public de-identified data:
IEDB, TESLA-derived evaluation sets, HiTIDE/NCI harmonised corpora, PDB, UniProt.

**Out of scope:** patient-identifiable or controlled-access data, wet-laboratory
work of any kind, and clinical or diagnostic claims. This is research software.
It is not a medical device, is not validated for clinical use, and candidate
ranking output carries a permanent research-use-only notice.

A longer-term programme toward personalised therapeutic pipelines exists, and is
gated on institutional collaboration, licensing and regulatory approval. Nothing
in this repository depends on it.

---

## Upstream work

- Zhang et al., *Cross-task interpretability through unified modeling reveals a
  universal shortcut bias in neoantigen prediction*, **Cell Genomics**, April 2026.
  Code: [WangLabTHU/ImmuBPI](https://github.com/WangLabTHU/ImmuBPI) ·
  data: [Zenodo 18610279](https://zenodo.org/records/18610279).
  *(The repository states MIT; the Zenodo record lists Apache-2.0. Check which
  applies before vendoring anything.)*
- Wells et al., *Key parameters of tumor epitope immunogenicity revealed through
  a consortium approach*, **Cell** 183(3), 2020. (TESLA)
- Reynisson et al., *NetMHCpan-4.1*, **Nucleic Acids Research** 48(W1), 2020.
- Albert et al., *BigMHC*, **Nature Machine Intelligence** 5, 2023.

## Licence

MIT — see [LICENSE](LICENSE). Third-party models, datasets and baselines remain
under their own terms.
