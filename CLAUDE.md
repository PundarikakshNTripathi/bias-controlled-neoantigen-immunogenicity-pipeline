# nbx — project instructions

A bias-controlled benchmark protocol for neoantigen immunogenicity prediction.
Final-year B.Tech. major project, Babu Banarasi Das University. Author:
Pundarikaksh Narayan Tripathi.

## What this project is actually about

The field has a documented failure: deep-learning immunogenicity predictors
exploit intra-HLA label imbalance as a statistical shortcut instead of learning
peptide-level determinants of T cell recognition. An HLA-only control — which
never sees the peptide and therefore cannot model recognition even in principle
— **outperformed fifteen published models** on the IEDB benchmark.

If a model structurally incapable of the task wins, the benchmark is not
measuring the task.

What we are building is a **complete computational pipeline** — L1 ingestion
through L6 release — that produces a trained immunogenicity predictor. What is
*novel* about it is the bias control: the diagnostics (O2), the adaptive
per-allele masking (O3), the geometric bias-resistant channel (O4), and the
evaluation protocol (O1) that ties them together.

The protocol is the **headline claim**, because it is the defensible one and it
ships first. It is not the whole system. Do not describe this project as "just a
benchmark", and do not describe the model as state of the art — neither is true.

This matters for how you work here: the entire intellectual content of the
project is about not fooling yourself with favourable-looking numbers. Code that
overclaims contradicts the thesis.

## Non-negotiable rules

These are enforced by hooks and tests. Do not work around them; if one blocks
something legitimate, say so and we change the rule deliberately.

1. **Never invent a number.** No metric value, record count, allele count or
   p-value may appear in code, docs, commit messages or chat unless it was
   computed by a command in this repo or cited to a source with a locator. If
   you need an example value in a docstring, label it as illustrative.

2. **Unmeasured renders as an em dash (—).** Never `0`, never `0.0`, never
   `null`, never a plausible-looking stand-in. Zero is a measurement; absence is
   not. `nbx.evaluate.gate.render_metric` is the only correct way to render.

3. **The gate is not bypassable.** No result may be displayed without its
   null-control margin. No rank may be assigned to a variant that has not
   exceeded the HLA-only control without raising shortcut dependence. Do not add
   `force=`, `override=` or `skip_controls=` to anything.

4. **Frozen splits are immutable.** Files under `artefacts/splits/` are never
   edited. A different criterion is a different experiment with a different
   hash.

5. **Null controls are mandatory and matched.** Every trained variant ships with
   an HLA-only and a peptide-only control at equal capacity, equal training
   budget and equal tuning effort. They cannot be disabled.

6. **Negative results are results.** If adaptive masking does not beat uniform
   masking, that is a finding to report, not a bug to fix. Never tune until a
   hypothesis passes.

7. **Report variance, not best runs.** Means across seeds with variance. Never
   the best seed.

## Scope boundary

Computational and machine-learning work only, on public de-identified data
(IEDB, TESLA-derived sets, HiTIDE/NCI harmonised corpora, PDB, UniProt).

**Never** add: patient-identifiable data, controlled-access datasets (no dbGaP),
wet-lab protocols, clinical or diagnostic claims. Candidate-ranking output
carries a permanent research-use-only notice.

If a task seems to need any of the above, stop and say so.

## Architecture

Six layers. The ordering is the point: **L2 sits upstream of L4**, so bias
measurement conditions training rather than annotating results afterwards.

| Layer | Module | Package | State |
|---|---|---|---|
| L1 ingestion | M1 Data harmonisation | `nbx.data.harmonise` | built for IEDB T cell export only, tested; other corpora not ingested |
| L2 bias quantification | M2 Bias diagnostics | `nbx.bias.diagnostics` | **built, tested** |
| L3 representation | M3 Representation service | `nbx.represent.cache` | contract only |
| L4 modelling | M4 Model core | `nbx.model.backbone` | contract only |
| L4 debiasing | M5 Masking policy | `nbx.bias.policy` | **built, tested** |
| L5 evaluation | M6 Evaluation harness | `nbx.evaluate` | gate, metrics, protocol and `nbx eval` built; tested on synthetic predictions only |
| L6 release | M7 Packaging | `nbx.release` | contract only |

Run `nbx status` to see current state. Do not trust this table over the command.

## The four quantities

Defined in `nbx/bias/diagnostics.py`, pinned by tests against analytic cases and
cross-checked against scipy and sklearn. Do not reimplement them elsewhere.

```
I(b;y)   = H(y) − H(y|b)              mutual information, bias candidate detector
ρ(a)     = log2(N⁺(a) / N⁻(a))        per-allele imbalance; saturates at ±10
κ(f)     = spearman(ρ(a), E[f(x)|a])  shortcut dependence — the headline number
M(a)     = {j : I(x_j;y|a) > τ}       adaptive per-allele masking policy
```

κ near zero is what a model that learned biology should show. A clearly positive
κ means the model scores an allele by how that allele was labelled in training.

## Working style

- **Tests first for anything load-bearing.** The bias math, the gate and the
  splits already work this way. Keep it.
- **Small commits, one concern each.** Conventional commit prefixes.
- **Before claiming a module is done**, run `python scripts/verify.py` and paste
  the actual output. Not "tests pass" — the output.
- **When reproducing a published result**, state the published number, its
  source, and your reproduced number. If they differ, document the discrepancy
  and evaluate against your own reproduction, clearly labelled.
- **Ask before adding a dependency.** The compute and the complexity budget are
  both small and fixed.

## Roadmap workflow

The remaining work is eight stages in `docs/ROADMAP.md`, which also holds the state. The
researcher types `/next`: plan mode for the current stage, then on approval build, verify,
document, commit, push and update the roadmap. A stage interrupted by a usage limit resumes
with the same command, from the checklist under "In progress". If the user says "continue",
"next" or "go", follow `.claude/skills/next/SKILL.md`.

## Upstream references

- Shortcut-bias paper (central reference): Zhang et al., *Cross-task
  interpretability through unified modeling reveals a universal shortcut bias in
  neoantigen prediction*, Cell Genomics, April 2026.
- Code and data: `github.com/WangLabTHU/ImmuBPI` (repo states MIT; the Zenodo
  record 18610279 lists Apache-2.0 — **check which applies before vendoring**).
  Data also on Google Drive per the repo README.
- TESLA consortium: Wells et al., Cell 2020.
- NetMHCpan-4.1 pseudo-sequence convention: 34 polymorphic residues within 4.0 Å
  of the peptide.

## Commands

```
nbx doctor       # verify the integrity core still behaves
nbx status       # which modules are built
python scripts/verify.py          # full verification loop
python scripts/verify.py --fast   # lint + core tests only
```
