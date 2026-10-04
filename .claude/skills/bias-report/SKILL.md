---
name: bias-report
description: Compute and present a bias profile for a corpus, in the project's reporting format. Use after harmonisation and before any training run.
argument-hint: [dataset-path]
allowed-tools: Read Write Bash(python:*) Bash(nbx:*) Grep Glob
---

## Bias profile for $ARGUMENTS

Produce the profile that conditions everything downstream. This runs **before**
any model is trained; its output is an input to training, not a commentary on it.

Report, in this order:

1. **Corpus shape** — records, alleles, per-source breakdown, peptide length
   distribution. Real counts from the frame.
2. **I(b;y) per candidate bias variable** — at minimum `hla_allele`,
   `source_dataset`, `peptide_length`, and position-2 residue. Report each as a
   value and as a share of H(y). Name which are bias *candidates*; do not call
   anything a shortcut on the strength of mutual information alone, because a
   genuine determinant also carries information.
3. **Per-allele ρ(a)** — the distribution, the count beyond the 10:1 threshold
   (|ρ| > log2(10) ≈ 3.32), and the count of alleles excluded for having a
   single sample.
4. **The derived masking policy M(a)** — how many alleles get a mask differing
   from the uniform P2/PΩ scheme, which is the quantity that makes adaptive
   masking worth testing at all.
5. **What this implies for the splits** — how many alleles can meet the balance
   criterion, and how many must be reported separately.

Use `nbx.bias.diagnostics`. Do not reimplement the mathematics; it is tested
against analytic cases and cross-checked against scipy and sklearn.

Write the profile to `artefacts/bias/` as JSON and state the path.

Every number comes from the computation you just ran. If something could not be
computed, render it as an em dash and say why.
