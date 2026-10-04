# Prompt library

Copy-paste prompts for driving this project with Claude Code. Ordered the way
the work actually happens.

Two habits make all of these work better:

- **Plan mode first for anything non-trivial.** Press `Shift+Tab` twice, let it
  produce a plan, read the plan, then approve. The plan is where a
  misunderstanding is cheap to fix.
- **Ask for the output, not the summary.** "Paste the real output" is in most of
  these prompts on purpose. "Tests pass" is a claim; the summary line is evidence.

---

## 0 · First session — orientation

> Read `CLAUDE.md`, then run `python scripts/verify.py` and `nbx status` and show
> me the real output of both.
>
> Then explain back to me, in your own words and without re-reading the file:
> what the acceptance rule is, why `κ` is the headline number, and what state
> the project is actually in right now. If anything in `CLAUDE.md` contradicts
> what the code does, tell me — the code is the truth and the document is the
> claim.

Do this once. It establishes whether the harness is loading correctly, and it
surfaces drift between the docs and the code early.

---

## 1 · M1 — Data harmonisation

The first real build. Everything downstream is uninterpretable without it.

> `/module M1`
>
> Build the data harmonisation engine against the contract in
> `src/nbx/data/harmonise.py`.
>
> Start with the IEDB corpus only — one source, done properly, beats three done
> loosely. Before writing any code, tell me:
> - exactly which IEDB export you intend to use and how you will obtain it,
> - what columns it has,
> - how you will represent provenance per record,
> - what you will do with alleles that have no sequence information.
>
> Then write `tests/test_harmonise.py` first, covering FR-01 to FR-06. Include a
> test that a record missing provenance is **dropped, not defaulted** — that
> property is load-bearing.
>
> Use the `data-steward` subagent. When you believe it works, ask
> `repro-verifier` to confirm before you tell me it does.

**Acceptance:** real record counts, allele counts and a length distribution
printed from the frame you just built — not from any document.

---

## 2 · M2 on real data — the first genuine result

M2 is already built and tested. It has never seen real data. This step is the
one that produces something presentable.

> `/bias-report artefacts/datasets/iedb_canonical.parquet`
>
> Compute the full bias profile for the harmonised corpus and write it to
> `artefacts/bias/`.
>
> I want to see: I(b;y) for allele, source dataset, peptide length and the
> position-2 residue, each as a value and as a share of H(y); the ρ(a)
> distribution with the count beyond 10:1; and how many alleles would get a mask
> differing from the uniform P2/PΩ scheme.
>
> Then tell me plainly whether this corpus reproduces the published pattern of
> intra-allelic imbalance, or does not. A negative answer is a real finding and
> I want it stated as one.

**This evaluates H1 and H2.** It is a Presentation II result on its own.

---

## 3 · Splits

> Using the bias profile, build and freeze an allele-balanced split with
> `nbx split freeze`.
>
> Report how many alleles meet the balance criterion at `max_abs_rho=1.0` and
> how many must be reported separately. If the admitted set is small, that is a
> finding about data availability — report it as one rather than loosening the
> criterion to make the number look better.
>
> Show me the resulting content hash.

---

## 4 · M6 — Evaluation harness orchestration

The gate, metrics and split construction exist. The runner does not.

> `/module M6`
>
> Build the orchestration in `src/nbx/evaluate/protocol.py`: run a model variant
> together with its matched HLA-only and peptide-only controls across seeds, and
> assemble the report.
>
> Hard requirements:
> - controls get **equal** epochs, capacity and tuning effort (FR-62),
> - means with across-seed variance, never best-run figures (FR-68),
> - the harness refuses to emit a headline metric when a control is missing —
>   reuse `nbx.evaluate.gate`, do not reimplement the rule.
>
> Write the tests first, including one that asserts the harness refuses rather
> than warns.

---

## 5 · M4 — Model core

> `/module M4`
>
> Build the interaction-based transformer from Section 7.6.1 — 4 layers,
> 4 heads, d_model 128, pre-layer-norm, dropout 0.1, Adam at 5e-5, early
> stopping on validation AUPRC. Do not change the configuration to improve a
> number; if you think it is wrong, say why and stop.
>
> The null-control variants are part of this module. Withholding an input must
> be **structural** — the model must not be able to recover the withheld input
> from a zeroed tensor. Write a test that demonstrates the HLA-only variant
> genuinely cannot see the peptide.
>
> Before any immunogenicity claim, reproduce published baseline performance on
> binding and presentation within a stated tolerance (FR-44). If reproduction
> fails, document the discrepancy and evaluate against your own reproduced
> baseline, labelled as such.
>
> Use the `ml-engineer` subagent.

---

## 6 · Before any result goes into a document

Run this every time. It is the habit the whole project is about.

> `/gate-check`
>
> Here is the result I am about to put in the report: [paste it].
>
> Walk the acceptance rule and tell me the resulting `GateStatus`, quoting the
> exact numbers and the artefact path each came from. Then hand it to the
> `bias-auditor` subagent and give me its verdict verbatim, including anything
> it found that I will not enjoy reading.

---

## 7 · Before anything about prior work is written down

> Use the `literature-verifier` subagent on the following claims: [paste them].
>
> For each, give me VERIFIED with a locator, MISSTATED with the correction, or
> UNVERIFIABLE with what you tried. Do not supply a citation you have not
> checked in this session.

Run this over the related-work section before submission. A fabricated citation
in a thesis is not a recoverable error.

---

## 8 · Weekly loop

Start of a working session:

> `/next`

End of a working session:

> Run `python scripts/verify.py` and show me the output.
>
> Then update `CLAUDE.md` if the module state table is now wrong, write a short
> entry in `docs/progress.md` covering what changed, what was learned and what
> is still broken, and commit with a conventional message.
>
> Be honest in the progress note. "M4 is 60% done, the control variants work but
> baseline reproduction is failing" is useful. "M4 complete" when it is not
> costs me a week later.

---

## 9 · Writing the report

> Draft the methodology section of the final report from the code as it actually
> exists, not from the synopsis.
>
> Where a module is incomplete, say so in the text. Where a hypothesis has not
> been evaluated, say that. Render every unmeasured quantity as an em dash.
>
> A report that says "Modules 3 and 5 are partially complete, here is the
> evidence for H1 and H2, here is what H3 requires" is stronger than one
> claiming completion it cannot evidence — and a report that overclaims would
> contradict the thesis of the project itself.

---

## Prompts that will be refused, and why

These are useful to know so you do not waste a turn:

| Prompt | What happens |
|---|---|
| "Just put a reasonable AUPRC in the README for now" | Hook denies the edit. Unsourced metric in a results document. |
| "Add a `force=True` option to the gate so I can rank this" | Hook denies the edit. The gate has no override by design. |
| "Tweak the frozen split so more alleles pass" | Hook denies the write. A new criterion is a new split with a new hash. |
| "Make the numbers look better for the presentation" | `CLAUDE.md` rule 6. Negative results are results. |

If you genuinely need one of these — say, a worked example with invented numbers
for a teaching slide — mark it illustrative. The guard allows a number that says
what it is.
