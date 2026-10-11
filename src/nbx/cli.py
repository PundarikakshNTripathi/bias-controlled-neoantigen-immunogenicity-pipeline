"""The nbx command line: the authoritative surface.

Every other interface in this project is a read-mostly projection of the
artefacts these commands write. A run is fully specified by a declarative
config plus a split hash; no interactive state determines a result.

Two behaviours are load-bearing and are enforced here rather than left to
discipline:

* `nbx eval` refuses to print a headline metric when a mandated null control
  is missing (FR-63). The refusal is a non-zero exit, not a warning.
* Any quantity that has not been measured renders as an em dash, never as
  zero and never as a plausible-looking stand-in (NFR-15).
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from nbx.evaluate.gate import UNMEASURED, GateStatus, render_metric

app = typer.Typer(
    name="nbx",
    help="A bias-controlled benchmark protocol for neoantigen immunogenicity "
         "prediction.",
    no_args_is_help=True,
    add_completion=False,
)
console = Console()

data_app = typer.Typer(
    help="Ingest and harmonise public corpora.", no_args_is_help=True
)
bias_app = typer.Typer(
    help="Quantify dataset bias before training.", no_args_is_help=True
)
split_app = typer.Typer(
    help="Construct and freeze allele-balanced splits.", no_args_is_help=True
)
repr_app = typer.Typer(
    help="Compute and cache fixed representations.", no_args_is_help=True
)
app.add_typer(data_app, name="data")
app.add_typer(bias_app, name="bias")
app.add_typer(split_app, name="split")
app.add_typer(repr_app, name="repr")

NOT_BUILT = (
    "[yellow]Not implemented yet.[/yellow] This command is a declared part of "
    "the interface specification (UI/UX Document, Table 5.1) and has a "
    "contract in the source. See PROMPTS.md for the prompt that builds it."
)


def _unbuilt(module: str) -> None:
    console.print(NOT_BUILT)
    console.print(f"[dim]module: {module}[/dim]")
    raise typer.Exit(code=2)


# --------------------------------------------------------------------------
@app.command()
def init(
    path: Path = typer.Argument(Path("."), help="Project directory to initialise."),
) -> None:
    """Scaffold a project: config directory, artefact store, environment lock."""
    store = path / "artefacts"
    for sub in ("datasets", "bias", "splits", "cache", "runs", "reports"):
        (store / sub).mkdir(parents=True, exist_ok=True)
    console.print(f"[green]Artefact store ready at[/green] {store}")


@data_app.command("harmonise")
def data_harmonise(
    iedb_tcell: Path = typer.Option(
        ...,
        "--iedb-tcell",
        help="IEDB T cell assay export (tcell_full_v3.zip or .csv).",
    ),
    mhc_pseudo: Path = typer.Option(
        ...,
        "--mhc-pseudo",
        help="Allele -> 34-residue pseudo-sequence table (NAME SEQUENCE per line).",
    ),
    out: Path = typer.Option(Path("artefacts/datasets"), "--out"),
    label_rule: str = typer.Option(
        None,
        "--label-rule",
        help="Pair label rule: majority (primary), strict or any_positive. "
             "Part of the hashed criterion, so each rule is its own artefact.",
    ),
    name: str = typer.Option(
        None,
        "--name",
        help="Artefact base name. Defaults to iedb_tcell_<label-rule>.",
    ),
) -> None:
    """Ingest public corpora into one canonical, content-hashed dataset.

    This pass ingests the IEDB T cell export only. One output row is one
    (peptide, allele) pair. Every number printed below is read off the frame
    that was just built; none is a constant.
    """
    from nbx.data.harmonise import (
        DEFAULT_LABEL_RULE,
        LABEL_RULES,
        HarmonisationError,
        harmonise,
    )

    rule = label_rule or DEFAULT_LABEL_RULE
    if rule not in LABEL_RULES:
        console.print(
            f"[red]unknown label rule:[/red] {rule!r}. "
            f"Choose one of {list(LABEL_RULES)}."
        )
        raise typer.Exit(code=1)

    try:
        corpus, paths = harmonise(
            {"iedb_tcell": iedb_tcell, "mhc_pseudo": mhc_pseudo},
            out_dir=out,
            name=name,
            label_rule=rule,
        )
    except HarmonisationError as exc:
        console.print(f"[red]harmonisation refused:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    report = corpus.report
    frame = corpus.frame

    summary = Table(title="canonical corpus", show_header=True, header_style="bold")
    summary.add_column("Quantity")
    summary.add_column("Value", justify="right")
    summary.add_row("label rule", report["label_rule"])
    summary.add_row("unit of analysis", report["unit_of_analysis"])
    summary.add_row("pairs (records)", f"{len(frame):,}")
    summary.add_row("peptides", f"{frame['peptide'].nunique():,}")
    summary.add_row("alleles", f"{frame['allele'].nunique():,}")
    summary.add_row("content hash", corpus.short_hash)
    for label, count in sorted(report["label_counts"].items()):
        summary.add_row(f"label {label}", f"{count:,}")
    summary.add_row("dropped (total)", f"{report['dropped_total']:,}")
    summary.add_row(
        "dropped: missing provenance", f"{report['dropped_for_missing_provenance']:,}"
    )
    summary.add_row(
        "excluded for length", f"{report['peptides_excluded_for_length']:,}"
    )
    summary.add_row(
        "assay records collapsed into pairs",
        f"{report['duplicate_records_merged']:,}",
    )
    summary.add_row(
        "pairs dropped by the label rule",
        f"{report['pairs_dropped_by_label_rule']:,}",
    )
    console.print(summary)

    patterns = report["pair_pattern_summary"]
    pattern_table = Table(
        title="pair disagreement pattern (every pair that reached the rule)",
        show_header=True,
        header_style="bold",
    )
    pattern_table.add_column("Pattern")
    pattern_table.add_column("Pairs", justify="right")
    for key in ("unanimous_positive", "unanimous_negative", "mixed", "tied"):
        pattern_table.add_row(key, f"{patterns[key]:,}")
    pattern_table.add_row("total", f"{patterns['pairs_total']:,}")
    console.print(pattern_table)
    console.print(f"[dim]{patterns['note']}[/dim]")

    rec = report["reconciliation"]
    recon = Table(title="reconciliation", show_header=True, header_style="bold")
    recon.add_column("Quantity")
    recon.add_column("Value", justify="right")
    for key in (
        "raw_row_count",
        "records_dropped_record_level",
        "records_dropped_pair_level",
        "records_collapsed_into_pairs",
        "pairs_labelled",
        "pairs_dropped",
        "pairs_total",
    ):
        recon.add_row(key, f"{rec[key]:,}")
    recon.add_row("record level balances", str(rec["record_level_balances"]))
    recon.add_row("pair level balances", str(rec["pair_level_balances"]))
    console.print(recon)

    lengths = Table(title="peptide length", show_header=True, header_style="bold")
    lengths.add_column("Length")
    lengths.add_column("Records", justify="right")
    for length in sorted(report["length_distribution"], key=int):
        lengths.add_row(length, f"{report['length_distribution'][length]:,}")
    console.print(lengths)

    sources_table = Table(title="per source", show_header=True, header_style="bold")
    sources_table.add_column("Source")
    sources_table.add_column("Records", justify="right")
    sources_table.add_column("Alleles", justify="right")
    sources_table.add_column("Positive", justify="right")
    sources_table.add_column("Negative", justify="right")
    for source, stats in sorted(report["per_source_breakdown"].items()):
        sources_table.add_row(
            source,
            f"{stats['record_count']:,}",
            f"{stats['allele_count']:,}",
            f"{stats['n_positive']:,}",
            f"{stats['n_negative']:,}",
        )
    console.print(sources_table)

    drops = Table(title="dropped by reason", show_header=True, header_style="bold")
    drops.add_column("Reason")
    drops.add_column("Records", justify="right")
    for reason, count in sorted(
        report["dropped_by_reason"].items(), key=lambda kv: (-kv[1], kv[0])
    ):
        drops.add_row(reason, f"{count:,}")
    console.print(drops)

    absent = report["alleles_absent_from_pseudo_sequence_table"]
    gap = Table(
        title=(
            f"alleles resolved but absent from the pseudo-sequence table "
            f"({absent['n_alleles']:,} alleles, {absent['n_records']:,} records; "
            f"top {absent['n_listed']:,} by records)"
        ),
        show_header=True,
        header_style="bold",
    )
    gap.add_column("Allele")
    gap.add_column("Records dropped", justify="right")
    # An empty table here is a measured empty set (the counts in the title say
    # so), not an unmeasured quantity, so it is left empty rather than filled
    # with em dashes.
    for entry in absent["alleles"]:
        gap.add_row(entry["allele"], f"{entry['records']:,}")
    console.print(gap)
    console.print(f"[dim]{absent['note']}[/dim]")

    comparison = report["reference_comparison"]
    refs = Table(
        title="published-statistic comparison (FR-06)",
        show_header=True,
        header_style="bold",
    )
    refs.add_column("Statistic")
    refs.add_column("Observed", justify="right")
    refs.add_column("Published", justify="right")
    refs.add_column("Discrepancy", justify="right")
    refs.add_column("Status")
    for entry in comparison["statistics"]:
        refs.add_row(
            entry["name"],
            entry["observed"],
            entry["reference"],
            entry["discrepancy"],
            entry["status"],
        )
    console.print(refs)
    if comparison["note"]:
        console.print(f"[dim]{comparison['note']}[/dim]")

    console.print(f"[green]dataset[/green] {paths['dataset']}")
    console.print(f"[green]report [/green] {paths['report']}")


@bias_app.command("profile")
def bias_profile(
    dataset: Path = typer.Argument(..., help="Canonical dataset to profile."),
    report: Path = typer.Option(Path("artefacts/bias/bias.json"), "--report"),
) -> None:
    """Compute I(b;y), per-allele rho and the derived masking policy.

    Runs BEFORE any model is trained. Its output conditions the modelling and
    evaluation layers rather than merely annotating their results.
    """
    import hashlib

    import pandas as pd

    from nbx.bias.profile import BiasReportError, build_bias_report, write_bias_report

    if not dataset.exists():
        console.print(f"[red]dataset not found:[/red] {dataset}")
        raise typer.Exit(code=1)
    try:
        frame = pd.read_parquet(dataset)
        result = build_bias_report(
            frame,
            input_path=str(dataset),
            input_sha256=hashlib.sha256(dataset.read_bytes()).hexdigest(),
        )
        written = write_bias_report(result, report)
    except BiasReportError as exc:
        console.print(f"[red]bias profile refused:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    mi = Table(title="I(b;y): bias candidates, not verdicts", header_style="bold")
    mi.add_column("Variable")
    mi.add_column("I (nats)", justify="right")
    mi.add_column("Share of H(y)", justify="right")
    mi.add_column("Note")
    for name, row in result["mutual_information"].items():
        note = ""
        if row["structurally_uninformative"]:
            note = "single-valued: zero by construction"
        mi.add_row(
            name,
            render_metric(row["nats"]),
            render_metric(row["share_of_label_entropy"]),
            note,
        )
    console.print(mi)

    rho = result["per_allele_rho"]
    pol = result["masking_policy"]
    split = result["split_admissibility"]
    summary = Table(title="corpus", header_style="bold")
    summary.add_column("Quantity")
    summary.add_column("Value", justify="right")
    summary.add_row("records", f"{result['n_records']:,}")
    summary.add_row("alleles", f"{result['n_alleles']:,}")
    summary.add_row("alleles with rho", f"{rho['n_alleles_with_rho']:,}")
    summary.add_row("alleles beyond 10:1", f"{rho['n_beyond_10_to_1']:,}")
    summary.add_row(
        "alleles excluded (single sample)", f"{rho['n_excluded_single_sample']:,}"
    )
    summary.add_row(
        "policy scope: 9-mers, alleles with enough records",
        f"{pol['n_alleles_in_scope']:,} of {pol['n_alleles_total']:,} alleles, "
        f"{pol['n_records_in_scope']:,} of {pol['n_records_total']:,} records",
    )
    summary.add_row(
        "policy: alleles with a non-empty mask", str(pol["n_alleles_nonempty_mask"])
    )
    summary.add_row(
        "split: admitted / separate",
        f"{split['n_admitted_alleles']:,} / {split['n_reported_separately']:,}",
    )
    console.print(summary)
    console.print(f"[green]report[/green] {written}")


@split_app.command("freeze")
def split_freeze(
    dataset: Path = typer.Argument(...),
    max_abs_rho: float = typer.Option(1.0, "--max-abs-rho"),
    min_samples: int = typer.Option(20, "--min-samples-per-allele"),
    seed: int = typer.Option(20260928, "--seed"),
) -> None:
    """Build and freeze an allele-balanced split (requirement R1)."""
    _unbuilt("nbx.splits.freeze (built; CLI wiring pending)")


@repr_app.command("cache")
def repr_cache(
    dataset: Path = typer.Argument(...),
    with_geometry: bool = typer.Option(False, "--with-geometry"),
) -> None:
    """Compute and cache language-model and structural representations."""
    _unbuilt("nbx.represent.cache")


@app.command()
def train(
    matrix: Path | None = typer.Option(None, "--matrix", help="Ablation config."),
    seeds: int = typer.Option(5, "--seeds"),
) -> None:
    """Run the ablation matrix, with the mandated null controls queued.

    The null controls are not optional and cannot be disabled from the
    command line. They are queued automatically alongside every matrix.
    """
    _unbuilt("nbx.model.backbone")


@app.command("eval")
def evaluate(
    run: Path = typer.Option(Path("artefacts/runs"), "--run"),
    split: str = typer.Option(..., "--split", help="Frozen split content hash."),
    protocol: str = typer.Option("bias-controlled", "--protocol"),
    out: Path | None = typer.Option(
        None, "--out", help="Also write the content-hashed JSON report here."
    ),
    records_digest: str | None = typer.Option(
        None,
        "--records-digest",
        help="Pin the records digest (sha256 of per-record allele, label, peptide). "
             "Without it the digest only detects an inconsistent split.json.",
    ),
) -> None:
    """Evaluate under the bias-controlled protocol.

    Verifies R1 to R4, then reports allele-stratified metrics, shortcut
    dependence and the margin over the HLA-only control for every variant.
    The bootstrap is always the pre-registered one (B and seed are not options).
    A refusal is a non-zero exit and prints no metric.
    """
    from nbx.evaluate import protocol as proto

    if protocol != "bias-controlled":
        console.print(f"unknown protocol: {protocol!r}", style="red", markup=False)
        raise typer.Exit(code=2)
    try:
        result = proto.run_protocol(
            run,
            split,
            expected_records_digest=records_digest,
            n_boot=proto.PRODUCTION_N_BOOT,
            bootstrap_seed=proto.DEFAULT_BOOTSTRAP_SEED,
        )
    except proto.ProtocolRefusal as exc:
        console.print(f"protocol refused: {exc}", style="red", markup=False)
        raise typer.Exit(code=1) from exc

    def note(text: str) -> None:
        console.print(text, style="dim", markup=False)

    table = Table(
        title="bias-controlled evaluation (mean, across-seed SD)",
        show_header=True,
        header_style="bold",
    )
    for col in ("Variant", "Role", "Gate", "Rank", "AUPRC", "Margin vs HLA-only",
                "kappa"):
        table.add_column(col)
    for name, v in result["variants"].items():
        shown = v["metrics"]["displayed"]
        auprc_cell = UNMEASURED
        margin_cell = v.get("margin_vs_hla_only", {}).get("rendered", "control row")
        kappa_cell = UNMEASURED
        if shown:
            m = v["metrics"]["pooled"]["auprc"]
            auprc_cell = f"{m['mean']} +/- {m['sd']}"
            kappa_cell = v["shortcut_dependence"]["kappa_mean"]
        table.add_row(name, v["role"], v["gate"]["status"], str(v["rank"]),
                      auprc_cell, margin_cell, kappa_cell)
    console.print(table)

    ts = result["test_set"]
    note(f"{ts['n_alleles']} test alleles, {ts['n_records']} test records")
    for name, v in result["variants"].items():
        if v["role"] not in ("model", "baseline"):
            continue
        note(f"{name}: {v['gate']['reason']}")
        if not v["metrics"]["displayed"]:
            note(f"{name}: metrics withheld (gate status pending)")
            continue
        sd = v["shortcut_dependence"]
        if v.get("reference_only"):
            note(f"{name}: reference_only; baseline kappa {sd['kappa_mean']}; "
                 "never ranked")
            continue
        h, p = v["margin_vs_hla_only"], v["margin_vs_peptide_only"]
        change = sd["abs_kappa_change"]
        note(f"{name}: kappa {sd['kappa_mean']}; baseline kappa {sd['baseline_kappa']}; "
             f"kappa change upper95 {change['upper95_one_sided']} "
             f"(non-inferiority margin {change['noninferiority_margin']})")
        note(f"{name}: HLA-only margin p {h['p_value']}, Holm-adjusted "
             f"{h.get('p_value_holm', UNMEASURED)}, 90% interval "
             f"[{h['ci90_low']}, {h['ci90_high']}]")
        note(f"{name}: peptide-only margin {p['rendered']}, p {p['p_value']}, "
             f"90% interval [{p['ci90_low']}, {p['ci90_high']}], "
             f"significantly worse: {p['significantly_worse']}")
    match = result["matching"]
    ver = match["structurally_verified"]
    note("matching: declared in config.json and compared, not verified: "
         + ", ".join(match["declared_metadata"]))
    note(f"matching: structurally verified: HLA-only constant within every allele "
         f"({ver['hla_only']['n_runs_checked']} runs); peptide-only identical on "
         f"shared peptides ({ver['peptide_only']['n_shared_peptide_groups']} shared "
         f"groups, {ver['peptide_only']['n_runs_checked']} runs)")
    rec = result["records"]
    note(f"split {result['split_hash'][:12]} (definition); records digest "
         f"{rec['sha256'][:12]} verified against split.json (allele, label, peptide); "
         f"pinned by caller: {'yes' if rec['pinned_by_caller'] else 'no'}")
    note(f"{result['n_model_variants_compared']} model variant(s) compared; Holm-"
         "adjusted p reported beside the unadjusted p")
    note(f"seeds {result['seeds']}, B={result['bootstrap']['B']}, "
         f"report {result['report_hash'][:12]}")
    if not result["preregistered_settings_used"]:
        note("settings NOT pre-registered in this report")
    for item in result["not_done"]:
        note(f"not done: {item}")
    note(result["notice"])
    if out is not None:
        console.print(f"report {proto.write_report(result, out)}", style="green",
                      markup=False)


@app.command()
def report(
    out: Path = typer.Option(Path("artefacts/reports"), "--out"),
) -> None:
    """Package the benchmark report: deliverable D1."""
    _unbuilt("nbx.release.package")


@app.command()
def doctor() -> None:
    """Check that the integrity core is present and behaving."""
    ok = True
    table = Table(title="nbx doctor", show_header=True, header_style="bold")
    table.add_column("Check")
    table.add_column("Result")

    try:
        import math

        from nbx.bias.diagnostics import entropy, mutual_information

        assert abs(entropy([0, 1]) - math.log(2)) < 1e-12
        assert abs(mutual_information([0, 0, 1, 1], ["a", "a", "b", "b"])
                   - math.log(2)) < 1e-12
        table.add_row("bias diagnostics", "[green]ok[/green]")
    except Exception as exc:  # pragma: no cover - diagnostic path
        table.add_row("bias diagnostics", f"[red]{exc}[/red]")
        ok = False

    try:
        from nbx.evaluate.gate import ControlMargin, evaluate_gate

        pending = evaluate_gate(
            headline_metric=0.9,
            hla_only_margin=ControlMargin(),
            peptide_only_margin=ControlMargin(delta=0.1, p_value=0.01),
            kappa_change_upper_bound=0.0,
            peptide_only_significantly_worse=False,
        )
        assert pending.status is GateStatus.PENDING
        assert not pending.displayable_metrics
        table.add_row("acceptance gate", "[green]ok[/green]")
    except Exception as exc:  # pragma: no cover
        table.add_row("acceptance gate", f"[red]{exc}[/red]")
        ok = False

    try:
        assert render_metric(None) == UNMEASURED
        table.add_row("unmeasured rendering", f"[green]ok ({UNMEASURED})[/green]")
    except Exception as exc:  # pragma: no cover
        table.add_row("unmeasured rendering", f"[red]{exc}[/red]")
        ok = False

    console.print(table)
    if not ok:
        raise typer.Exit(code=1)


@app.command()
def status() -> None:
    """Show which modules are built and which are contracts only."""
    modules = [
        ("M1", "Data harmonisation", "nbx.data.harmonise"),
        ("M2", "Bias diagnostics", "nbx.bias.diagnostics"),
        ("M3", "Representation service", "nbx.represent.cache"),
        ("M4", "Model core", "nbx.model.backbone"),
        ("M5", "Debiasing controller", "nbx.bias.policy"),
        ("M6", "Evaluation harness", "nbx.evaluate.protocol"),
        ("M7", "Release packaging", "nbx.release.package"),
    ]
    table = Table(show_header=True, header_style="bold")
    table.add_column("Module")
    table.add_column("Name")
    table.add_column("State")

    import importlib

    for tag, name, mod in modules:
        try:
            m = importlib.import_module(mod)
            src = Path(m.__file__ or "").read_text(encoding="utf-8")
            built = "NotImplementedError" not in src
        except Exception:
            built = False
        table.add_row(
            tag, name,
            "[green]built[/green]" if built else f"[dim]contract only {UNMEASURED}[/dim]"
        )
    console.print(table)


if __name__ == "__main__":  # pragma: no cover
    app()
