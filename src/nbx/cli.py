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
    name: str = typer.Option("iedb_tcell", "--name", help="Artefact base name."),
) -> None:
    """Ingest public corpora into one canonical, content-hashed dataset.

    This pass ingests the IEDB T cell export only. Every number printed below
    is read off the frame that was just built; none is a constant.
    """
    from nbx.data.harmonise import HarmonisationError, harmonise

    try:
        corpus, paths = harmonise(
            {"iedb_tcell": iedb_tcell, "mhc_pseudo": mhc_pseudo},
            out_dir=out,
            name=name,
        )
    except HarmonisationError as exc:
        console.print(f"[red]harmonisation refused:[/red] {exc}")
        raise typer.Exit(code=1) from exc

    report = corpus.report
    frame = corpus.frame

    summary = Table(title="canonical corpus", show_header=True, header_style="bold")
    summary.add_column("Quantity")
    summary.add_column("Value", justify="right")
    summary.add_row("records", f"{len(frame):,}")
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
        "duplicates merged", f"{report['duplicate_records_merged']:,}"
    )
    summary.add_row(
        "label-conflict pairs dropped",
        f"{report['label_conflict_pairs_dropped']:,}",
    )
    console.print(summary)

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
    _unbuilt("nbx.bias.profile (diagnostics built; CLI wiring pending)")


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
) -> None:
    """Evaluate under the bias-controlled protocol.

    Verifies R1 to R4, then reports allele-stratified metrics, shortcut
    dependence and the margin over the HLA-only control for every variant.
    """
    _unbuilt("nbx.evaluate.protocol")


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
            kappa=0.0,
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
