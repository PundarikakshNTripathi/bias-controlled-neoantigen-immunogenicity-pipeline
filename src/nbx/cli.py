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
    sources: str = typer.Option(..., "--sources", help="Comma-separated corpus names."),
    out: Path = typer.Option(Path("artefacts/datasets"), "--out"),
) -> None:
    """Ingest public corpora into one canonical, content-hashed dataset."""
    _unbuilt("nbx.data.harmonise")


@bias_app.command("profile")
def bias_profile(
    dataset: Path = typer.Argument(..., help="Canonical dataset to profile."),
    report: Path = typer.Option(Path("artefacts/bias/bias.json"), "--report"),
) -> None:
    """Compute I(b;y), per-allele rho and the derived masking policy.

    Runs BEFORE any model is trained. Its output conditions the modelling and
    evaluation layers rather than merely annotating their results.
    """
    _unbuilt("nbx.data.harmonise (required first)")


@split_app.command("freeze")
def split_freeze(
    dataset: Path = typer.Argument(...),
    max_abs_rho: float = typer.Option(1.0, "--max-abs-rho"),
    min_samples: int = typer.Option(20, "--min-samples-per-allele"),
    seed: int = typer.Option(20260928, "--seed"),
) -> None:
    """Build and freeze an allele-balanced split (requirement R1)."""
    _unbuilt("nbx.data.harmonise (required first)")


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
