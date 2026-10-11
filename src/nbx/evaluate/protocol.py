"""Module 6: Bias-Controlled Evaluation Harness (FR-60 to FR-69).

The harness is model-agnostic. Training lives in M4; this module consumes the
PREDICTIONS of runs and produces the evaluation: metrics with across-seed
variance, the margin over the matched HLA-only control, shortcut dependence,
and the acceptance verdict of `nbx.evaluate.gate`. It does not reimplement the
acceptance rule or the metrics; it feeds them.

Status
------
Built: input contract and on-disk layout, refusal conditions (FR-62, FR-63),
structural validation of the controls, pooled and allele-stratified metrics
(FR-68), paired allele-clustered bootstrap for the margin over the HLA-only
control (analysis plan section 5), the kappa-change bound that the gate reads,
the content-hashed report, and the `nbx eval` command.

NOT done: FR-65 (distribution-shifted control sets, in which the intra-allelic
imbalance of the test data is inverted relative to training), the module's
acceptance reproduction against a published model on a published benchmark, and
pinning the baseline's identity in the analysis plan (it is still read from the
run directory).

On-disk layout (what `run_protocol` reads and `write_run_dir` writes)
---------------------------------------------------------------------
    <run_dir>/split.json                 the frozen split and the corpus records
    <run_dir>/<variant>/config.json      one per variant
    <run_dir>/<variant>/seed_<n>.parquet columns: idx (int), score (float)

`split.json` carries the split's full definition so its content hash can be
recomputed, plus the per-record allele, label and peptide of the corpus (row i of
the corpus is index i) and a digest of those records:

    {"content_hash": str, "criterion": {...SplitCriterion canonical fields...},
     "n_records": int, "train_idx": [...], "val_idx": [...], "test_idx": [...],
     "quarantined_idx": [...], "admitted_alleles": [...],
     "reported_separately": [...],
     "records": {"allele": [...], "label": [...], "peptide": [...]},
     "records_sha256": sha256 of the canonical JSON of "records"}

The FrozenSplit hash does not cover labels, alleles or peptides. The records
digest does, and the harness recomputes it. A digest that lives next to the data
it vouches for detects inconsistency, not a coordinated edit of both; a caller who
wants the second guarantee pins the digest (`expected_records_digest`, or
`--records-digest` on the command line), and the report says whether it did.

`config.json`:

    {"variant": str (must equal the directory name),
     "role": "model" | "baseline" | "hla_only" | "peptide_only",
     "matched_to": null for model/baseline, else the variant it is a control of,
     "run_config": {"epochs": int, "n_parameters": int,
                    "early_stopping_metric": str, "tuning_trials": int}}

A score file must hold exactly one score for every test record and nothing else.
The trained seeds of every variant must be exactly EXPECTED_SEEDS.

Refusal
-------
`run_protocol` raises `ProtocolRefusal` and emits nothing when the inputs cannot
support a bias-controlled comparison. A refusal is not a warning: there is no
partial report. The conditions are listed in `_validate`.

Decisions this module makes that the analysis plan does not specify
-------------------------------------------------------------------
* Ties are handled by a seed-free definition (see `nbx.evaluate.stats`): AUPRC is
  the expected average precision over random orderings within ties, and the ranking
  metrics are expectations the same way. `bootstrap_seed` seeds the resampling only.
* A model with no baseline run cannot be shown not to have raised shortcut
  dependence, so the kappa interval is withheld from the gate and the variant is
  PENDING. Every model also needs its own matched controls; so does the baseline.
* The baseline row is evaluated like a model (it needs controls, it gets a margin)
  but its kappa is compared with itself, so the shortcut condition is vacuous for
  it; it is flagged `reference_only` and never ranked.
* MCC needs hard predictions. The rule is: predict positive for the n_pos highest
  scores, where n_pos is the number of positives in the evaluated set.
* The parameter-count tolerance for matched capacity defaults to exact and is
  written into the report whenever it is used.
* The structural HLA-only check demands EXACT equality within an allele, so a
  trainer must emit identical floats for identical inputs.
* The gate reads the unadjusted one-sided p. The Holm-adjusted p across the model
  variants of a run is reported beside it and is not used by the gate.
"""

from __future__ import annotations

import functools
import hashlib
import json
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from nbx.bias.diagnostics import per_allele_imbalance, shortcut_dependence
from nbx.evaluate import stats
from nbx.evaluate.gate import (
    NONINFERIORITY_MARGIN_KAPPA,
    UNMEASURED,
    ControlMargin,
    GateResult,
    evaluate_gate,
    render_metric,
)
from nbx.evaluate.metrics import allele_stratified, auroc
from nbx.evaluate.stats import ALPHA
from nbx.splits.freeze import FrozenSplit, SplitCriterion

__all__ = [
    "MIN_SEEDS",
    "EXPECTED_SEEDS",
    "ALPHA",
    "PRODUCTION_N_BOOT",
    "DEFAULT_BOOTSTRAP_SEED",
    "DEFAULT_RANKING_K",
    "ROLES",
    "ProtocolRefusal",
    "RunConfig",
    "VariantRun",
    "SplitRecord",
    "write_run_dir",
    "load_run_dir",
    "evaluate_runs",
    "run_protocol",
    "report_to_json",
    "write_report",
]

EXPECTED_SEEDS: tuple[int, ...] = (1, 2, 3, 4, 5)
"""Pre-registered (analysis plan section 5): every trained variant uses exactly
these seeds. Training more and submitting a subset is a forking path, so it is
refused, as is any other set of five."""

MIN_SEEDS = len(EXPECTED_SEEDS)
"""The analysis plan's minimum, kept as a named constant (at least five)."""

ROLES = ("model", "baseline", "hla_only", "peptide_only")
CONTROL_ROLES = ("hla_only", "peptide_only")

PRODUCTION_N_BOOT = 10_000
"""Analysis plan section 5. The command line always uses this value."""

DEFAULT_BOOTSTRAP_SEED = 20261011
"""Arbitrary but fixed in advance; it is recorded in every report. The command
line always uses this value."""

DEFAULT_RANKING_K = 100
"""k for PPV@k and Recall@k: the TESLA top-100 convention used by fraction ranked."""

DIGITS = 6
FR63 = "No headline metric may be displayed without both null controls (FR-63)."

_NOT_DONE = (
    "FR-65: distribution-shifted control sets (intra-allelic imbalance of the test "
    "data inverted relative to training) are not evaluated by this harness.",
    "M6 acceptance reproduction against a published model on a published benchmark "
    "has not been done.",
    "Condition (iii) of the gate, an absolute ceiling on |kappa|, is not implemented: "
    "its value is the researcher's decision.",
    "The baseline's identity is read from the run directory, not pinned in the "
    "analysis plan.",
)


class ProtocolRefusal(Exception):  # noqa: N818 - the name is part of the contract
    """The inputs cannot support a bias-controlled evaluation. Nothing is emitted."""


# ==========================================================================
# Inputs
# ==========================================================================
@dataclass(frozen=True)
class RunConfig:
    """The training budget of a run. Matched controls must equal their model."""

    epochs: int
    n_parameters: int
    early_stopping_metric: str
    tuning_trials: int


@dataclass(frozen=True)
class VariantRun:
    """One variant trained with one seed, reduced to its test predictions."""

    variant: str
    role: str
    seed: int
    scores: Mapping[int, float]
    config: RunConfig
    matched_to: str | None = None
    """For a control, the model (or baseline) variant it is matched to."""


def _digest_records(
    alleles: Sequence[str], labels: Sequence[int], peptides: Sequence[str]
) -> str:
    payload = json.dumps(
        {
            "allele": [str(a) for a in alleles],
            "label": [int(y) for y in labels],
            "peptide": [str(p) for p in peptides],
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class SplitRecord:
    """A frozen split plus the corpus allele, label and peptide of every record."""

    content_hash: str
    criterion: dict[str, Any]
    n_records: int
    train_idx: tuple[int, ...]
    val_idx: tuple[int, ...]
    test_idx: tuple[int, ...]
    quarantined_idx: tuple[int, ...]
    admitted_alleles: tuple[str, ...]
    reported_separately: tuple[str, ...]
    alleles: tuple[str, ...]
    labels: tuple[int, ...]
    peptides: tuple[str, ...]
    records_sha256: str
    """The digest DECLARED alongside the records; `_validate` recomputes it."""

    @classmethod
    def from_frozen(
        cls,
        split: FrozenSplit,
        alleles: Sequence[Any],
        labels: Sequence[int],
        peptides: Sequence[Any],
    ) -> SplitRecord:
        if not (len(alleles) == len(labels) == len(peptides) == split.n_records):
            raise ValueError("alleles, labels and peptides must cover every record")
        a = tuple(str(x) for x in alleles)
        y = tuple(int(x) for x in labels)
        p = tuple(str(x) for x in peptides)
        return cls(
            content_hash=split.content_hash,
            criterion=json.loads(split.criterion.canonical_json()),
            n_records=split.n_records,
            train_idx=tuple(split.train_idx),
            val_idx=tuple(split.val_idx),
            test_idx=tuple(split.test_idx),
            quarantined_idx=tuple(split.quarantined_idx),
            admitted_alleles=tuple(split.admitted_alleles),
            reported_separately=tuple(split.reported_separately),
            alleles=a,
            labels=y,
            peptides=p,
            records_sha256=_digest_records(a, y, p),
        )

    def with_records(
        self,
        alleles: Sequence[Any] | None = None,
        labels: Sequence[int] | None = None,
        peptides: Sequence[Any] | None = None,
    ) -> SplitRecord:
        """A copy with edited records and a freshly computed, consistent digest.

        The split hash is unchanged, because it does not cover the records. Used to
        build coordinated edits and structural scenarios in tests.
        """
        a = tuple(str(x) for x in (self.alleles if alleles is None else alleles))
        y = tuple(int(x) for x in (self.labels if labels is None else labels))
        p = tuple(str(x) for x in (self.peptides if peptides is None else peptides))
        return replace(
            self, alleles=a, labels=y, peptides=p,
            records_sha256=_digest_records(a, y, p),
        )

    def recomputed_records_digest(self) -> str:
        return _digest_records(self.alleles, self.labels, self.peptides)

    def recomputed_hash(self) -> str:
        frozen = FrozenSplit(
            criterion=SplitCriterion(**self.criterion),
            train_idx=self.train_idx,
            val_idx=self.val_idx,
            test_idx=self.test_idx,
            admitted_alleles=self.admitted_alleles,
            reported_separately=self.reported_separately,
            n_records=self.n_records,
            quarantined_idx=self.quarantined_idx,
        )
        return frozen.content_hash

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "content_hash": self.content_hash,
            "criterion": self.criterion,
            "n_records": self.n_records,
            "train_idx": list(self.train_idx),
            "val_idx": list(self.val_idx),
            "test_idx": list(self.test_idx),
            "quarantined_idx": list(self.quarantined_idx),
            "admitted_alleles": list(self.admitted_alleles),
            "reported_separately": list(self.reported_separately),
            "records": {
                "allele": list(self.alleles),
                "label": list(self.labels),
                "peptide": list(self.peptides),
            },
            "records_sha256": self.records_sha256,
        }

    @classmethod
    def from_json_dict(cls, doc: Mapping[str, Any]) -> SplitRecord:
        try:
            records = doc["records"]
            return cls(
                content_hash=str(doc["content_hash"]),
                criterion=dict(doc["criterion"]),
                n_records=int(doc["n_records"]),
                train_idx=tuple(int(i) for i in doc["train_idx"]),
                val_idx=tuple(int(i) for i in doc["val_idx"]),
                test_idx=tuple(int(i) for i in doc["test_idx"]),
                quarantined_idx=tuple(int(i) for i in doc.get("quarantined_idx", [])),
                admitted_alleles=tuple(str(a) for a in doc["admitted_alleles"]),
                reported_separately=tuple(str(a) for a in doc["reported_separately"]),
                alleles=tuple(str(a) for a in records["allele"]),
                labels=tuple(int(y) for y in records["label"]),
                peptides=tuple(str(p) for p in records["peptide"]),
                records_sha256=str(doc["records_sha256"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ProtocolRefusal(f"split.json is malformed: {exc!r}") from exc


# ==========================================================================
# Disk layout
# ==========================================================================
def write_run_dir(run_dir: Path, record: SplitRecord, runs: Sequence[VariantRun]) -> Path:
    """Write the layout `load_run_dir` reads. Used by M4 and by the test suite."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "split.json").write_text(
        json.dumps(record.to_json_dict(), sort_keys=True), encoding="utf-8"
    )
    configs: dict[str, dict[str, Any]] = {}
    for run in runs:
        doc = {
            "variant": run.variant,
            "role": run.role,
            "matched_to": run.matched_to,
            "run_config": asdict(run.config),
        }
        if configs.setdefault(run.variant, doc) != doc:
            raise ValueError(f"variant {run.variant!r} has differing configs by seed")
        vdir = run_dir / run.variant
        vdir.mkdir(exist_ok=True)
        idx = sorted(run.scores)
        frame = pd.DataFrame(
            {
                "idx": np.array(idx, dtype=np.int64),
                "score": np.array([run.scores[i] for i in idx], dtype=np.float64),
            }
        )
        frame.to_parquet(vdir / f"seed_{run.seed}.parquet", index=False)
    for variant, doc in configs.items():
        (run_dir / variant / "config.json").write_text(
            json.dumps(doc, sort_keys=True), encoding="utf-8"
        )
    return run_dir


_SEED_FILE = re.compile(r"seed_(\d+)\.parquet")


def load_run_dir(run_dir: Path) -> tuple[SplitRecord, list[VariantRun]]:
    """Read the layout. Anything malformed is a refusal, never a partial load."""
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise ProtocolRefusal(f"run directory not found: {run_dir.name!r}")
    split_path = run_dir / "split.json"
    if not split_path.is_file():
        raise ProtocolRefusal("run directory has no split.json")
    try:
        record = SplitRecord.from_json_dict(
            json.loads(split_path.read_text(encoding="utf-8"))
        )
    except json.JSONDecodeError as exc:
        raise ProtocolRefusal(f"split.json is not valid JSON: {exc}") from exc

    runs: list[VariantRun] = []
    for vdir in sorted(p for p in run_dir.iterdir() if p.is_dir()):
        cfg_path = vdir / "config.json"
        if not cfg_path.is_file():
            raise ProtocolRefusal(f"variant directory {vdir.name!r} has no config.json")
        try:
            doc = json.loads(cfg_path.read_text(encoding="utf-8"))
            config = RunConfig(**doc["run_config"])
            role, matched, name = str(doc["role"]), doc.get("matched_to"), doc["variant"]
        except (json.JSONDecodeError, KeyError, TypeError) as exc:
            raise ProtocolRefusal(
                f"config.json of {vdir.name!r} is malformed: {exc!r}"
            ) from exc
        if name != vdir.name:
            raise ProtocolRefusal(
                f"config.json names variant {name!r} inside directory {vdir.name!r}"
            )
        seed_files = sorted(vdir.glob("*.parquet"))
        if not seed_files:
            raise ProtocolRefusal(f"variant {vdir.name!r} has no seed files")
        for path in seed_files:
            m = _SEED_FILE.fullmatch(path.name)
            if m is None:
                raise ProtocolRefusal(f"unrecognised file {vdir.name}/{path.name}")
            frame = pd.read_parquet(path)
            if set(frame.columns) != {"idx", "score"}:
                raise ProtocolRefusal(
                    f"{vdir.name}/{path.name} must have exactly columns idx and score"
                )
            idx = frame["idx"].to_numpy()
            if idx.dtype.kind not in "iu" or frame["score"].dtype.kind not in "fiu":
                raise ProtocolRefusal(f"{vdir.name}/{path.name} has non-numeric columns")
            if len(set(idx.tolist())) != len(idx):
                raise ProtocolRefusal(f"{vdir.name}/{path.name} repeats a record index")
            runs.append(
                VariantRun(
                    variant=vdir.name,
                    role=role,
                    seed=int(m.group(1)),
                    scores={
                        int(i): float(s)
                        for i, s in zip(idx.tolist(), frame["score"].tolist())
                    },
                    config=config,
                    matched_to=None if matched is None else str(matched),
                )
            )
    return record, runs


# ==========================================================================
# Validation: every refusal condition lives here
# ==========================================================================
@dataclass
class _Plan:
    by_variant: dict[str, dict[int, VariantRun]]
    role: dict[str, str]
    seeds: list[int]
    control_of: dict[tuple[str, str], str]      # (subject, role) -> control variant
    subjects: list[str]                          # models first, then the baseline
    baseline: str | None
    structural: dict[str, dict[str, Any]]


def _short(h: str) -> str:
    return h[:12] if h else "(none)"


def _validate(
    record: SplitRecord,
    runs: Sequence[VariantRun],
    split_hash: str,
    n_parameters_rel_tol: float,
    expected_records_digest: str | None,
) -> _Plan:
    """Raise ProtocolRefusal unless a bias-controlled comparison is possible.

    Conditions, in order: the split is the requested one and is internally
    consistent; its records (allele, label, peptide) match their digest, and the
    caller's pinned digest if one was given; no peptide occurs in more than one of
    train, val and test; the runs are well formed; every model (and the baseline)
    has an HLA-only and a peptide-only control; the seed sets of a model and its
    controls are equal; the seeds are exactly EXPECTED_SEEDS; the controls' budget
    matches the model's (FR-62); every run scores exactly the split's test records;
    every score is finite; the controls satisfy their structural definition; the
    test set has both classes.
    """
    if record.content_hash != split_hash:
        raise ProtocolRefusal(
            f"split hash differs: the run directory holds split "
            f"{_short(record.content_hash)} but {_short(split_hash)} was requested"
        )
    try:
        recomputed = record.recomputed_hash()
    except (TypeError, ValueError) as exc:
        raise ProtocolRefusal(f"split.json cannot be reconstructed: {exc}") from exc
    if recomputed != record.content_hash:
        raise ProtocolRefusal(
            "split.json does not match its declared content hash "
            f"(recomputed {_short(recomputed)}, declared {_short(record.content_hash)})"
        )
    n = record.n_records
    if not (len(record.alleles) == len(record.labels) == len(record.peptides) == n):
        raise ProtocolRefusal("split.json records do not cover n_records")
    if record.recomputed_records_digest() != record.records_sha256:
        raise ProtocolRefusal(
            "records digest does not match: an allele, label or peptide in split.json "
            "differs from the records it was computed over (tampered or corrupted)"
        )
    if expected_records_digest is not None and (
        expected_records_digest != record.records_sha256
    ):
        raise ProtocolRefusal(
            f"records digest {_short(record.records_sha256)} differs from the pinned "
            f"digest {_short(expected_records_digest)}"
        )
    every = (*record.train_idx, *record.val_idx, *record.test_idx)
    if any(i < 0 or i >= n for i in every):
        raise ProtocolRefusal("split.json indexes records outside the corpus")
    if any(y not in (0, 1) for y in record.labels):
        raise ProtocolRefusal("split.json labels must be binary")
    if not record.test_idx or not record.train_idx:
        raise ProtocolRefusal("split has an empty train or test partition")

    seen: dict[str, str] = {}
    for part, indices in (
        ("train", record.train_idx),
        ("val", record.val_idx),
        ("test", record.test_idx),
    ):
        for peptide in {record.peptides[i] for i in indices}:
            if peptide in seen and seen[peptide] != part:
                raise ProtocolRefusal(
                    f"a peptide occurs in more than one partition ({seen[peptide]} "
                    f"and {part}); the split leaks and no result on it is valid"
                )
            seen[peptide] = part

    admitted = set(record.admitted_alleles)
    if any(record.alleles[i] not in admitted for i in record.test_idx):
        raise ProtocolRefusal("a test record belongs to an allele that was not admitted")
    test_labels = {record.labels[i] for i in record.test_idx}
    if test_labels != {0, 1}:
        raise ProtocolRefusal("the test set contains a single class; AUPRC is undefined")

    if not runs:
        raise ProtocolRefusal("no runs were supplied")
    by_variant: dict[str, dict[int, VariantRun]] = {}
    for run in runs:
        seeds = by_variant.setdefault(run.variant, {})
        if run.seed in seeds:
            raise ProtocolRefusal(f"variant {run.variant!r} repeats seed {run.seed}")
        seeds[run.seed] = run

    role: dict[str, str] = {}
    for name, seeds in by_variant.items():
        first = next(iter(seeds.values()))
        for r in seeds.values():
            if (r.role, r.config, r.matched_to) != (
                first.role, first.config, first.matched_to
            ):
                raise ProtocolRefusal(
                    f"variant {name!r} has a role, config or match that differs by seed"
                )
        if first.role not in ROLES:
            raise ProtocolRefusal(f"variant {name!r} has unknown role {first.role!r}")
        if (first.role in CONTROL_ROLES) != (first.matched_to is not None):
            raise ProtocolRefusal(
                f"variant {name!r}: only controls name the variant they are matched to"
            )
        role[name] = first.role

    subjects = sorted(n_ for n_, r in role.items() if r == "model")
    baselines = sorted(n_ for n_, r in role.items() if r == "baseline")
    if not subjects:
        raise ProtocolRefusal("no model variant was supplied")
    if len(baselines) > 1:
        raise ProtocolRefusal(f"at most one baseline variant is allowed; got {baselines}")
    everything = subjects + baselines

    control_of: dict[tuple[str, str], str] = {}
    for name, kind in sorted(role.items()):
        if kind not in CONTROL_ROLES:
            continue
        target = next(iter(by_variant[name].values())).matched_to
        if target not in everything:
            raise ProtocolRefusal(
                f"control {name!r} is matched to {target!r}, which is not a model "
                "or baseline variant"
            )
        key = (str(target), kind)
        if key in control_of:
            raise ProtocolRefusal(
                f"{target!r} has more than one {kind} control: "
                f"{control_of[key]!r}, {name!r}"
            )
        control_of[key] = name

    for subject in everything:
        for needed in CONTROL_ROLES:
            if (subject, needed) not in control_of:
                raise ProtocolRefusal(
                    f"variant {subject!r} has no {needed} control; the mandated null "
                    "controls cannot be skipped (FR-62, FR-63)"
                )

    for subject in everything:
        want = sorted(by_variant[subject])
        for needed in CONTROL_ROLES:
            control = control_of[(subject, needed)]
            got = sorted(by_variant[control])
            if got != want:
                raise ProtocolRefusal(
                    f"seed sets differ: {subject!r} reports seeds {want} but its "
                    f"{needed} control {control!r} reports seeds {got}"
                )
    seed_sets = {tuple(sorted(by_variant[s])) for s in everything}
    if len(seed_sets) != 1:
        raise ProtocolRefusal(
            f"seed sets differ between model and baseline variants: {sorted(seed_sets)}"
        )
    seeds_list = list(next(iter(seed_sets)))
    if seeds_list != sorted(EXPECTED_SEEDS):
        raise ProtocolRefusal(
            f"seeds {seeds_list} supplied; the pre-registered seeds are exactly "
            f"{sorted(EXPECTED_SEEDS)} (analysis plan section 5)"
        )

    for subject in everything:
        mc = next(iter(by_variant[subject].values())).config
        for needed in CONTROL_ROLES:
            control = control_of[(subject, needed)]
            cc = next(iter(by_variant[control].values())).config
            for fld in ("epochs", "early_stopping_metric", "tuning_trials"):
                if getattr(mc, fld) != getattr(cc, fld):
                    raise ProtocolRefusal(
                        f"control {control!r} does not match {subject!r} in {fld} "
                        f"({getattr(cc, fld)!r} vs {getattr(mc, fld)!r}); FR-62"
                    )
            slack = n_parameters_rel_tol * max(mc.n_parameters, 1)
            if abs(cc.n_parameters - mc.n_parameters) > slack:
                raise ProtocolRefusal(
                    f"control {control!r} does not match {subject!r} in n_parameters "
                    f"({cc.n_parameters} vs {mc.n_parameters}, relative tolerance "
                    f"{n_parameters_rel_tol}); FR-62"
                )

    test_set = set(record.test_idx)
    for name, seeds in sorted(by_variant.items()):
        for seed, run in sorted(seeds.items()):
            have = set(run.scores)
            missing, extra = test_set - have, have - test_set
            if missing and not extra:
                raise ProtocolRefusal(
                    f"variant {name!r} seed {seed} is missing a score for "
                    f"{len(missing)} test record(s)"
                )
            if missing or extra:
                raise ProtocolRefusal(
                    f"test indices differ from the split's test set in variant "
                    f"{name!r} seed {seed}: {len(missing)} missing, {len(extra)} outside"
                )
            values = np.fromiter((run.scores[i] for i in record.test_idx), dtype=float)
            if not np.all(np.isfinite(values)):
                raise ProtocolRefusal(
                    f"variant {name!r} seed {seed} has non-finite (NaN or inf) "
                    "predictions"
                )

    structural = _check_structure(record, by_variant, role)

    return _Plan(
        by_variant=by_variant,
        role=role,
        seeds=seeds_list,
        control_of=control_of,
        subjects=everything,
        baseline=baselines[0] if baselines else None,
        structural=structural,
    )


def _check_structure(
    record: SplitRecord,
    by_variant: Mapping[str, Mapping[int, VariantRun]],
    role: Mapping[str, str],
) -> dict[str, dict[str, Any]]:
    """Verify what the control labels claim, from the predictions themselves.

    hla_only: predictions are constant within every test allele (exact equality),
    for every seed. A model that never sees the peptide cannot do otherwise. A
    noise model labelled hla_only fails here.

    peptide_only: records of the test set that share a peptide get identical
    scores, for every seed. A model that never sees the allele cannot do
    otherwise. The check is vacuous when no two test records share a peptide, and
    the report says so.

    Neither check proves the control was built from the right inputs: an HLA-only
    lookup built from anything, or a peptide-only model that happens to ignore the
    allele on the shared peptides, passes. The matching block states this.
    """
    test_idx = list(record.test_idx)
    alleles = [record.alleles[i] for i in test_idx]
    names = sorted(set(alleles))
    pos = {a: j for j, a in enumerate(names)}
    allele_idx = np.array([pos[a] for a in alleles], dtype=np.int64)
    order = np.argsort(allele_idx, kind="stable")
    starts = np.flatnonzero(np.r_[True, np.diff(allele_idx[order]) != 0])

    by_peptide: dict[str, list[int]] = {}
    for j, i in enumerate(test_idx):
        by_peptide.setdefault(record.peptides[i], []).append(j)
    shared_groups = [js for js in by_peptide.values() if len(js) > 1]

    n_hla = n_pep = 0
    for name, seeds in sorted(by_variant.items()):
        kind = role[name]
        if kind not in CONTROL_ROLES:
            continue
        for seed, run in sorted(seeds.items()):
            values = np.fromiter((run.scores[i] for i in test_idx), dtype=float)
            if kind == "hla_only":
                n_hla += 1
                sv = values[order]
                hi = np.maximum.reduceat(sv, starts)
                lo = np.minimum.reduceat(sv, starts)
                varies = hi != lo
                if varies.any():
                    first_bad = int(np.argmax(varies))
                    allele = names[int(allele_idx[order][starts][first_bad])]
                    raise ProtocolRefusal(
                        f"structural check failed: {name!r} seed {seed} is labelled "
                        f"hla_only but its predictions vary within allele {allele!r}"
                    )
            else:
                n_pep += 1
                for js in shared_groups:
                    if np.ptp(values[js]) != 0.0:
                        raise ProtocolRefusal(
                            f"structural check failed: {name!r} seed {seed} is labelled "
                            "peptide_only but scores records that share a peptide "
                            "differently"
                        )
    note = (
        ""
        if shared_groups
        else "Vacuous: no two test records share a peptide, so this check cannot fail."
    )
    return {
        "hla_only": {
            "rule": "predictions are constant within every test allele (exact "
                    "equality), for every seed",
            "n_runs_checked": n_hla,
            "n_alleles": len(names),
        },
        "peptide_only": {
            "rule": "test records that share a peptide receive identical "
                    "predictions, for every seed",
            "n_runs_checked": n_pep,
            "n_shared_peptide_groups": len(shared_groups),
            "note": note,
        },
    }


# ==========================================================================
# Metrics (seed-free tie handling; see nbx.evaluate.stats)
# ==========================================================================
METRIC_NAMES = (
    "auprc",
    "auroc",
    "ppv_at_k",
    "recall_at_k",
    "fraction_ranked",
    "top20_immunogenic_fraction",
    "mcc",
)

_MetricFn = Callable[[Sequence[float], Sequence[int]], float]


def _arr(scores: Sequence[float], labels: Sequence[int]) -> tuple[np.ndarray, np.ndarray]:
    return np.asarray(scores, dtype=float), np.asarray(labels, dtype=np.int64)


def _metric_functions(k: int) -> dict[str, _MetricFn]:
    def auprc_(s: Sequence[float], y: Sequence[int]) -> float:
        return stats.expected_average_precision(*_arr(s, y))

    def auroc_(s: Sequence[float], y: Sequence[int]) -> float:
        return auroc(list(s), list(y))          # ties already average-ranked

    def ppv(s: Sequence[float], y: Sequence[int]) -> float:
        return stats.expected_precision_at_k(*_arr(s, y), k)

    def recall(s: Sequence[float], y: Sequence[int]) -> float:
        return stats.expected_recall_at_k(*_arr(s, y), k)

    def fraction_ranked_(s: Sequence[float], y: Sequence[int]) -> float:
        return stats.expected_recall_at_k(*_arr(s, y), 100)

    def top20(s: Sequence[float], y: Sequence[int]) -> float:
        a, b = _arr(s, y)
        return stats.expected_precision_at_k(a, b, min(20, a.size))

    def mcc(s: Sequence[float], y: Sequence[int]) -> float:
        return stats.expected_mcc_top_npos(*_arr(s, y))

    return {
        "auprc": auprc_,
        "auroc": auroc_,
        "ppv_at_k": ppv,
        "recall_at_k": recall,
        "fraction_ranked": fraction_ranked_,
        "top20_immunogenic_fraction": top20,
        "mcc": mcc,
    }


def _try(fn: Callable[[], float]) -> float | None:
    try:
        return float(fn())
    except ValueError:
        return None                       # undefined here: unmeasured, never zero


@dataclass
class _RunResult:
    pooled: dict[str, float | None]
    weighted: dict[str, float | None]
    macro: dict[str, float | None]
    n_alleles: dict[str, int]
    kappa: float | None
    shared_means: np.ndarray | None


def _score_run(
    raw: np.ndarray,
    labels: list[int],
    alleles: list[str],
    train_rho: Mapping[Any, float],
    shared: list[str],
    k: int,
) -> _RunResult:
    fns = _metric_functions(k)
    s_list = raw.tolist()
    pooled: dict[str, float | None] = {}
    for n, f in fns.items():
        pooled[n] = _try(functools.partial(f, s_list, labels))
    weighted: dict[str, float | None] = {}
    macro: dict[str, float | None] = {}
    n_alleles: dict[str, int] = {}
    for n, f in fns.items():
        weighted[n], macro[n], n_alleles[n] = None, None, 0
        try:
            per, agg = allele_stratified(f, s_list, labels, alleles, True)
            _, agg_macro = allele_stratified(f, s_list, labels, alleles, False)
        except ValueError:
            continue
        weighted[n], macro[n], n_alleles[n] = agg, agg_macro, len(per)

    kappa: float | None
    try:
        kappa = shortcut_dependence(train_rho, alleles, s_list)
    except ValueError:
        kappa = None
    means: np.ndarray | None = None
    if kappa is not None:
        sums: dict[str, float] = {}
        counts: dict[str, int] = {}
        for a, v in zip(alleles, s_list):
            if a in train_rho:
                sums[a] = sums.get(a, 0.0) + v
                counts[a] = counts.get(a, 0) + 1
        means = np.array([sums[a] / counts[a] for a in shared])
    return _RunResult(pooled, weighted, macro, n_alleles, kappa, means)


def _nan_last(value: float | None) -> float:
    """Sort helper: an unmeasured headline sorts below every measured one."""
    return -1.0 if value is None else value


def _mean_sd(values: Sequence[float | None]) -> tuple[float | None, float | None]:
    if not values or any(v is None for v in values):
        return None, None
    arr = np.array(values, dtype=float)
    return float(arr.mean()), float(arr.std(ddof=1)) if arr.size > 1 else None


def _r(value: float | None) -> str:
    return render_metric(value, DIGITS)


# ==========================================================================
# Orchestration
# ==========================================================================
def _config_hash(config: RunConfig) -> str:
    payload = json.dumps(asdict(config), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _interval_block(values: np.ndarray) -> dict[str, Any]:
    lo, hi, n = stats.percentile_interval(values, 0.90)
    return {
        "level": _r(0.90),
        "low": _r(lo),
        "high": _r(hi),
        "n_valid_resamples": n,
    }


def evaluate_runs(
    record: SplitRecord,
    runs: Sequence[VariantRun],
    *,
    split_hash: str,
    expected_records_digest: str | None = None,
    n_boot: int = PRODUCTION_N_BOOT,
    bootstrap_seed: int = DEFAULT_BOOTSTRAP_SEED,
    ranking_k: int = DEFAULT_RANKING_K,
    n_parameters_rel_tol: float = 0.0,
) -> dict[str, Any]:
    """Validate, then evaluate in-memory runs against a split. Returns the report.

    `n_boot` below the analysis plan's 10,000 is for tests; the report then says
    `preregistered_settings_used: false`. The command line cannot do this.
    `n_parameters_rel_tol` is the relative tolerance on the parameter count of a
    control against its model; zero (exact) unless a justification exists, and it
    is written into the report whenever it is non-zero.
    """
    if n_boot < 1 or ranking_k < 1:
        raise ValueError("n_boot and ranking_k must be positive")
    if not 0.0 <= n_parameters_rel_tol < 1.0:
        raise ValueError("n_parameters_rel_tol must lie in [0, 1)")
    plan = _validate(
        record, runs, split_hash, n_parameters_rel_tol, expected_records_digest
    )

    # ---- the test set --------------------------------------------------
    test_idx = list(record.test_idx)
    labels = [record.labels[i] for i in test_idx]
    alleles = [record.alleles[i] for i in test_idx]
    allele_names = sorted(set(alleles))
    allele_pos = {a: j for j, a in enumerate(allele_names)}
    allele_idx = np.array([allele_pos[a] for a in alleles], dtype=np.int64)
    labels_arr = np.array(labels, dtype=np.int64)
    counts = stats.cluster_counts(len(allele_names), n_boot, bootstrap_seed)

    train_rho = per_allele_imbalance(
        [record.alleles[i] for i in record.train_idx],
        [record.labels[i] for i in record.train_idx],
    )
    shared = sorted((a for a in allele_names if a in train_rho), key=str)
    rho_shared = np.array([train_rho[a] for a in shared], dtype=float)
    shared_cols = [allele_pos[a] for a in shared]
    counts_shared = counts[:, shared_cols]

    # ---- per-run scoring -----------------------------------------------
    scores: dict[tuple[str, int], np.ndarray] = {}
    results: dict[tuple[str, int], _RunResult] = {}
    for name, seeds in plan.by_variant.items():
        for seed, run in seeds.items():
            raw = np.array([run.scores[i] for i in test_idx], dtype=float)
            scores[(name, seed)] = raw
            results[(name, seed)] = _score_run(
                raw, labels, alleles, train_rho, shared, ranking_k
            )

    def agg_metric(
        name: str, group: str, metric: str
    ) -> tuple[float | None, float | None]:
        return _mean_sd(
            [getattr(results[(name, s)], group)[metric] for s in plan.seeds]
        )

    def kappa_stats(name: str) -> tuple[float | None, float | None]:
        return _mean_sd([results[(name, s)].kappa for s in plan.seeds])

    baseline_kappa: float | None = None
    if plan.baseline is not None:
        baseline_kappa = kappa_stats(plan.baseline)[0]

    def kappa_boot(name: str) -> np.ndarray | None:
        means = [results[(name, s)].shared_means for s in plan.seeds]
        if any(m is None for m in means) or len(shared) < 3:
            return None
        return stats.kappa_resampled(
            rho_shared, [m for m in means if m is not None], counts_shared
        )

    baseline_boot = kappa_boot(plan.baseline) if plan.baseline is not None else None

    # ---- margins -------------------------------------------------------
    margins: dict[tuple[str, str], stats.MarginEstimate] = {}
    for name in plan.subjects:
        for ctl_role in CONTROL_ROLES:
            ctl = plan.control_of[(name, ctl_role)]
            margins[(name, ctl_role)] = stats.margin_over_control(
                [scores[(name, s)] for s in plan.seeds],
                [scores[(ctl, s)] for s in plan.seeds],
                labels_arr,
                allele_idx,
                counts,
            )

    model_names = sorted(n for n in plan.subjects if plan.role[n] == "model")
    p_raw = [margins[(n, "hla_only")].p_value for n in model_names]
    holm: dict[str, float | None] = {n: None for n in model_names}
    if all(p is not None for p in p_raw):
        adjusted = stats.holm_adjust([float(p) for p in p_raw if p is not None])
        holm = dict(zip(model_names, adjusted))

    # ---- gate and rows -------------------------------------------------
    variants: dict[str, Any] = {}
    gate_results: dict[str, GateResult] = {}
    mean_auprc: dict[str, float | None] = {}

    for name in sorted(plan.by_variant):
        role = plan.role[name]
        cfg = next(iter(plan.by_variant[name].values())).config
        entry: dict[str, Any] = {
            "role": role,
            "config": asdict(cfg),
            "config_hash": _config_hash(cfg),
            "seeds": list(plan.seeds),
        }
        mean_auprc[name] = agg_metric(name, "pooled", "auprc")[0]
        kap_mean, kap_sd = kappa_stats(name)

        if role in CONTROL_ROLES:
            gate = evaluate_gate(is_control=True)
            gate_results[name] = gate
            entry["matched_to"] = next(iter(plan.by_variant[name].values())).matched_to
            entry["gate"] = _gate_dict(gate)
            entry["rank"] = UNMEASURED
            entry["metrics"] = _metrics_block(agg_metric, name, "control row")
            entry["shortcut_dependence"] = {
                "kappa_mean": _r(kap_mean),
                "kappa_sd": _r(kap_sd),
                "note": "null control; shortcut dependence of a control is context only",
            }
            variants[name] = entry
            continue

        est_h = margins[(name, "hla_only")]
        est_p = margins[(name, "peptide_only")]
        cm_h = ControlMargin(delta=est_h.delta, p_value=est_h.p_value, alpha=ALPHA)
        cm_p = ControlMargin(delta=est_p.delta, p_value=est_p.p_value, alpha=ALPHA)

        boot = kappa_boot(name)
        change: np.ndarray | None = None
        upper: float | None = None
        harness_note = ""
        if role == "baseline":
            upper = 0.0                           # compared with itself: vacuous
        elif plan.baseline is None or baseline_boot is None:
            harness_note = (
                "No baseline run (or its kappa is unmeasured): the change in "
                "shortcut dependence cannot be estimated, so the kappa interval is "
                "withheld from the gate and no rank can be assigned."
            )
        elif boot is not None:
            change = np.abs(boot) - np.abs(baseline_boot)
            upper = stats.one_sided_bound(change, "upper")
        gate = evaluate_gate(
            headline_metric=mean_auprc[name],
            hla_only_margin=cm_h,
            peptide_only_margin=cm_p,
            kappa_change_upper_bound=upper,
            peptide_only_significantly_worse=est_p.significantly_worse,
        )
        gate_results[name] = gate
        entry["gate"] = _gate_dict(gate)
        entry["rank"] = UNMEASURED                # assigned below, only if PASSED
        if harness_note:
            entry["harness_note"] = harness_note
        if role == "baseline":
            entry["reference_only"] = True
            entry["rank_note"] = (
                "Reference baseline: its kappa is compared with itself, so the "
                "shortcut condition is vacuous for this row. Never ranked."
            )

        if not gate.displayable_metrics:
            entry["metrics"] = {"displayed": False, "note": FR63}
            entry["margin_vs_hla_only"] = {"rendered": UNMEASURED,
                                           "note": "withheld while PENDING"}
            entry["margin_vs_peptide_only"] = {"rendered": UNMEASURED,
                                               "note": "withheld while PENDING"}
            entry["shortcut_dependence"] = {"displayed": False}
            variants[name] = entry
            continue

        entry["metrics"] = _metrics_block(agg_metric, name, cm_h.render())
        entry["margin_vs_hla_only"] = _margin_block(est_h, cm_h)
        if role == "model":
            entry["margin_vs_hla_only"]["p_value_holm"] = _r(holm[name])
        entry["margin_vs_peptide_only"] = _margin_block(est_p, cm_p)
        entry["margin_vs_peptide_only"].update(
            {
                "p_value_worse": _r(est_p.p_value_worse),
                "significantly_worse": est_p.significantly_worse,
                "note": "gates only through condition (iv): the model is "
                        "significantly worse than the peptide-only control",
            }
        )
        kappa_block: dict[str, Any] = {
            "kappa_mean": _r(kap_mean),
            "kappa_sd": _r(kap_sd),
            "baseline_kappa": _r(baseline_kappa),
        }
        if boot is not None:
            kappa_block["kappa_ci90"] = _interval_block(boot)
        if change is not None:
            lo, hi, n_valid = stats.percentile_interval(change, 0.90)
            kappa_block["abs_kappa_change"] = {
                "upper95_one_sided": _r(upper),
                "ci90_low": _r(lo),
                "ci90_high": _r(hi),
                "noninferiority_margin": _r(NONINFERIORITY_MARGIN_KAPPA),
                "n_valid_resamples": n_valid,
                "definition": "|kappa_model| - |kappa_baseline|, allele-clustered "
                              "bootstrap, seed-averaged within each resample",
            }
        elif role == "baseline":
            kappa_block["abs_kappa_change"] = {
                "upper95_one_sided": _r(upper),
                "note": "baseline compared with itself",
            }
        else:
            kappa_block["abs_kappa_change"] = {
                "upper95_one_sided": UNMEASURED,
                "ci90_low": UNMEASURED,
                "ci90_high": UNMEASURED,
                "noninferiority_margin": _r(NONINFERIORITY_MARGIN_KAPPA),
                "n_valid_resamples": 0,
            }
        entry["shortcut_dependence"] = kappa_block
        variants[name] = entry

    # ---- ranking: only variants the gate made rankable ------------------
    ranked = sorted(
        (n for n, g in gate_results.items() if g.rankable and plan.role[n] == "model"),
        key=lambda n: (-_nan_last(mean_auprc[n]), n),
    )
    for position, n in enumerate(ranked, start=1):
        variants[n]["rank"] = position

    # ---- alleles outside the admitted test set --------------------------
    admitted_set = set(record.admitted_alleles)
    outside: dict[str, dict[str, int]] = {}
    for a, y in zip(record.alleles, record.labels):
        if a not in admitted_set:
            c = outside.setdefault(a, {"records": 0, "positive": 0})
            c["records"] += 1
            c["positive"] += y

    preregistered = (
        n_boot == PRODUCTION_N_BOOT
        and bootstrap_seed == DEFAULT_BOOTSTRAP_SEED
        and ranking_k == DEFAULT_RANKING_K
        and n_parameters_rel_tol == 0.0
    )

    report: dict[str, Any] = {
        "protocol": "bias-controlled",
        "schema_version": 2,
        "notice": "Research use only. Not for clinical or diagnostic use.",
        "split_hash": record.content_hash,
        "records": {
            "sha256": record.records_sha256,
            "verified_against_split_json": True,
            "pinned_by_caller": expected_records_digest is not None,
            "covers": "per-record allele, label and peptide of the corpus",
            "note": "The split hash alone does not cover labels, alleles or peptides. "
                    "The digest was recomputed from split.json; it detects "
                    "inconsistency, and detects a coordinated edit only when pinned.",
        },
        "seeds": list(plan.seeds),
        "preregistered_settings_used": preregistered,
        "headline_metric": "auprc",
        "alpha": _r(ALPHA),
        "bootstrap": {
            "B": n_boot,
            "seed": bootstrap_seed,
            "unit": "test alleles (clusters), resampled with replacement; the same "
                    "resamples are used for every comparison and every seed",
            "test": "paired by seed and by resample; one-sided against delta <= 0",
            "p_value": "(1 + #{b : delta_b <= 0}) / (B_valid + 1)",
            "interval": "two-sided 90% percentile interval (5th and 95th "
                        "percentiles), the interval that matches the one-sided test "
                        "at alpha 0.05; one_sided_lower95 is its lower end",
            "rng": f"numpy.random.default_rng (PCG64), numpy {np.__version__}",
            "tie_handling": "expected average precision over uniformly random "
                            "orderings within groups of tied scores (exact, equal to "
                            "metrics.auprc when nothing ties); ranking metrics are "
                            "expectations the same way; no seed is involved, and the "
                            "bootstrap seed seeds the resampling only",
        },
        "matching": {
            "declared_metadata": ["epochs", "n_parameters", "early_stopping_metric",
                                  "tuning_trials", "role", "matched_to"],
            "declared_metadata_note": "Compared as written in config.json. They are "
                                      "claims by the trainer, not observations.",
            "structurally_verified": plan.structural,
            "not_verified": "Capacity, epochs, early stopping and tuning effort are "
                            "declared and compared but not verified against the "
                            "training itself. The structural rules catch a control "
                            "whose predictions break its definition; they do not "
                            "prove it was built from the right inputs.",
            "n_parameters_rel_tol": _r(n_parameters_rel_tol),
        },
        "weighting": {
            "pooled": "all admitted test records pooled; every record has equal "
                      "weight, so alleles with more test records count for more",
            "allele_stratified_support_weighted": (
                "per-allele metric averaged with weight proportional to the allele's "
                "number of test records (support); alleles where the metric is "
                "undefined are omitted, see n_alleles_defined"
            ),
            "allele_stratified_macro": (
                "per-allele metric averaged with equal weight per allele; alleles "
                "where the metric is undefined are omitted, see n_alleles_defined"
            ),
        },
        "ranking_k": ranking_k,
        "mcc_rule": "predict positive for the n_pos highest scores, n_pos being the "
                    "number of positives in the evaluated set; computed from the "
                    "expected confusion counts when scores tie",
        "aggregation": "mean and across-seed standard deviation (ddof=1); the best "
                       "seed is never shown",
        "n_model_variants_compared": len(model_names),
        "multiplicity": (
            f"{len(model_names)} model variant(s) compared in this run. The HLA-only "
            "margin p-value of each is also reported Holm-adjusted across that family "
            "(p_value_holm). The gate reads the unadjusted one-sided p at alpha 0.05; "
            "any claim made across several variants should use the adjusted value."
        ),
        "kappa_gate_rule": (
            "Gate rule in force (decision D4), burden of proof on the model: (i) "
            "significant margin over the HLA-only control; (ii) non-inferiority: the "
            "upper 95% bound of |kappa_model| - |kappa_baseline| is at most "
            f"{_r(NONINFERIORITY_MARGIN_KAPPA)}; (iv) not significantly worse than "
            "the peptide-only control. (iii), an absolute ceiling on |kappa|, is NOT "
            "implemented: its value is the researcher's decision."
        ),
        "test_set": {
            "n_records": len(test_idx),
            "n_alleles": len(allele_names),
            "alleles": allele_names,
        },
        "not_admitted": {
            "alleles": [
                {"allele": a, "n_records": c["records"], "n_positive": c["positive"]}
                for a, c in sorted(outside.items())
            ],
            "metrics": UNMEASURED,
            "note": "Alleles outside the admitted test set are listed for "
                    "transparency; no metric is computed on them and they are never "
                    "pooled into any figure above.",
        },
        "variants": variants,
        "ranking": ranked,
        "ranking_note": "Rank order is by mean AUPRC among variants that passed the "
                        "gate; differences between ranked variants are not tested.",
        "not_done": list(_NOT_DONE),
    }
    report["report_hash"] = _hash_report(report)
    return report


def _gate_dict(gate: GateResult) -> dict[str, Any]:
    return {
        "status": gate.status.value,
        "reason": gate.reason,
        "rankable": gate.rankable,
        "displayable_metrics": gate.displayable_metrics,
    }


def _margin_block(est: stats.MarginEstimate, margin: ControlMargin) -> dict[str, Any]:
    return {
        "rendered": margin.render(),
        "significant": margin.significant,
        "delta": _r(est.delta),
        "p_value": _r(est.p_value),
        "ci90_low": _r(est.ci90_low),
        "ci90_high": _r(est.ci90_high),
        "one_sided_lower95": _r(est.one_sided_lower95),
        "per_seed_delta_mean": _r(est.per_seed_mean),
        "per_seed_delta_sd": _r(est.per_seed_sd),
        "n_valid_resamples": est.n_valid_resamples,
        "n_seeds": est.n_seeds,
    }


def _metrics_block(
    agg: Callable[..., Any], name: str, margin_text: str
) -> dict[str, Any]:
    block: dict[str, Any] = {"displayed": True}
    groups = (
        ("pooled", "pooled"),
        ("allele_stratified_support_weighted", "weighted"),
        ("allele_stratified_macro", "macro"),
    )
    for out_key, attr in groups:
        rows: dict[str, Any] = {}
        for metric in METRIC_NAMES:
            mean, sd = agg(name, attr, metric)
            rows[metric] = {
                "mean": _r(mean),
                "sd": _r(sd),
                "margin_over_hla_only": margin_text,
            }
        block[out_key] = rows
    return block


def _hash_report(report: Mapping[str, Any]) -> str:
    body = {k: v for k, v in report.items() if k != "report_hash"}
    payload = json.dumps(body, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def report_to_json(report: Mapping[str, Any]) -> str:
    """Canonical text of a report: sorted keys, two-space indent, UTF-8 (em dashes)."""
    return json.dumps(report, sort_keys=True, indent=2, ensure_ascii=False) + "\n"


def write_report(report: Mapping[str, Any], path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(report_to_json(report).encode("utf-8"))
    return path


def run_protocol(
    run_dir: Path,
    split_hash: str,
    *,
    expected_records_digest: str | None = None,
    n_boot: int = PRODUCTION_N_BOOT,
    bootstrap_seed: int = DEFAULT_BOOTSTRAP_SEED,
    ranking_k: int = DEFAULT_RANKING_K,
    n_parameters_rel_tol: float = 0.0,
) -> dict[str, Any]:
    """Execute the full bias-controlled evaluation for one ablation matrix.

    Reads `<run_dir>` (see the module docstring), refuses with `ProtocolRefusal`
    unless every precondition holds, and returns the content-hashed report.
    """
    record, runs = load_run_dir(Path(run_dir))
    return evaluate_runs(
        record,
        runs,
        split_hash=split_hash,
        expected_records_digest=expected_records_digest,
        n_boot=n_boot,
        bootstrap_seed=bootstrap_seed,
        ranking_k=ranking_k,
        n_parameters_rel_tol=n_parameters_rel_tol,
    )
