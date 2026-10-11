"""Corpus-level bias report (layer L2 output; FR-20, FR-21, FR-22, FR-24).

This module only assembles the report. The mathematics (mutual information,
per-allele imbalance, per-position information, the masking policy, split
construction) lives in `nbx.bias.diagnostics`, `nbx.bias.policy` and
`nbx.splits.freeze`, where it is pinned by analytic tests. Nothing is
recomputed here.

Three properties are deliberate:

* A quantity that cannot be computed is stored as `UNMEASURED` (an em dash),
  never as zero and never as null.
* Mutual information is reported as a bias CANDIDATE signal. A variable that
  carries information about the label is not thereby a shortcut; a genuine
  determinant also carries information.
* The split admissibility block is computed in memory and is never frozen here.
  Freezing is `nbx split freeze`, and a frozen split is immutable.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from nbx.bias.diagnostics import entropy, profile_dataset
from nbx.bias.policy import (
    adaptive_policy,
    per_position_mi,
    select_threshold,
    uniform_policy,
)
from nbx.evaluate.gate import UNMEASURED
from nbx.splits.freeze import SplitCriterion, build_split

__all__ = ["BiasReportError", "build_bias_report", "write_bias_report"]

REQUIRED_COLUMNS = (
    "peptide",
    "peptide_length",
    "allele",
    "label",
    "source_dataset",
    "assay_type",
)
POLICY_PEPTIDE_LENGTH = 9
POLICY_TAU_QUANTILE = 0.90
TAU_SELECTED_ON = "full_corpus_no_validation_split"


class BiasReportError(ValueError):
    """The corpus cannot be profiled, or the report cannot be written."""


def _validate(df: pd.DataFrame) -> None:
    if len(df) == 0:
        raise BiasReportError("the corpus is empty; there is nothing to profile")
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise BiasReportError(f"required column(s) missing: {', '.join(missing)}")
    if not set(df["label"].unique()) <= {0, 1}:
        raise BiasReportError("labels must be binary 0/1")
    if df["label"].nunique() < 2:
        raise BiasReportError(
            "the corpus must contain both classes; a single-class corpus has "
            "no label entropy to explain"
        )


def _masking_policy_block(df: pd.DataFrame, min_records: int) -> dict:
    scope = df[df["peptide_length"] == POLICY_PEPTIDE_LENGTH]
    counts = scope.groupby("allele").size()
    scope = scope[scope["allele"].isin(counts[counts >= min_records].index)]
    block: dict = {
        "peptide_length": POLICY_PEPTIDE_LENGTH,
        "min_records_per_allele": min_records,
        "tau_quantile": POLICY_TAU_QUANTILE,
        "tau_selected_on": TAU_SELECTED_ON,
        "note": (
            "tau is chosen on the corpus being profiled because no frozen "
            "validation split exists yet; plug-in mutual information is inflated "
            "for small allele subsets; treat this block as exploratory"
        ),
        "n_alleles_in_scope": int(scope["allele"].nunique()),
        "n_alleles_total": int(df["allele"].nunique()),
        "n_records_in_scope": len(scope),
        "n_records_total": len(df),
    }
    if scope.empty:
        block.update(
            tau_nats=UNMEASURED,
            n_alleles_nonempty_mask=UNMEASURED,
            n_alleles_mask_equals_uniform=UNMEASURED,
            n_alleles_mask_includes_position_2=UNMEASURED,
            n_alleles_mask_includes_c_terminus=UNMEASURED,
        )
        return block

    mi_by_allele = per_position_mi(
        [list(p) for p in scope["peptide"]],
        scope["label"].tolist(),
        scope["allele"].tolist(),
        POLICY_PEPTIDE_LENGTH,
    )
    tau = select_threshold(mi_by_allele, POLICY_TAU_QUANTILE)
    adaptive = adaptive_policy(mi_by_allele, tau, POLICY_PEPTIDE_LENGTH)
    uniform = uniform_policy(POLICY_PEPTIDE_LENGTH).default
    p2, c_term = min(uniform), max(uniform)
    masks = list(adaptive.per_allele.values())
    # A bare "differs from uniform" count is not reported: an empty mask never
    # equals the uniform mask, so that count is inflated whenever tau is high.
    block.update(
        tau_nats=tau,
        n_alleles_nonempty_mask=sum(len(m) > 0 for m in masks),
        n_alleles_mask_equals_uniform=sum(m == uniform for m in masks),
        n_alleles_mask_includes_position_2=sum(p2 in m for m in masks),
        n_alleles_mask_includes_c_terminus=sum(c_term in m for m in masks),
    )
    return block


def build_bias_report(
    df: pd.DataFrame,
    *,
    input_path: str | None = None,
    input_sha256: str | None = None,
    min_records_per_allele: int = 20,
    split_criterion: SplitCriterion | None = None,
    extra_variables: Sequence[str] = (),
) -> dict:
    """Profile a harmonised corpus. Every figure is computed from `df`."""
    _validate(df)
    labels = df["label"].astype(int).tolist()

    variables = {
        "hla_allele": df["allele"].tolist(),
        "source_dataset": df["source_dataset"].tolist(),
        "assay_type": df["assay_type"].tolist(),
        "peptide_length": df["peptide_length"].tolist(),
        "position_2_residue": df["peptide"].str[1].tolist(),
        "c_terminal_residue": df["peptide"].str[-1].tolist(),
    }
    for col in extra_variables:
        if col not in df.columns:
            raise BiasReportError(f"required column(s) missing: {col}")
        variables[col] = df[col].tolist()

    profile = profile_dataset(labels, variables, allele_key="hla_allele")
    h_y = entropy(labels)

    mi_block = {}
    for name, nats in profile.mutual_information.items():
        n_distinct = len(set(variables[name]))
        mi_block[name] = {
            "nats": nats,
            "share_of_label_entropy": nats / h_y,
            "n_distinct_values": n_distinct,
            "structurally_uninformative": n_distinct == 1,
        }

    rho = sorted(profile.per_allele_rho.values())
    rho_summary: dict[str, float | str]
    if rho:
        rho_summary = {
            "min": rho[0],
            "median": rho[len(rho) // 2]
            if len(rho) % 2
            else (rho[len(rho) // 2 - 1] + rho[len(rho) // 2]) / 2,
            "max": rho[-1],
        }
    else:
        rho_summary = {"min": UNMEASURED, "median": UNMEASURED, "max": UNMEASURED}
    per_allele_block = {
        "n_alleles_with_rho": len(rho),
        "n_excluded_single_sample": profile.n_alleles_excluded_single_sample,
        "n_beyond_10_to_1": len(profile.imbalanced_alleles),
        "n_saturated": sum(abs(v) >= 10.0 for v in rho),
        "n_positive_rho": sum(v > 0 for v in rho),
        "n_negative_rho": sum(v < 0 for v in rho),
        **rho_summary,
    }

    criterion = split_criterion or SplitCriterion(
        min_samples_per_allele=min_records_per_allele
    )
    # Peptide-grouped, as pre-registered (docs/analysis_plan.md section 3): every record
    # of a peptide lies in one partition, and records withheld to guarantee that are
    # quarantined and counted, never moved.
    split = build_split(
        df["allele"].tolist(), labels, criterion, groups=df["peptide"].tolist()
    )
    split_block = {
        "criterion": json.loads(split.criterion.canonical_json()),
        "group_by": split.criterion.group_by,
        "n_quarantined": len(split.quarantined_idx),
        "n_admitted_alleles": len(split.admitted_alleles),
        "n_reported_separately": len(split.reported_separately),
        "admitted_alleles": list(split.admitted_alleles),
        "reported_separately": list(split.reported_separately),
        "n_test_records": len(split.test_idx),
        "frozen": False,
    }

    return {
        "input_path": input_path if input_path is not None else UNMEASURED,
        "input_file_sha256": input_sha256 if input_sha256 is not None else UNMEASURED,
        "n_records": len(df),
        "n_alleles": int(df["allele"].nunique()),
        "label_counts": {
            str(k): int(v) for k, v in df["label"].value_counts().sort_index().items()
        },
        "source_breakdown": {
            str(k): int(v) for k, v in df["source_dataset"].value_counts().items()
        },
        "length_distribution": {
            str(k): int(v)
            for k, v in df["peptide_length"].value_counts().sort_index().items()
        },
        "label_entropy_nats": h_y,
        "mutual_information": mi_block,
        "bias_candidates_gt_25pct_of_label_entropy": sorted(profile.bias_candidates()),
        "per_allele_rho": per_allele_block,
        "masking_policy": _masking_policy_block(df, min_records_per_allele),
        "split_admissibility": split_block,
        "caveats": [
            "mutual information identifies bias candidates, not shortcuts",
            "plug-in mutual information estimates are biased upward for small groups",
            "a single-valued variable has zero mutual information by construction",
        ],
    }


def write_bias_report(report: dict, path: Path) -> Path:
    """Write the report as canonical JSON; never into a frozen-splits directory."""
    path = Path(path)
    if "splits" in {p.lower() for p in path.parts}:
        raise BiasReportError(
            "refusing to write under a 'splits' directory: frozen splits are immutable"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path
