"""Tests for the corpus-level bias report (FR-20, FR-21, FR-22, FR-24, FR-60/61).

All fixtures are small synthetic frames built here; none of their values is a
measurement. The mathematics itself is pinned in `test_bias_diagnostics.py`; these
tests pin that the report USES it, accounts for every allele, and refuses
degenerate input.
"""

from __future__ import annotations

import json
import math

import pandas as pd
import pytest

from nbx.bias.diagnostics import entropy, mutual_information
from nbx.bias.profile import (
    BiasReportError,
    build_bias_report,
    write_bias_report,
)
from nbx.evaluate.gate import UNMEASURED


def _row(peptide, allele, label, source="S1", assay="A1"):
    return {
        "peptide": peptide,
        "peptide_length": len(peptide),
        "allele": allele,
        "label": label,
        "source_dataset": source,
        "assay_type": assay,
    }


def _frame(per_allele_counts: dict[str, tuple[int, int]], source="S1") -> pd.DataFrame:
    """One allele per entry: (n_positive, n_negative). Nine-mers with varied residues."""
    rows, k = [], 0
    for allele, (npos, nneg) in per_allele_counts.items():
        for label, n in ((1, npos), (0, nneg)):
            for _ in range(n):
                aa = "ACDEFGHIKLMNPQRSTVWY"
                pep = "".join(aa[(k * 7 + j * 3) % 20] for j in range(9))
                rows.append(_row(pep, allele, label, source=source))
                k += 1
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# FR-20: mutual information, as a value and as a share of H(y)
# --------------------------------------------------------------------------
def test_fr20_report_uses_the_tested_mutual_information_function():
    df = _frame(
        {"HLA-A*01:01": (30, 10), "HLA-A*02:01": (5, 35), "HLA-B*07:02": (20, 20)}
    )
    rep = build_bias_report(df)
    expected = mutual_information(df["label"].tolist(), df["allele"].tolist())
    got = rep["mutual_information"]["hla_allele"]
    assert got["nats"] == pytest.approx(expected)
    assert got["share_of_label_entropy"] == pytest.approx(
        expected / entropy(df["label"].tolist())
    )


def test_fr20_a_variable_that_determines_the_label_carries_all_of_it():
    # Analytic case: allele fully determines the label, so I(b;y) = H(y).
    df = _frame({"HLA-A*01:01": (40, 0), "HLA-A*02:01": (0, 40)})
    rep = build_bias_report(df)
    assert rep["mutual_information"]["hla_allele"]["share_of_label_entropy"] == (
        pytest.approx(1.0)
    )


def test_fr20_constant_variable_is_flagged_not_silently_reported_as_clean():
    df = _frame({"HLA-A*01:01": (30, 10), "HLA-A*02:01": (5, 35)}, source="ONLY")
    rep = build_bias_report(df)
    sd = rep["mutual_information"]["source_dataset"]
    assert sd["nats"] == pytest.approx(0.0)
    # A single-valued variable cannot show bias; the report must say that the
    # zero is structural, so nobody reads it as "no source bias measured".
    assert sd["n_distinct_values"] == 1
    assert sd["structurally_uninformative"] is True


def test_fr20_candidates_are_never_called_shortcuts():
    df = _frame({"HLA-A*01:01": (40, 0), "HLA-A*02:01": (0, 40)})
    text = json.dumps(build_bias_report(df)).lower()
    assert "bias_candidates" in text
    assert "is_a_shortcut" not in text and '"shortcut"' not in text


# --------------------------------------------------------------------------
# FR-21: per-allele imbalance, saturation, single-sample exclusion
# --------------------------------------------------------------------------
def test_fr21_counts_alleles_beyond_ten_to_one_and_saturated():
    df = _frame(
        {
            "HLA-A*01:01": (30, 30),    # rho = 0
            "HLA-A*02:01": (44, 4),     # 11:1  -> beyond 10:1
            "HLA-B*07:02": (30, 0),     # only positives -> saturated
            "HLA-B*08:01": (0, 25),     # only negatives -> saturated
        }
    )
    d = build_bias_report(df)["per_allele_rho"]
    assert d["n_alleles_with_rho"] == 4
    assert d["n_beyond_10_to_1"] == 3
    assert d["n_saturated"] == 2


def test_fr21_single_sample_alleles_are_excluded_and_counted():
    df = _frame({"HLA-A*01:01": (30, 30), "HLA-C*01:02": (1, 0)})
    d = build_bias_report(df)["per_allele_rho"]
    assert d["n_alleles_with_rho"] == 1
    assert d["n_excluded_single_sample"] == 1


# --------------------------------------------------------------------------
# FR-22 and FR-50/51: the masking policy and how much of the corpus it covers
# --------------------------------------------------------------------------
def test_fr22_policy_scope_is_stated_and_counts_come_from_the_frame():
    df = _frame({"HLA-A*01:01": (30, 30), "HLA-A*02:01": (25, 25), "HLA-B*07:02": (4, 4)})
    pol = build_bias_report(df, min_records_per_allele=20)["masking_policy"]
    assert pol["peptide_length"] == 9
    assert pol["n_alleles_in_scope"] == 2            # B*07:02 has only 8 records
    assert pol["n_alleles_total"] == 3
    assert pol["n_records_in_scope"] == 110 and pol["n_records_total"] == 118
    assert pol["tau_selected_on"] == "full_corpus_no_validation_split"
    n = pol["n_alleles_in_scope"]
    for key in (
        "n_alleles_nonempty_mask",
        "n_alleles_mask_equals_uniform",
        "n_alleles_mask_includes_position_2",
        "n_alleles_mask_includes_c_terminus",
    ):
        assert 0 <= pol[key] <= n, key
    assert "n_alleles_mask_differs_from_uniform" not in pol


def test_fr22_no_alleles_in_policy_scope_is_reported_unmeasured_not_zero():
    df = _frame({"HLA-A*01:01": (4, 4)})
    pol = build_bias_report(df, min_records_per_allele=20)["masking_policy"]
    assert pol["n_alleles_in_scope"] == 0
    assert pol["n_alleles_nonempty_mask"] == UNMEASURED
    assert pol["n_alleles_mask_equals_uniform"] == UNMEASURED
    assert pol["tau_nats"] == UNMEASURED


# --------------------------------------------------------------------------
# FR-60/61: admissibility for balanced splits, computed but never frozen here
# --------------------------------------------------------------------------
def test_fr60_split_admissibility_is_reported_and_is_not_a_frozen_split():
    # R1 bounds the ratio in the TEST set, which build_split meets by drawing a
    # balanced block, so a heavily imbalanced allele can still be admitted. An allele
    # is reported separately when it is too small or has only one class.
    df = _frame(
        {
            "HLA-A*01:01": (40, 40),    # admitted
            "HLA-A*02:01": (60, 3),     # admitted: balanced test block can be drawn
            "HLA-B*07:02": (10, 5),     # 15 records < min_samples_per_allele
            "HLA-B*08:01": (45, 0),     # single class: no ratio can be bounded
        }
    )
    sa = build_bias_report(df)["split_admissibility"]
    assert sa["frozen"] is False
    assert sa["n_admitted_alleles"] + sa["n_reported_separately"] == 4
    assert set(sa["reported_separately"]) == {"HLA-B*07:02", "HLA-B*08:01"}
    assert set(sa["admitted_alleles"]) == {"HLA-A*01:01", "HLA-A*02:01"}


# --------------------------------------------------------------------------
# Degenerate and malformed input must raise
# --------------------------------------------------------------------------
def test_degenerate_empty_frame_raises():
    with pytest.raises(BiasReportError, match="empty"):
        build_bias_report(_frame({}))


def test_degenerate_single_class_corpus_raises():
    df = _frame({"HLA-A*01:01": (30, 0), "HLA-A*02:01": (20, 0)})
    with pytest.raises(BiasReportError, match="both classes"):
        build_bias_report(df)


def test_missing_required_column_raises():
    df = _frame({"HLA-A*01:01": (30, 30)}).drop(columns=["assay_type"])
    with pytest.raises(BiasReportError, match="assay_type"):
        build_bias_report(df)


def test_non_binary_label_raises():
    df = _frame({"HLA-A*01:01": (30, 30)})
    df.loc[0, "label"] = 2
    with pytest.raises(BiasReportError, match="binary"):
        build_bias_report(df)


# --------------------------------------------------------------------------
# Artefact: deterministic, provenance-bearing, never written into frozen splits
# --------------------------------------------------------------------------
def test_report_records_input_identity_and_is_deterministic(tmp_path):
    df = _frame({"HLA-A*01:01": (30, 30), "HLA-A*02:01": (25, 25)})
    r1 = build_bias_report(df, input_sha256="deadbeef")
    r2 = build_bias_report(df, input_sha256="deadbeef")
    assert r1 == r2
    assert r1["input_file_sha256"] == "deadbeef"
    p1 = write_bias_report(r1, tmp_path / "bias.json")
    first = p1.read_bytes()
    p2 = write_bias_report(r2, tmp_path / "bias.json")
    assert p2.read_bytes() == first


def test_writer_refuses_to_write_under_frozen_splits(tmp_path):
    df = _frame({"HLA-A*01:01": (30, 30)})
    rep = build_bias_report(df)
    with pytest.raises(BiasReportError, match="splits"):
        write_bias_report(rep, tmp_path / "artefacts" / "splits" / "bias.json")


def test_report_numbers_come_from_the_frame_not_constants():
    a = build_bias_report(_frame({"HLA-A*01:01": (30, 30), "HLA-A*02:01": (25, 25)}))
    b = build_bias_report(_frame({"HLA-A*01:01": (30, 30), "HLA-A*02:01": (25, 25),
                                   "HLA-B*07:02": (20, 20)}))
    assert a["n_records"] != b["n_records"]
    assert a["n_alleles"] == 2 and b["n_alleles"] == 3
    assert math.isfinite(a["label_entropy_nats"])


# --------------------------------------------------------------------------
# The command: nbx bias profile
# --------------------------------------------------------------------------
def test_cli_profile_writes_report_and_records_input_hash(tmp_path):
    from typer.testing import CliRunner

    from nbx.cli import app

    df = _frame({"HLA-A*01:01": (30, 30), "HLA-A*02:01": (25, 25)})
    src = tmp_path / "corpus.parquet"
    df.to_parquet(src)
    out = tmp_path / "bias.json"
    result = CliRunner().invoke(app, ["bias", "profile", str(src), "--report", str(out)])
    assert result.exit_code == 0, result.output
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["n_records"] == len(df)
    assert len(data["input_file_sha256"]) == 64


def test_cli_profile_refuses_missing_dataset_with_nonzero_exit(tmp_path):
    from typer.testing import CliRunner

    from nbx.cli import app

    result = CliRunner().invoke(
        app, ["bias", "profile", str(tmp_path / "nope.parquet"),
              "--report", str(tmp_path / "b.json")]
    )
    assert result.exit_code == 1
    assert not (tmp_path / "b.json").exists()


def test_cli_profile_refuses_degenerate_corpus_and_writes_nothing(tmp_path):
    from typer.testing import CliRunner

    from nbx.cli import app

    src = tmp_path / "one_class.parquet"
    _frame({"HLA-A*01:01": (30, 0), "HLA-A*02:01": (20, 0)}).to_parquet(src)
    out = tmp_path / "b.json"
    result = CliRunner().invoke(app, ["bias", "profile", str(src), "--report", str(out)])
    assert result.exit_code == 1
    assert not out.exists()
