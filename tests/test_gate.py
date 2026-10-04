"""Truth table for the acceptance rule.

This file is the executable form of the sentence in Section 7.7:

    A result is not considered positive unless it improves discrimination
    without increasing shortcut dependence.

If a change to the codebase makes any of these pass when it should fail, the
project has lost the property it exists to demonstrate.
"""

from __future__ import annotations

from nbx.evaluate.gate import (
    UNMEASURED,
    ControlMargin,
    GateStatus,
    evaluate_gate,
    render_metric,
)

SIGNIFICANT = ControlMargin(delta=0.05, p_value=0.001)
NOT_SIGNIFICANT = ControlMargin(delta=0.004, p_value=0.42)
NEGATIVE_BUT_SIGNIFICANT = ControlMargin(delta=-0.05, p_value=0.001)
UNRUN = ControlMargin()


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------
def test_unmeasured_renders_as_em_dash_not_zero():
    assert render_metric(None) == UNMEASURED
    assert render_metric(None) != "0"
    assert render_metric(None) != ""


def test_zero_is_a_measurement_and_renders_as_zero():
    assert render_metric(0.0) == "0.0000"


def test_nan_and_inf_render_as_unmeasured():
    assert render_metric(float("nan")) == UNMEASURED
    assert render_metric(float("inf")) == UNMEASURED


# --------------------------------------------------------------------------
# The four states
# --------------------------------------------------------------------------
def test_missing_null_control_is_pending_and_hides_metrics():
    result = evaluate_gate(
        headline_metric=0.91,
        hla_only_margin=UNRUN,
        peptide_only_margin=SIGNIFICANT,
        kappa=0.02,
    )
    assert result.status is GateStatus.PENDING
    assert not result.rankable
    assert not result.displayable_metrics, (
        "FR-63: a headline metric must not be displayed without its controls"
    )
    assert "hla_only" in result.reason


def test_missing_peptide_only_control_is_also_pending():
    result = evaluate_gate(
        headline_metric=0.91,
        hla_only_margin=SIGNIFICANT,
        peptide_only_margin=UNRUN,
        kappa=0.02,
    )
    assert result.status is GateStatus.PENDING
    assert "peptide_only" in result.reason


def test_controls_present_but_margin_not_significant():
    result = evaluate_gate(
        headline_metric=0.91,
        hla_only_margin=NOT_SIGNIFICANT,
        peptide_only_margin=SIGNIFICANT,
        kappa=0.02,
        baseline_kappa=0.30,
    )
    assert result.status is GateStatus.NOT_SIGNIFICANT
    assert not result.rankable
    assert result.displayable_metrics, "metrics are shown; only the rank is withheld"


def test_a_significant_but_negative_margin_is_not_an_improvement():
    """Beating the control means beating it, not merely differing from it."""
    result = evaluate_gate(
        headline_metric=0.50,
        hla_only_margin=NEGATIVE_BUT_SIGNIFICANT,
        peptide_only_margin=SIGNIFICANT,
        kappa=0.02,
        baseline_kappa=0.30,
    )
    assert result.status is GateStatus.NOT_SIGNIFICANT
    assert not result.rankable


def test_margin_significant_but_shortcut_dependence_rose():
    result = evaluate_gate(
        headline_metric=0.93,
        hla_only_margin=SIGNIFICANT,
        peptide_only_margin=SIGNIFICANT,
        kappa=0.45,
        baseline_kappa=0.30,
    )
    assert result.status is GateStatus.FAILED_SHORTCUT
    assert not result.rankable
    assert "not surviving" in result.reason


def test_shortcut_comparison_uses_magnitude_not_sign():
    """kappa = -0.6 is as dependent as +0.6, in the opposite direction."""
    result = evaluate_gate(
        headline_metric=0.93,
        hla_only_margin=SIGNIFICANT,
        peptide_only_margin=SIGNIFICANT,
        kappa=-0.60,
        baseline_kappa=0.30,
    )
    assert result.status is GateStatus.FAILED_SHORTCUT


def test_both_conditions_met_passes_and_is_rankable():
    result = evaluate_gate(
        headline_metric=0.93,
        hla_only_margin=SIGNIFICANT,
        peptide_only_margin=SIGNIFICANT,
        kappa=0.05,
        baseline_kappa=0.30,
    )
    assert result.status is GateStatus.PASSED
    assert result.rankable
    assert result.displayable_metrics


def test_unmeasured_kappa_blocks_a_pass():
    result = evaluate_gate(
        headline_metric=0.93,
        hla_only_margin=SIGNIFICANT,
        peptide_only_margin=SIGNIFICANT,
        kappa=None,
        baseline_kappa=0.30,
    )
    assert result.status is GateStatus.PENDING
    assert not result.rankable


def test_null_controls_are_labelled_and_never_ranked():
    result = evaluate_gate(is_control=True)
    assert result.status is GateStatus.CONTROL
    assert not result.rankable


def test_gate_has_no_override_parameter():
    """There is deliberately no way to force a pass."""
    import inspect

    params = set(inspect.signature(evaluate_gate).parameters)
    for forbidden in ("force", "override", "skip_controls", "ignore_kappa"):
        assert forbidden not in params


# --------------------------------------------------------------------------
# ControlMargin semantics
# --------------------------------------------------------------------------
def test_margin_is_not_significant_until_both_parts_are_measured():
    assert not ControlMargin(delta=0.05).measured
    assert not ControlMargin(p_value=0.01).measured
    assert not ControlMargin(delta=0.05).significant


def test_unmeasured_margin_renders_as_em_dash():
    assert UNRUN.render() == UNMEASURED


def test_alpha_is_respected():
    borderline = ControlMargin(delta=0.05, p_value=0.04)
    assert borderline.significant
    assert not ControlMargin(delta=0.05, p_value=0.04, alpha=0.01).significant
