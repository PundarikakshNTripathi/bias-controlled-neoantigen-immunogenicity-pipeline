"""Truth table for the acceptance rule.

This file is the executable form of the sentence in Section 7.7:

    A result is not considered positive unless it improves discrimination
    without increasing shortcut dependence.

The rule puts the burden of proof on the MODEL (decision D4). PASSED needs
(i) a significant margin over the HLA-only control, (ii) the upper 95% bound of
[|kappa_model| - |kappa_baseline|] at most NONINFERIORITY_MARGIN_KAPPA, and
(iv) the model not significantly worse than the peptide-only control. Condition
(iii), an absolute ceiling on |kappa_model|, is deliberately not implemented: its
value is the researcher's decision and has not been given.

Every number below is illustrative: it pins the logic, it is not a measurement.

If a change to the codebase makes any of these pass when it should fail, the
project has lost the property it exists to demonstrate.
"""

from __future__ import annotations

import inspect

from nbx.evaluate.gate import (
    NONINFERIORITY_MARGIN_KAPPA,
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
# Helpers: a fully specified, passing set of inputs that each test perturbs
# --------------------------------------------------------------------------
def gate(**overrides):
    inputs = dict(
        headline_metric=0.93,
        hla_only_margin=SIGNIFICANT,
        peptide_only_margin=NOT_SIGNIFICANT,
        kappa_change_upper_bound=-0.10,
        peptide_only_significantly_worse=False,
    )
    inputs.update(overrides)
    return evaluate_gate(**inputs)


def test_the_noninferiority_margin_is_the_preregistered_default():
    assert NONINFERIORITY_MARGIN_KAPPA == 0.0


# --------------------------------------------------------------------------
# The states
# --------------------------------------------------------------------------
def test_missing_hla_only_control_is_pending_and_hides_metrics():
    result = gate(hla_only_margin=UNRUN)
    assert result.status is GateStatus.PENDING
    assert not result.rankable
    assert not result.displayable_metrics, (
        "FR-63: a headline metric must not be displayed without its controls"
    )
    assert "hla_only" in result.reason


def test_missing_peptide_only_control_is_also_pending():
    result = gate(peptide_only_margin=UNRUN)
    assert result.status is GateStatus.PENDING
    assert "peptide_only" in result.reason


def test_unmeasured_headline_metric_is_pending():
    assert gate(headline_metric=None).status is GateStatus.PENDING


def test_controls_present_but_margin_not_significant():
    result = gate(hla_only_margin=NOT_SIGNIFICANT, kappa_change_upper_bound=0.5)
    assert result.status is GateStatus.NOT_SIGNIFICANT
    assert not result.rankable
    assert result.displayable_metrics, "metrics are shown; only the rank is withheld"


def test_a_significant_but_negative_margin_is_not_an_improvement():
    """Beating the control means beating it, not merely differing from it."""
    result = gate(hla_only_margin=NEGATIVE_BUT_SIGNIFICANT)
    assert result.status is GateStatus.NOT_SIGNIFICANT
    assert not result.rankable


def test_significant_margin_and_kappa_change_bound_above_the_margin_fails():
    result = gate(kappa_change_upper_bound=NONINFERIORITY_MARGIN_KAPPA + 0.03)
    assert result.status is GateStatus.FAILED_SHORTCUT
    assert not result.rankable
    assert "(ii)" in result.reason and "not surviving" in result.reason
    assert "(iv)" not in result.reason


def test_significant_margin_kappa_bound_within_margin_and_peptide_not_worse_passes():
    result = gate(kappa_change_upper_bound=NONINFERIORITY_MARGIN_KAPPA - 0.20)
    assert result.status is GateStatus.PASSED
    assert result.rankable
    assert result.displayable_metrics


def test_a_bound_exactly_at_the_margin_is_within_it():
    result = gate(kappa_change_upper_bound=NONINFERIORITY_MARGIN_KAPPA)
    assert result.status is GateStatus.PASSED


def test_significant_margin_but_significantly_worse_than_peptide_only_fails():
    result = gate(peptide_only_significantly_worse=True)
    assert result.status is GateStatus.FAILED_SHORTCUT
    assert not result.rankable
    assert "(iv)" in result.reason and "(ii)" not in result.reason


def test_failure_of_both_conditions_names_both():
    result = gate(kappa_change_upper_bound=0.2, peptide_only_significantly_worse=True)
    assert result.status is GateStatus.FAILED_SHORTCUT
    assert "(ii)" in result.reason and "(iv)" in result.reason


def test_an_equally_shortcut_dependent_baseline_does_not_let_the_model_pass():
    """The audit's measured failure (illustrative numbers, not a measurement).

    Model and baseline both score alleles by their training label rate, so
    |kappa| is near its maximum for both and the change in |kappa| is centred on
    zero: its 95% interval is, say, [-0.02, +0.03]. A rule that only asked
    whether kappa rose, on the point estimate, could pass this model. Here the
    upper bound is positive, the burden of proof is on the model, and it fails.
    """
    result = gate(kappa_change_upper_bound=0.03)
    assert result.status is GateStatus.FAILED_SHORTCUT
    assert not result.rankable


def test_unmeasured_kappa_interval_blocks_a_pass():
    result = gate(kappa_change_upper_bound=None)
    assert result.status is GateStatus.PENDING
    assert not result.rankable
    assert "kappa" in result.reason


def test_unmeasured_peptide_only_comparison_blocks_a_pass():
    result = gate(peptide_only_significantly_worse=None)
    assert result.status is GateStatus.PENDING
    assert not result.rankable
    assert "peptide" in result.reason


def test_not_significant_is_reported_before_missing_shortcut_inputs():
    result = gate(
        hla_only_margin=NOT_SIGNIFICANT,
        kappa_change_upper_bound=None,
        peptide_only_significantly_worse=None,
    )
    assert result.status is GateStatus.NOT_SIGNIFICANT


def test_null_controls_are_labelled_and_never_ranked():
    result = evaluate_gate(is_control=True)
    assert result.status is GateStatus.CONTROL
    assert not result.rankable


def test_gate_has_no_override_or_tolerance_parameter():
    """There is deliberately no way to force a pass, to loosen the margin per call,
    or to compare against a point-estimate tolerance."""
    params = set(inspect.signature(evaluate_gate).parameters)
    for forbidden in ("force", "override", "skip_controls", "ignore_kappa"):
        assert forbidden not in params
    for removed in ("kappa_tolerance", "baseline_kappa", "kappa"):
        assert removed not in params
    assert not any("noninferiority" in name for name in params)


def test_no_absolute_kappa_ceiling_is_implemented_and_the_gap_is_documented():
    """Condition (iii) is the researcher's decision. It is a documented gap, not a
    parameter and not a default."""
    params = set(inspect.signature(evaluate_gate).parameters)
    assert not any("ceiling" in p or "max_abs" in p for p in params)
    module_doc = inspect.getmodule(evaluate_gate).__doc__ or ""
    doc = (evaluate_gate.__doc__ or "") + module_doc
    assert "(iii)" in doc and "not implemented" in doc.lower()


def test_peptide_only_failure_is_reported_under_failed_shortcut():
    """Design choice flagged for the researcher: the five statuses are unchanged,
    so condition (iv) shares FAILED_SHORTCUT with condition (ii)."""
    assert {s.name for s in GateStatus} == {
        "PENDING", "NOT_SIGNIFICANT", "FAILED_SHORTCUT", "PASSED", "CONTROL"
    }
    assert "(iv)" in (GateStatus.__doc__ or "")


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
