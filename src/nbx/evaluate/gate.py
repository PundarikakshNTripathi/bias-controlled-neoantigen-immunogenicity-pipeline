"""The acceptance rule, implemented as code rather than as a convention.

Section 7.7 of the synopsis states the rule in one sentence:

    A result is not considered positive unless it improves discrimination
    WITHOUT increasing shortcut dependence.

Everything in this module exists so that no interface, report or figure in
this project can display a result that has not passed that rule. The gate is
evaluated at the data-access layer, so a view added later cannot bypass it by
forgetting to check.

The burden of proof is on the model (decision D4). PASSED requires all of:

    (i)   a significant margin over the matched HLA-only control;
    (ii)  the upper 95% bound of [ |kappa_model| - |kappa_baseline| ] is at most
          NONINFERIORITY_MARGIN_KAPPA, so a model must be shown NOT to be more
          shortcut-dependent than the baseline, not merely not shown to be;
    (iv)  the model is not significantly worse than its peptide-only control.

Condition (iii), an absolute ceiling on |kappa_model|, is NOT implemented. Its
value is a scientific judgement that belongs to the researcher and has not been
given; there is no parameter and no default for it. Until it exists, a model
exactly as shortcut-dependent as its baseline, with a degenerate zero-width
interval, is not stopped by (ii) alone.

The five reachable states are enumerated in `GateStatus` and tested as a truth
table. There is deliberately no "force" or "override" argument, and the
non-inferiority margin is a module constant, not a per-call argument.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

__all__ = [
    "NONINFERIORITY_MARGIN_KAPPA",
    "GateStatus",
    "ControlMargin",
    "GateResult",
    "evaluate_gate",
    "UNMEASURED",
    "render_metric",
]

NONINFERIORITY_MARGIN_KAPPA: float = 0.0
"""Pre-registered default (analysis plan, decision D4): the largest upper 95% bound
of the change in |kappa| relative to the baseline that still counts as "not more
shortcut-dependent". Zero means the interval must lie at or below no increase.
Changing it changes the rule and needs the researcher's approval."""

# What an unmeasured quantity renders as, everywhere, without exception.
UNMEASURED = "—"  # em dash


def render_metric(value: float | None, digits: int = 4) -> str:
    """Render a metric for display.

    A quantity that has not been computed renders as an em dash. It is never
    zero, never an empty string, and never a plausible-looking stand-in.
    Zero is a measurement; absence is not.
    """
    if value is None:
        return UNMEASURED
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return UNMEASURED
    return f"{value:.{digits}f}"


class GateStatus(str, Enum):
    """The states a model variant can be in.

    Condition (iv), significantly worse than the peptide-only control, is
    reported under FAILED_SHORTCUT together with condition (ii). That is a
    design choice for the researcher to ratify; the reason string names which
    condition failed.
    """

    PENDING = "pending"
    """A mandated null control is missing, or a required metric is unmeasured."""

    NOT_SIGNIFICANT = "not_significant"
    """Controls present; the margin over the HLA-only control is not significant."""

    FAILED_SHORTCUT = "failed_shortcut"
    """Margin is significant, but condition (ii) (shortcut dependence not shown to
    be no higher than the baseline's) or condition (iv) (significantly worse than
    the peptide-only control) failed. Reported, not ranked. DESIGN CHOICE for the
    researcher to ratify: (iv) is reported under this status rather than under a
    sixth one, and the reason string names which condition failed."""

    PASSED = "passed"
    """Discrimination improved and shortcut dependence was shown not to rise."""

    CONTROL = "control"
    """This row is itself a null control; it is never ranked."""


@dataclass(frozen=True)
class ControlMargin:
    """A model's margin over a null control, with its significance verdict.

    `p_value` is whatever the evaluation harness computed (a paired test
    across seeds, by default). `delta` is the point difference in the headline
    metric. Both may be None before the controls have run.
    """

    delta: float | None = None
    p_value: float | None = None
    alpha: float = 0.05

    @property
    def measured(self) -> bool:
        return self.delta is not None and self.p_value is not None

    @property
    def significant(self) -> bool:
        """Significant means measured, positive and below alpha. All three."""
        delta, p_value = self.delta, self.p_value
        if delta is None or p_value is None:
            return False
        return delta > 0.0 and p_value < self.alpha

    def render(self) -> str:
        if not self.measured:
            return UNMEASURED
        mark = "significant" if self.significant else "n.s."
        return f"{self.delta:+.4f} ({mark})"


@dataclass(frozen=True)
class GateResult:
    """The outcome of applying the acceptance rule to one model variant."""

    status: GateStatus
    reason: str
    rankable: bool

    @property
    def displayable_metrics(self) -> bool:
        """Whether numeric metrics may be shown at all for this variant.

        A variant with no null control may not display a headline metric:
        requirement FR-63. It may still display its configuration.
        """
        return self.status is not GateStatus.PENDING


def evaluate_gate(
    *,
    is_control: bool = False,
    headline_metric: float | None = None,
    hla_only_margin: ControlMargin | None = None,
    peptide_only_margin: ControlMargin | None = None,
    kappa_change_upper_bound: float | None = None,
    peptide_only_significantly_worse: bool | None = None,
) -> GateResult:
    """Apply the acceptance rule to one model variant.

    Parameters
    ----------
    is_control
        True when this row is itself a null control. Controls are reported
        inline and labelled, and are never ranked.
    headline_metric
        The discrimination metric (AUPRC by default). None means unmeasured.
    hla_only_margin
        Margin over the mandatory HLA-only control. Required. Condition (i).
    peptide_only_margin
        Margin over the peptide-only control. Required to be present as proof
        that the control was run; whether the model is significantly WORSE is
        carried by `peptide_only_significantly_worse`.
    kappa_change_upper_bound
        Upper 95% bound of |kappa_model| - |kappa_baseline|, computed by the
        evaluation harness (an allele-clustered bootstrap). Condition (ii) is
        `kappa_change_upper_bound <= NONINFERIORITY_MARGIN_KAPPA`. None means
        the interval was not measured, which blocks a pass.
    peptide_only_significantly_worse
        True when the harness found the model significantly worse than the
        peptide-only control (one-sided, alpha 0.05, opposite tail). Condition
        (iv). None means the comparison was not made, which blocks a pass.

    Condition (iii), an absolute ceiling on |kappa_model|, is not implemented;
    see the module docstring.

    Returns
    -------
    GateResult
        Status, a human-readable reason, and whether a rank may be assigned.
    """
    if is_control:
        return GateResult(
            status=GateStatus.CONTROL,
            reason="Null control. Reported inline and never ranked.",
            rankable=False,
        )

    # --- R2: both mandated null controls must have been run ---------------
    missing = []
    if hla_only_margin is None or not hla_only_margin.measured:
        missing.append("hla_only")
    if peptide_only_margin is None or not peptide_only_margin.measured:
        missing.append("peptide_only")
    if missing or hla_only_margin is None:
        return GateResult(
            status=GateStatus.PENDING,
            reason=(
                "Mandated null control(s) missing: "
                + ", ".join(missing)
                + ". No headline metric may be displayed (FR-63)."
            ),
            rankable=False,
        )

    if headline_metric is None:
        return GateResult(
            status=GateStatus.PENDING,
            reason="Headline metric has not been measured.",
            rankable=False,
        )

    # --- (i) the margin must be significant -------------------------------
    hla_margin: ControlMargin = hla_only_margin
    if not hla_margin.significant:
        return GateResult(
            status=GateStatus.NOT_SIGNIFICANT,
            reason=(
                "Does not exceed the HLA-only control by a statistically "
                f"significant margin ({hla_margin.render()}). "
                "Metrics are shown; no rank is assigned."
            ),
            rankable=False,
        )

    # --- the shortcut conditions need their inputs ------------------------
    if kappa_change_upper_bound is None:
        return GateResult(
            status=GateStatus.PENDING,
            reason=(
                "The interval for the change in shortcut dependence kappa has "
                "not been measured."
            ),
            rankable=False,
        )
    if peptide_only_significantly_worse is None:
        return GateResult(
            status=GateStatus.PENDING,
            reason="The comparison with the peptide-only control has not been made.",
            rankable=False,
        )

    failed = []
    if kappa_change_upper_bound > NONINFERIORITY_MARGIN_KAPPA:
        failed.append(
            "condition (ii): the upper 95% bound of the change in |kappa| over "
            f"the baseline is {kappa_change_upper_bound:+.4f}, above the "
            f"non-inferiority margin {NONINFERIORITY_MARGIN_KAPPA:+.4f}; shortcut "
            "dependence is not shown to be no higher than the baseline's"
        )
    if peptide_only_significantly_worse:
        failed.append(
            "condition (iv): the model is significantly worse than the "
            "peptide-only control"
        )
    if failed:
        return GateResult(
            status=GateStatus.FAILED_SHORTCUT,
            reason=(
                f"Discrimination improved ({hla_margin.render()}) but "
                + "; and ".join(failed)
                + ". Reported as not surviving the controls."
            ),
            rankable=False,
        )

    return GateResult(
        status=GateStatus.PASSED,
        reason=(
            f"Exceeds the HLA-only control ({hla_margin.render()}); the upper "
            "95% bound of the change in |kappa| over the baseline is "
            f"{kappa_change_upper_bound:+.4f} (margin "
            f"{NONINFERIORITY_MARGIN_KAPPA:+.4f}); not significantly worse than "
            "the peptide-only control."
        ),
        rankable=True,
    )
