"""Module 6: Bias-Controlled Evaluation Harness (FR-60 to FR-69).

PARTIALLY IMPLEMENTED. The acceptance rule (gate.py), the metrics
(metrics.py) and the split construction (nbx/splits/freeze.py) are complete
and tested. What remains is the orchestration that runs a model variant and
its mandated controls across seeds and assembles the report.

Required behaviour still to build
---------------------------------
FR-62  Execute an HLA-only and a peptide-only control of MATCHED capacity and
       training budget for every reported model. Matched means matched: equal
       epochs, equal tuning effort, equal early-stopping criterion.
FR-65  Distribution-shifted control sets, in which the intra-allelic imbalance
       of the test data is inverted relative to training.
FR-68  Report means with across-seed variance, never best-run figures.

The acceptance test for this module, stated in Section 9 of the SRS, is the
demanding one: run the harness on a PUBLISHED model and a PUBLISHED benchmark,
reproduce that model's published figure under pooled evaluation, and
additionally report the allele-stratified and control-margin figures the
original evaluation omitted.
"""

from __future__ import annotations

from pathlib import Path


def run_protocol(run_dir: Path, split_hash: str):
    """Execute the full bias-controlled evaluation for one ablation matrix."""
    raise NotImplementedError(
        "Protocol orchestration is not implemented. See PROMPTS.md section "
        "'M6 - Evaluation Harness'."
    )
