"""Module 7: Release and Reproducibility Packaging (FR-80 to FR-84).

NOT YET IMPLEMENTED. Contract only.

Packages deliverable D1 -- the bias-controlled benchmark protocol -- as
executable code rather than as a document.

Required behaviour
------------------
FR-80  Package configurations, environment specification, split hashes, seeds
       and result tables for public release.
FR-81  Generate the report figures and tables from the SAME artefacts the
       results came from, so a figure cannot drift from the number it depicts.
FR-82  A clean-environment run from the published repository reproduces the
       headline result tables from raw public data without manual intervention.
       This is the acceptance test, and it is demanding on purpose.
FR-83  Open licence; third-party baselines used strictly under their own terms.
FR-84  A preprint documenting methodology and results, INCLUDING negative
       results.

The protocol ships before any model result exists. D1 is the principal
contribution; a model result is a consumer of it.
"""

from __future__ import annotations

from pathlib import Path


def package_release(artefact_dir: Path, out: Path):
    """Assemble the public release bundle."""
    raise NotImplementedError(
        "Module 7 is not implemented. See PROMPTS.md."
    )
