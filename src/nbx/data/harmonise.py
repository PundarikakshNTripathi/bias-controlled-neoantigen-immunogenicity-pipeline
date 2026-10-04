"""Module 1: Data Harmonisation Engine (FR-01 to FR-06).

NOT YET IMPLEMENTED. This file is the contract Claude Code implements against.

Harmonisation is a substantive step, not a clerical one: errors introduced
here propagate silently into every downstream measurement. The acceptance
test is in `tests/test_harmonise.py`, which must be written alongside the
implementation.

Required behaviour
------------------
FR-01  Ingest IEDB, TESLA-derived and HiTIDE corpora into one canonical frame.
FR-02  Reconcile allele nomenclature; assign 34-residue contact pseudo-sequences.
FR-03  Filter peptides to 8-14 residues; deduplicate across overlapping sources.
FR-04  Retain source dataset, assay type and negative-sample provenance per record.
FR-05  Emit a content-hashed, frozen canonical dataset.
FR-06  Reproduce published dataset statistics where a reference exists; document
       any discrepancy rather than quietly matching.

Hard rules
----------
* A record without provenance is dropped, never defaulted. Positive-to-negative
  ratios are artefacts of experimental design, so an analysis that pools across
  sources without provenance is uninterpretable.
* An allele lacking sequence information is either removed or mapped to a
  same-two-digit-typing representative, and the choice is recorded PER RECORD.
* Record counts reported anywhere must come from the frame, never from a
  constant in the source.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PEPTIDE_MIN_LENGTH = 8
PEPTIDE_MAX_LENGTH = 14
PSEUDO_SEQUENCE_LENGTH = 34   # NetMHCpan convention: polymorphic residues
                              # within 4.0 A of the peptide


@dataclass(frozen=True)
class Provenance:
    """Where a record came from. Mandatory on every record."""

    source_dataset: str
    assay_type: str
    negative_sample_origin: str
    allele_mapping_note: str = ""


def harmonise(sources: dict[str, Path], out_dir: Path):
    """Build the canonical dataset. See FR-01 to FR-06."""
    raise NotImplementedError(
        "Module 1 is not implemented. Run /nbx:module M1 in Claude Code, or see "
        "PROMPTS.md section 'M1 - Data Harmonisation Engine'."
    )
