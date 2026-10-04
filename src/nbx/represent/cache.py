"""Module 3: Representation Service (FR-30 to FR-35).

NOT YET IMPLEMENTED. Contract only.

Representations are invariant across the ablation matrix, so computing them
once and caching converts the dominant cost from a per-run cost into a
one-time cost. That is the whole reason this module exists separately.

Required behaviour
------------------
FR-30  Tokenised peptide and pseudo-sequence tensors, with positional
       encodings applied separately to the two segments so that positional
       semantics are not conflated across the peptide-allele boundary.
FR-31  Optional fixed protein-language-model embeddings (feature extractor,
       never fine-tuned - that is outside the compute budget).
FR-32  Optional geometric interface graph from an open AlphaFold3-class
       structure predictor.
FR-33  Per-residue predicted-structure confidence carried through as an
       EXPLICIT feature. A low-confidence prediction must never be silently
       treated as equivalent to a high-confidence one.
FR-34  Deterministic given a fixed input and cache key.
FR-35  Where structural coverage is budget-restricted, report coverage
       explicitly and evaluate the affected hypothesis on the covered subset.
"""

from __future__ import annotations

from pathlib import Path


def build_cache(dataset_path: Path, cache_dir: Path, with_geometry: bool = False):
    """Compute and cache the fixed representation channels."""
    raise NotImplementedError(
        "Module 3 is not implemented. See PROMPTS.md section "
        "'M3 - Representation Service'."
    )
