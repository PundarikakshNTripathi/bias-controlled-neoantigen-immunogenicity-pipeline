"""Module 4: Model Core (FR-40 to FR-44).

NOT YET IMPLEMENTED. Contract only.

Architecture fixed by Section 7.6.1 of the synopsis
---------------------------------------------------
Interaction-based transformer encoder, chosen over a dual-encoder design so
that cross-attention weights between peptide and allele tokens can be read as
position-wise interaction evidence, which the interpretability analysis needs.

    4 encoder layers, 4 attention heads, hidden dimension 128,
    pre-layer-normalisation, dropout 0.1, two-layer projection head,
    Adam at learning rate 5e-5, early stopping on validation AUPRC.

This configuration is deliberately small. It has been shown sufficient to
match specialised architectures on binding and presentation, which is what
establishes that any deficiency observed on immunogenicity is attributable to
data and bias rather than to model underfitting.

FR-43 is the one that matters most: the null-control variants, in which the
peptide input or the allele input is withheld entirely, are part of the model
core and not an afterthought. The evaluation harness refuses to report a
headline metric without them.

FR-44  Before ANY immunogenicity claim: reproduce published baseline
       performance on binding and presentation within a stated tolerance.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BackboneConfig:
    """The configuration fixed by Section 7.6.1. Change it only with a reason."""

    n_layers: int = 4
    n_heads: int = 4
    d_model: int = 128
    dropout: float = 0.1
    learning_rate: float = 5e-5
    pre_layer_norm: bool = True


class InputVariant:
    """Which inputs a variant is permitted to see."""

    FULL = "peptide+hla"
    HLA_ONLY = "hla_only"            # mandatory null control
    PEPTIDE_ONLY = "peptide_only"    # mandatory null control


def build_model(config: BackboneConfig, variant: str = InputVariant.FULL):
    """Construct the encoder for one input variant."""
    raise NotImplementedError(
        "Module 4 is not implemented. See PROMPTS.md section 'M4 - Model Core'."
    )
