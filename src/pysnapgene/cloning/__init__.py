"""Cloning simulation: type IIS digestion, Golden Gate, and homology-based assembly.

>>> from Bio.Restriction import BsaI
>>> from pysnapgene.cloning import golden_gate, translate_orf
>>> product = golden_gate(vector_seq, insert_seq, BsaI)
>>> translate_orf(product, atg_position)
'MVTGYRLFEEILGGG...'

>>> from pysnapgene.cloning import infusion
>>> assembly = infusion([linearised_backbone, synthesised_fragment])
>>> len(assembly), [len(j) for j in assembly.junctions]
(5859, [15, 15])
"""

from __future__ import annotations

from .goldengate import Fragment, digest, excise, golden_gate, translate_orf
from .homology import (
    INFUSION_MIN_HOMOLOGY,
    Assembly,
    Junction,
    ambiguous_homology,
    homology_assemble,
    infusion,
)

__all__ = [
    "Fragment",
    "digest",
    "excise",
    "golden_gate",
    "translate_orf",
    "Assembly",
    "Junction",
    "homology_assemble",
    "infusion",
    "ambiguous_homology",
    "INFUSION_MIN_HOMOLOGY",
]
