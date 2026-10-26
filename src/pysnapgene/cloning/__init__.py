"""Cloning simulation: type IIS digestion and Golden Gate assembly.

>>> from Bio.Restriction import BsaI
>>> from pysnapgene.cloning import golden_gate, translate_orf
>>> product = golden_gate(vector_seq, insert_seq, BsaI)
>>> translate_orf(product, atg_position)
'MVTGYRLFEEILGGG...'
"""

from __future__ import annotations

from .goldengate import Fragment, digest, excise, golden_gate, translate_orf

__all__ = ["Fragment", "digest", "excise", "golden_gate", "translate_orf"]
