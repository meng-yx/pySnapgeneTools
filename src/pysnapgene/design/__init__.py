"""Sequence design: codon optimisation and repeat removal.

>>> from pysnapgene.design import CodonUsage, codon_optimize, longest_repeat
>>> usage = CodonUsage.from_csv("ecoli_codon_usage.csv")
>>> first = codon_optimize("GGGSGGSDGG", usage, forbidden=["GGTCTC"])
>>> second = codon_optimize("GGGSGGSDGG", usage, avoid=[first.dna], forbidden=["GGTCTC"])
>>> second.shared        # longest substring the two spellings share
6
"""

from __future__ import annotations

from .codons import CodonChoice, CodonUsage, codon_optimize, site_free, translate
from .repeats import (
    DEFAULT_REPEAT_LIMIT,
    Repeat,
    RepeatEdits,
    break_repeats,
    longest_common_substring,
    longest_repeat,
    longest_shared_substring,
)

__all__ = [
    "CodonUsage",
    "CodonChoice",
    "codon_optimize",
    "translate",
    "site_free",
    "Repeat",
    "RepeatEdits",
    "DEFAULT_REPEAT_LIMIT",
    "longest_repeat",
    "longest_common_substring",
    "longest_shared_substring",
    "break_repeats",
]
