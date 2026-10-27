"""Sequence design: codon optimisation and repeat removal.

>>> from pysnapgene.design import CodonUsage, codon_optimize, longest_repeat
>>> usage = CodonUsage.from_csv("ecoli_codon_usage.csv")
>>> first = codon_optimize("GGGSGGSDGG", usage, forbidden=["GGTCTC"])
>>> second = codon_optimize("GGGSGGSDGG", usage, avoid=[first.dna], forbidden=["GGTCTC"])
>>> second.shared        # longest substring the two spellings share
6
"""

from __future__ import annotations

from .build import (
    Design,
    Part,
    build_dna,
    fix_complexity,
    frozen,
    identity,
    identity_floor,
    protein,
)
from .codons import CodonChoice, CodonUsage, codon_optimize, site_free, translate
from .complexity import (
    TWIST_GENE_FRAGMENT,
    Homology,
    Issue,
    Rules,
    complexity_report,
    duplex_tm,
    exact_repeats,
    gc_fraction,
    homologous_regions,
    homopolymers,
    repeat_density,
    score,
    summarise,
)
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
    "Part",
    "protein",
    "frozen",
    "Design",
    "build_dna",
    "fix_complexity",
    "identity",
    "identity_floor",
    "Rules",
    "TWIST_GENE_FRAGMENT",
    "Issue",
    "Homology",
    "complexity_report",
    "summarise",
    "score",
    "repeat_density",
    "exact_repeats",
    "homologous_regions",
    "duplex_tm",
    "gc_fraction",
    "homopolymers",
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
