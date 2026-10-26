"""Codon optimisation against a usage table.

The point of the ``avoid`` argument is the one that matters for plasmid design: when an
element appears twice in the same construct, the two copies must not be spelled the same,
or a single recombination across the repeat deletes everything between them. So a spelling
is chosen to share as short a substring as possible with the copies already placed.
"""

from __future__ import annotations

import random
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

from Bio.Seq import Seq

__all__ = ["CodonUsage", "CodonChoice", "codon_optimize", "translate", "site_free"]


def translate(nucleotides: str) -> str:
    """Protein of ``nucleotides``, trimmed to whole codons.

    Trimming rather than raising is deliberate: design code slices sub-sequences out of a
    construct constantly, and a trailing partial codon is normal rather than an error.
    """
    usable = nucleotides[: len(nucleotides) // 3 * 3]
    return str(Seq(usable).translate()) if usable else ""


def site_free(sequence: str, motifs: Iterable[str], *, both_strands: bool = True) -> bool:
    """True when none of ``motifs`` occurs in ``sequence``.

    With ``both_strands`` the reverse complement of each motif is checked too, which is
    what you want for a restriction site: an enzyme does not care which strand it reads.
    """
    sequence = sequence.upper()
    for motif in motifs:
        motif = motif.upper()
        if motif in sequence:
            return False
        if both_strands and str(Seq(motif).reverse_complement()) in sequence:
            return False
    return True


class CodonUsage:
    """Per-organism codon usage: which codons encode an amino acid, and how often.

    Built from a table of ``codon, amino_acid, fraction`` rows, where *fraction* is the
    share of that codon among the synonyms for its amino acid.
    """

    def __init__(self, fractions: Mapping[str, tuple[str, float]]):
        self._fraction: dict[str, float] = {}
        self._by_aa: dict[str, list[tuple[str, float]]] = {}
        for codon, (amino_acid, fraction) in fractions.items():
            codon = codon.upper().replace("U", "T")
            self._fraction[codon] = float(fraction)
            self._by_aa.setdefault(amino_acid, []).append((codon, float(fraction)))
        for codons in self._by_aa.values():
            codons.sort(key=lambda pair: pair[1], reverse=True)
        if not self._by_aa:
            raise ValueError("codon usage table is empty")

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        codon_column: str = "codon",
        amino_acid_column: str = "amino_acid",
        fraction_column: str = "fraction",
    ) -> CodonUsage:
        """Read a usage table from CSV, without requiring pandas."""
        import csv

        fractions: dict[str, tuple[str, float]] = {}
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                fractions[row[codon_column]] = (
                    row[amino_acid_column],
                    float(row[fraction_column]),
                )
        return cls(fractions)

    def fraction(self, codon: str) -> float:
        return self._fraction[codon.upper().replace("U", "T")]

    def amino_acids(self) -> list[str]:
        return sorted(self._by_aa)

    def codons(self, amino_acid: str, *, min_fraction: float = 0.0) -> list[str]:
        """Codons for ``amino_acid``, commonest first, at or above ``min_fraction``.

        Never returns empty for a known amino acid: if the threshold excludes everything,
        the commonest codon is still offered, because a protein has to be encodable.
        """
        try:
            codons = self._by_aa[amino_acid]
        except KeyError:
            raise KeyError(f"no codons for amino acid {amino_acid!r}") from None
        kept = [codon for codon, fraction in codons if fraction >= min_fraction]
        return kept or [codons[0][0]]

    def best(self, amino_acid: str) -> str:
        """The single commonest codon for ``amino_acid``."""
        return self._by_aa[amino_acid][0][0]

    def synonyms(self, codon: str, *, min_fraction: float = 0.0) -> list[str]:
        """Other codons for the same amino acid, commonest first."""
        codon = codon.upper().replace("U", "T")
        amino_acid = str(Seq(codon).translate())
        return [c for c in self.codons(amino_acid, min_fraction=min_fraction) if c != codon]


class CodonChoice(NamedTuple):
    """A chosen spelling, and the worst substring it shares with the sequences avoided."""

    dna: str
    shared: int


def _longest_common(a: str, b: str) -> int:
    from .repeats import longest_common_substring

    return longest_common_substring(a, b)


def codon_optimize(
    protein: str,
    usage: CodonUsage,
    *,
    rng: random.Random | None = None,
    avoid: Sequence[str] = (),
    forbidden: Sequence[str] = (),
    min_fraction: float = 0.10,
    tries: int = 500,
    target: int = 10,
    argmax: bool = False,
) -> CodonChoice:
    """Encode ``protein`` as DNA.

    ``argmax`` takes the commonest codon everywhere, which is deterministic but gives one
    spelling per protein -- no use when the same element appears twice. Otherwise codons are
    sampled in proportion to their usage, and the best of ``tries`` attempts is kept:
    *best* meaning the shortest substring shared with anything in ``avoid``. Sampling stops
    early once a spelling shares no more than ``target`` bases.

    ``forbidden`` motifs (restriction sites, typically) are rejected on both strands.
    """
    rng = rng or random.Random()
    if argmax:
        dna = "".join(usage.best(amino_acid) for amino_acid in protein)
        if forbidden and not site_free(dna, forbidden):
            raise ValueError("the argmax spelling carries a forbidden motif")
        return CodonChoice(dna, max((_longest_common(dna, a) for a in avoid), default=0))

    best: CodonChoice | None = None
    for _ in range(tries):
        codons = []
        for amino_acid in protein:
            choices = usage.codons(amino_acid, min_fraction=min_fraction)
            weights = [usage.fraction(codon) for codon in choices]
            codons.append(rng.choices(choices, weights=weights)[0])
        dna = "".join(codons)
        if forbidden and not site_free(dna, forbidden):
            continue
        shared = max((_longest_common(dna, a) for a in avoid), default=0)
        if best is None or shared < best.shared:
            best = CodonChoice(dna, shared)
        if best.shared <= target:
            break
    if best is None:
        raise ValueError(
            f"no spelling of {protein[:20]!r} avoided the forbidden motifs in {tries} tries"
        )
    return best
