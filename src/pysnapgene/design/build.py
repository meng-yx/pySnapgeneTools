"""Protein to DNA, deterministically, in a form a synthesis vendor will accept.

One function turns a parts list into a sequence: parts that are proteins get codons, parts
that are already DNA are kept verbatim, and the result is then edited -- synonymously, one
codon at a time -- until it stops breaking the vendor's rules.

Three things make it usable for ordering rather than just for encoding:

*Nothing is random.* The same parts and the same usage table give the same DNA, on any
machine, in any order. Codon choice is a total order over the candidates -- fewest bases
shared with the copies you asked it to differ from, then fewest bases repeated from what is
already built, then commonest codon, then alphabetical -- so there is always exactly one
winner and no seed to record.

*Some parts cannot move.* A ccdB cassette, a cloning overhang, a promoter: pass them as
``frozen`` and no repair will touch them, which is the difference between a tool you can
use on a real construct and one you can use on a toy.

*Two copies of one protein get different spellings.* ``unlike`` tells a part to share as
little as possible with a part already built. That matters more than it sounds: two
independently optimised copies of a 157 aa protein are still ~79 % identical base for base,
which vendors read as one long repeat and cells read as a recombination substrate.

    >>> plan = [
    ...     frozen("overhang", "AATTGGTCTCAAGGA"),
    ...     protein("tag", "VTGYRLFEEIL"),
    ...     protein("tag_2", "VTGYRLFEEIL", unlike=["tag"]),
    ... ]
    >>> design = build_dna(plan, usage, forbidden=["GGTCTC"])   # doctest: +SKIP

What it does not do: it never changes the protein, never reorders parts, and never inserts
or deletes a base, so a construct whose complexity comes from its architecture -- two copies
of the same long protein, say -- can come back still flagged. The report says which rule is
still broken, because the remedy for that is a design decision, not a codon.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import NamedTuple

from Bio.Seq import Seq

from .codons import CodonUsage, site_free
from .complexity import Issue, Rules, TWIST_GENE_FRAGMENT, complexity_report, score

__all__ = [
    "Part",
    "protein",
    "frozen",
    "Design",
    "build_dna",
    "fix_complexity",
    "identity",
    "identity_floor",
]

#: how far back a candidate codon is checked for repeating what is already built
_REPEAT_LOOKBACK = 24


class Part(NamedTuple):
    """One piece of a construct: a protein to encode, or DNA to keep as it is."""

    name: str
    #: the protein to encode, or None for a DNA part
    peptide: str | None
    #: the DNA, for a part given as DNA
    dna: str | None
    #: names of earlier parts this one should share as little sequence as possible with
    unlike: tuple[str, ...] = ()
    #: whether repairs may change this part's codons
    editable: bool = True

    @property
    def is_protein(self) -> bool:
        return self.peptide is not None


def protein(name: str, peptide: str, *, unlike: Sequence[str] = ()) -> Part:
    """A part given as protein, to be encoded.

    ``unlike`` names earlier parts whose spelling this one should avoid -- use it for the
    second copy of anything that appears twice.
    """
    peptide = peptide.upper()
    if not peptide:
        raise ValueError(f"part {name!r} has an empty protein sequence")
    if "*" in peptide[:-1]:
        raise ValueError(f"part {name!r} has an internal stop codon")
    return Part(name, peptide, None, tuple(unlike), True)


def frozen(name: str, sequence: str, *, editable: bool = False) -> Part:
    """A part given as DNA, kept verbatim.

    This is how a cassette that must not change -- ccdB, an overhang, a promoter, an RBS --
    survives a repair pass. Pass ``editable=True`` for DNA that *may* be re-spelled, in
    which case it must be a whole number of codons in its own frame.
    """
    sequence = sequence.upper()
    if editable and len(sequence) % 3:
        raise ValueError(
            f"part {name!r} is editable but {len(sequence)} bp is not a whole number of "
            f"codons; an edit swaps a codon, so it needs a frame"
        )
    return Part(name, None, sequence, (), editable)


class Design(NamedTuple):
    """What :func:`build_dna` and :func:`fix_complexity` return."""

    sequence: str
    #: name -> half-open span in ``sequence``; empty for :func:`fix_complexity`
    spans: dict[str, tuple[int, int]]
    #: rules still broken, worst first; empty means the vendor's checks all pass
    issues: list[Issue]
    #: bases changed by the repair pass
    edits: int
    #: repair rounds used
    rounds: int

    @property
    def clean(self) -> bool:
        return not self.issues

    def part(self, name: str) -> str:
        start, end = self.spans[name]
        return self.sequence[start:end]

    def report(self) -> str:
        lines = [f"{len(self.sequence)} bp, {self.edits} base(s) changed in "
                 f"{self.rounds} repair round(s)"]
        if self.issues:
            lines.append(f"{len(self.issues)} rule(s) still broken:")
            lines += [f"  {issue}" for issue in self.issues]
        else:
            lines.append("no rule broken")
        return "\n".join(lines)


def _codon_identity(candidate: str, references: Sequence[str]) -> int:
    """Bases this codon shares with the codon at the same place in each reference."""
    return sum(
        1
        for reference in references
        for a, b in zip(candidate, reference)
        if a == b
    )


def _tail_repeat(sequence: str, seed: int) -> int:
    """Length of the longest suffix of ``sequence`` that also occurs earlier in it.

    Capped at :data:`_REPEAT_LOOKBACK`: the point is to rank candidate codons, and a
    repeat longer than that is already as bad as the ranking needs to know.
    """
    longest = min(_REPEAT_LOOKBACK, len(sequence) - 1)
    for size in range(longest, seed - 1, -1):
        if sequence.find(sequence[-size:], 0, len(sequence) - 1) >= 0:
            return size
    return 0


def _creates_forbidden(sequence: str, forbidden: Sequence[str], added: int) -> bool:
    """Whether appending the last ``added`` bases has produced a forbidden motif.

    Only placements that overlap the new bases count. Checking a fixed window either side
    instead would reject every codon placed next to a frozen part that contains the motif
    on purpose -- which is exactly where a ccdB cassette's own BsaI sites sit.
    """
    for motif in forbidden:
        width = len(motif) + added - 1
        if not site_free(sequence[-width:], [motif]):
            return True
    return False


def _sites(sequence: str, forbidden: Sequence[str]) -> int:
    """Occurrences of every forbidden motif, both strands."""
    total = 0
    for motif in forbidden:
        total += sequence.count(motif)
        total += sequence.count(str(Seq(motif).reverse_complement()))
    return total


def build_dna(
    parts: Iterable[Part],
    usage: CodonUsage,
    *,
    forbidden: Sequence[str] = (),
    rules: Rules = TWIST_GENE_FRAGMENT,
    min_fraction: float = 0.10,
    repair: bool = True,
    rounds: int = 400,
    ignore: Iterable[str] = (),
) -> Design:
    """Encode ``parts`` into one DNA sequence, then repair it against ``rules``.

    ``forbidden`` motifs -- restriction sites, usually -- are kept out of the encoded parts
    on both strands. They are *not* removed from frozen parts: a cassette flanked by its
    own BsaI sites is the normal case, and silently rewriting it would be wrong.

    ``min_fraction`` drops codons rarer than that share of their amino acid's usage, which
    is what keeps the divergence pass from reaching for codons *E. coli* barely translates.

    Raises if a protein part cannot be spelled at all without a forbidden motif.
    """
    parts = list(parts)
    names = [part.name for part in parts]
    if len(set(names)) != len(names):
        raise ValueError("two parts share a name")
    forbidden = tuple(motif.upper() for motif in forbidden)

    spans: dict[str, tuple[int, int]] = {}
    editable: list[tuple[int, int]] = []
    sequence = ""
    for part in parts:
        start = len(sequence)
        if part.is_protein:
            references = [spans[name] for name in part.unlike if name in spans]
            missing = [name for name in part.unlike if name not in spans]
            if missing:
                raise ValueError(
                    f"part {part.name!r} is to be unlike {missing}, which "
                    f"{'is' if len(missing) == 1 else 'are'} not built before it"
                )
            reference_dna = [sequence[a:b] for a, b in references]
            assert part.peptide is not None
            for index, amino_acid in enumerate(part.peptide):
                options = usage.codons(amino_acid, min_fraction=min_fraction)
                at = 3 * index
                here = [ref[at : at + 3] for ref in reference_dna if at + 3 <= len(ref)]
                ranked = []
                for codon in options:
                    trial = sequence + codon
                    if _creates_forbidden(trial, forbidden, 3):
                        continue
                    ranked.append((
                        _codon_identity(codon, here),
                        _tail_repeat(trial, rules.repeat_seed),
                        -usage.fraction(codon),
                        codon,
                    ))
                if not ranked:
                    raise ValueError(
                        f"no codon for {amino_acid} at position {index + 1} of "
                        f"{part.name!r} avoids {list(forbidden)}"
                    )
                sequence += min(ranked)[3]
        else:
            assert part.dna is not None
            sequence += part.dna
        spans[part.name] = (start, len(sequence))
        if part.editable:
            editable.append((start, len(sequence)))

    if not repair:
        return Design(sequence, spans, complexity_report(sequence, rules, ignore=ignore), 0, 0)

    fixed = fix_complexity(
        sequence, editable, usage, forbidden=forbidden, rules=rules,
        min_fraction=min_fraction, rounds=rounds, ignore=ignore,
    )
    return Design(fixed.sequence, spans, fixed.issues, fixed.edits, fixed.rounds)


def fix_complexity(
    sequence: str,
    editable: Sequence[tuple[int, int]],
    usage: CodonUsage,
    *,
    forbidden: Sequence[str] = (),
    rules: Rules = TWIST_GENE_FRAGMENT,
    min_fraction: float = 0.10,
    rounds: int = 400,
    ignore: Iterable[str] = (),
) -> Design:
    """Make synonymous codon swaps until ``sequence`` stops breaking ``rules``.

    ``editable`` lists half-open ``(start, end)`` spans that may change; each must begin on
    a codon boundary of its own reading frame, because an edit swaps a whole codon.
    Everything outside them is protected.

    Each round takes the worst remaining issue, walks the codons inside it (and inside its
    partner region, for a repeat), and takes the **first** swap that lowers the total score
    -- first in a fixed order: by position, then by codon usage, then alphabetically. First
    rather than best is deliberate: it keeps the edit count down and the run time linear in
    the number of edits, and because the order is total the result is still reproducible.

    An issue no single codon can improve is left alone and reported, rather than raising:
    "two copies of one protein are 61 % identical" is a fact about the construct, not a bug
    to be worked around.
    """
    sequence = sequence.upper()
    forbidden = tuple(motif.upper() for motif in forbidden)
    original = sequence

    def codons_in(start: int, end: int) -> list[int]:
        """Codon starts inside an editable span that overlap the half-open region."""
        out: list[int] = []
        for first, last in editable:
            for position in range(first, last - 2, 3):
                if position + 3 > start and position < end:
                    out.append(position)
        return sorted(set(out))

    current = score(sequence, rules, ignore=ignore)
    stuck: set[tuple[str, int]] = set()
    # Where the walk through each region got to. Without it a 500 bp region is re-walked
    # from its left end after every accepted edit, which turns a linear repair quadratic.
    cursor: dict[tuple[str, int], int] = {}
    used = 0
    while used < rounds:
        issues = [issue for issue in complexity_report(sequence, rules, ignore=ignore)
                  if (issue.kind, issue.start) not in stuck]
        if not issues:
            break
        used += 1
        issue = issues[0]
        key = (issue.kind, issue.start)
        regions = [(issue.start, issue.end)]
        if issue.partner >= 0:
            regions.append((issue.partner, issue.partner + (issue.end - issue.start)))
        positions = [p for start, end in regions for p in codons_in(start, end)]
        improved = False
        for index in range(cursor.get(key, 0), len(positions)):
            position = positions[index]
            codon = sequence[position : position + 3]
            for alternative in usage.synonyms(codon, min_fraction=min_fraction):
                trial = sequence[:position] + alternative + sequence[position + 3 :]
                if forbidden:
                    lo, hi = max(0, position - 12), position + 15
                    if _sites(trial[lo:hi], forbidden) > _sites(sequence[lo:hi], forbidden):
                        continue      # the swap would create a site, not merely sit by one
                candidate = score(trial, rules, ignore=ignore)
                if candidate < current - 1e-12:
                    sequence, current, improved = trial, candidate, True
                    cursor[key] = index + 1
                    break
            if improved:
                break
        if not improved:
            stuck.add(key)

    edits = sum(1 for before, after in zip(original, sequence) if before != after)
    remaining = complexity_report(sequence, rules, ignore=ignore)
    return Design(sequence, {}, remaining, edits, used)


def identity_floor(
    peptide: str, usage: CodonUsage, *, min_fraction: float = 0.10
) -> float:
    """The lowest DNA identity two spellings of ``peptide`` could possibly have.

    Identity is positional, so the best case is reached codon by codon: for each residue,
    the two synonyms that share the fewest bases. Methionine and tryptophan have one codon
    each and contribute 3 matches whatever you do, and in most families the first two bases
    are fixed by the amino acid, so the floor for a typical protein is around 60 % -- not 0.

    Worth computing before trying to drive two copies further apart: if the measured
    identity is already at the floor, no codon choice will improve it and the remedy is
    architectural.
    """
    matches = 0
    for amino_acid in peptide.upper():
        options = usage.codons(amino_acid, min_fraction=min_fraction)
        matches += min(
            sum(1 for x, y in zip(first, second) if x == y)
            for first in options
            for second in options
        )
    return matches / (3 * len(peptide)) if peptide else 0.0


def identity(a: str, b: str) -> float:
    """Fraction of bases that match, position by position, over the shorter of the two.

    The measure vendors call homology and exact-repeat searches cannot see: two spellings
    of one protein share no long substring and are still ~79 % identical.
    """
    if not a or not b:
        return 0.0
    pairs = list(zip(a.upper(), b.upper()))
    return sum(1 for x, y in pairs if x == y) / len(pairs)
