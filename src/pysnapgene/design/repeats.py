"""Find direct repeats in a sequence, and remove them with synonymous edits.

Why this matters: a plasmid carrying the same sequence twice can recombine across the two
copies and delete everything between them. When what sits between them is the counter-
selection marker, the deletion product survives selection and grows like a real clone, so
the repeat has to be designed out rather than screened out.

The useful threshold is the length at which a repeat stops being expected by chance. A
~6 kb plasmid contains on the order of 7.6 coincidental 10-mers but only ~0.1 coincidental
13-mers, so 13 bp is about the shortest limit that is actually reachable by design.
"""

from __future__ import annotations

import random
from collections.abc import Iterable, Sequence
from typing import NamedTuple

from .codons import CodonUsage, site_free

__all__ = [
    "Repeat",
    "RepeatEdits",
    "DEFAULT_REPEAT_LIMIT",
    "longest_common_substring",
    "longest_repeat",
    "break_repeats",
]

#: The length above which a repeat is no longer expected by chance in a plasmid.
DEFAULT_REPEAT_LIMIT = 13


class Repeat(NamedTuple):
    """A sequence occurring at least twice, and the two positions found."""

    motif: str
    first: int
    second: int

    def __len__(self) -> int:
        return len(self.motif)


class RepeatEdits(NamedTuple):
    """Result of :func:`break_repeats`."""

    sequence: str
    unfixable: list[str]
    edits: int


def longest_common_substring(a: str, b: str) -> int:
    """Length of the longest substring shared by ``a`` and ``b``."""
    if not a or not b:
        return 0
    previous = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        current = [0] * (len(b) + 1)
        a_i = a[i - 1]
        for j in range(1, len(b) + 1):
            if a_i == b[j - 1]:
                current[j] = previous[j - 1] + 1
                best = max(best, current[j])
        previous = current
    return best


def longest_repeat(
    sequence: str,
    min_length: int = DEFAULT_REPEAT_LIMIT + 1,
    max_length: int = 200,
    *,
    circular: bool = False,
) -> Repeat | None:
    """The longest substring of ``sequence`` occurring at least twice, or None.

    Only repeats of at least ``min_length`` are reported, longest first, so
    ``longest_repeat(seq, limit + 1)`` answers "is there a repeat longer than *limit*".
    ``max_length`` bounds the search; raise it if repeats longer than 200 bp are expected.

    With ``circular`` the search wraps the origin, which is what a plasmid requires: two
    halves of a repeat can sit either side of position 0.
    """
    if min_length < 1:
        raise ValueError("min_length must be at least 1")
    length = len(sequence)
    # a linear repeat needs two disjoint copies, so it cannot exceed half the sequence;
    # a circular one may overlap itself as a tandem repeat
    widest = min(max_length, length - 1 if circular else length // 2)
    for size in range(widest, min_length - 1, -1):
        haystack = sequence + sequence[: size - 1] if circular else sequence
        seen: dict[str, int] = {}
        for i in range(len(haystack) - size + 1):
            motif = haystack[i : i + size]
            start = seen.get(motif)
            if start is None:
                seen[motif] = i
                continue
            # in a circular search the same occurrence reappears one turn later
            if circular and (i - start) % length == 0:
                continue
            return Repeat(motif, start % length, i % length)
    return None


def break_repeats(
    sequence: str,
    editable: Sequence[tuple[int, int]],
    usage: CodonUsage,
    *,
    rng: random.Random | None = None,
    max_length: int = DEFAULT_REPEAT_LIMIT,
    forbidden: Iterable[str] = (),
    min_fraction: float = 0.08,
    rounds: int = 400,
    circular: bool = False,
) -> RepeatEdits:
    """Make synonymous edits until no repeat longer than ``max_length`` remains.

    ``editable`` lists half-open ``(start, end)`` spans that may be changed; each must be
    **in frame**, i.e. start on a codon boundary of the reading frame it belongs to, since
    the edits swap whole codons. Anything outside those spans is left alone, which is how a
    preserved cassette or a backbone is protected.

    Edits are chosen to be minimal in two senses: one codon at a time, and starting from the
    middle of the repeat outwards, so the repeat is split rather than shortened from an end.
    When a repeat has no editable position in either copy it is reported in ``unfixable``
    rather than raising, because that is a design decision for the caller: the two copies
    are both in fixed sequence and something else has to give.
    """
    rng = rng or random.Random()
    forbidden = tuple(forbidden)

    def span_of(position: int) -> tuple[int, int] | None:
        return next(((a, b) for a, b in editable if a <= position < b), None)

    original = sequence
    bases = list(sequence)
    unfixable: list[str] = []
    for _ in range(rounds):
        hit = longest_repeat("".join(bases), max_length + 1, circular=circular)
        if hit is None:
            break
        done = False
        for occurrence in (hit.second, hit.first):  # prefer editing the later copy
            for offset in range(len(hit.motif) // 2, -1, -1):
                for position in (
                    occurrence + offset,
                    occurrence + len(hit.motif) - 1 - offset,
                ):
                    span = span_of(position)
                    if span is None:
                        continue
                    start = span[0] + ((position - span[0]) // 3) * 3
                    codon = "".join(bases[start : start + 3])
                    alternatives = usage.synonyms(codon, min_fraction=min_fraction)
                    rng.shuffle(alternatives)
                    for alternative in alternatives:
                        trial = bases[:]
                        trial[start : start + 3] = list(alternative)
                        candidate = "".join(trial)
                        if forbidden and not site_free(candidate, forbidden):
                            continue
                        now = longest_repeat(candidate, max_length + 1, circular=circular)
                        if now is None or len(now.motif) < len(hit.motif):
                            bases = trial
                            done = True
                            break
                    if done:
                        break
                if done:
                    break
            if done:
                break
        if not done:
            unfixable.append(hit.motif)
            break

    cleaned = "".join(bases)
    edits = sum(1 for before, after in zip(original, cleaned) if before != after)
    return RepeatEdits(cleaned, unfixable, edits)
