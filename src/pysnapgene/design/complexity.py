"""Score a sequence the way a synthesis vendor does, before paying to find out.

A gene synthesis quote comes back with a complexity verdict, and the complaints are rarely
about the exact repeats that :mod:`pysnapgene.design.repeats` removes. The expensive ones
are statistical: what *fraction* of the sequence sits in some short repeat, and whether two
regions are similar enough to recombine even though no long substring is shared. Two copies
of the same protein, codon-differentiated so that they share no run longer than 12 bp, are
still around 79 % identical base for base -- which a vendor reads as one long repeat.

The rules here are modelled on the complexity warnings Twist returns for gene fragments:

    High repeat density detected -- more than 15 % of your sequence is composed of small
    repeats (9 bp or longer).

    Direct repeat with >=68 C Tm -- direct repeats with a melting temperature at or over
    68 C increase the risk of recombination.

    Long, low homology repeat region detected.

They are an approximation. The vendor's algorithm is not published, the thresholds here are
read off its messages, and a sequence that passes this check can still be quoted as complex.
What the check is good for is comparing two designs: a sequence that scores better here is
better by the same measures the vendor names.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from typing import NamedTuple

__all__ = [
    "Issue",
    "Rules",
    "TWIST_GENE_FRAGMENT",
    "complexity_report",
    "repeat_density",
    "repeat_coverage",
    "exact_repeats",
    "homologous_regions",
    "Homology",
    "duplex_tm",
    "gc_fraction",
    "gc_extremes",
    "homopolymers",
    "score",
    "summarise",
]


class Issue(NamedTuple):
    """One thing a vendor would complain about.

    ``start`` and ``end`` are 0-based, half-open, and bound the whole offending region --
    for a repeat, the first copy; ``partner`` gives the start of the second one, or -1
    where the issue is not a pair. ``value`` and ``limit`` are in the units named by
    ``kind``, so an issue can be ranked by how far over the line it is.
    """

    kind: str
    detail: str
    start: int
    end: int
    value: float
    limit: float
    partner: int = -1

    @property
    def excess(self) -> float:
        """How far past the limit, as a fraction of the limit. Negative means inside it."""
        if self.limit == 0:
            return self.value
        return (self.value - self.limit) / abs(self.limit)

    def __str__(self) -> str:
        where = f"{self.start + 1}-{self.end}"
        if self.partner >= 0:
            where += f" and {self.partner + 1}-{self.partner + (self.end - self.start)}"
        return f"{self.kind} at {where}: {self.detail}"


class Rules(NamedTuple):
    """Thresholds a synthesised sequence is held to.

    The defaults are :data:`TWIST_GENE_FRAGMENT`. Build your own to match another vendor,
    or to tighten one knob: ``TWIST_GENE_FRAGMENT._replace(max_repeat_density=0.12)``.
    """

    name: str = "Twist gene fragment"
    #: length the vendor will synthesise as one fragment
    min_length: int = 300
    max_length: int = 5000
    #: a repeat this long or longer counts towards the density figure
    repeat_seed: int = 9
    #: fraction of bases allowed to sit inside such a repeat
    max_repeat_density: float = 0.15
    #: longest exact direct repeat allowed, in bp -- the recombination limit, not a vendor one
    max_exact_repeat: int = 13
    #: an inexact repeat is scored over a window this wide
    homology_window: int = 100
    #: and flagged at or above this identity within the window
    max_homology: float = 0.70
    #: a repeat duplex at or above this melting temperature is flagged
    max_repeat_tm: float = 68.0
    #: GC over the whole sequence
    min_gc: float = 0.25
    max_gc: float = 0.65
    #: and over a sliding window, where vendors allow far more local swing
    gc_window: int = 50
    min_window_gc: float = 0.10
    max_window_gc: float = 0.90
    #: longest run of a single base
    max_homopolymer: int = 13


TWIST_GENE_FRAGMENT = Rules()


def gc_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def duplex_tm(length: int, gc: float, identity: float = 1.0, *, salt: float = 0.05) -> float:
    """Melting temperature of a duplex ``length`` bp long, in Celsius.

    The long-duplex formula ``81.5 + 16.6 log10[Na+] + 0.41 %GC - 500/length``, with one
    degree subtracted per percent of mismatch, which is the usual rule of thumb for an
    imperfect duplex. It is an estimate: above ~50 bp a nearest-neighbour calculation is
    not applicable either, and vendors do not publish which formula they use.
    """
    if length <= 0:
        return float("nan")
    return (
        81.5
        + 16.6 * math.log10(salt)
        + 0.41 * (100 * gc)
        - 500 / length
        - 100 * (1 - identity)
    )


def repeat_coverage(sequence: str, seed: int) -> set[int]:
    """Every position that sits inside some exact direct repeat of at least ``seed`` bp.

    Both copies count, which is what "composed of small repeats" means: the whole of a
    duplicated stretch is repeat, not just its second appearance.
    """
    sequence = sequence.upper()
    if seed <= 0 or len(sequence) < seed:
        return set()
    positions: dict[str, list[int]] = {}
    for i in range(len(sequence) - seed + 1):
        positions.setdefault(sequence[i : i + seed], []).append(i)
    covered: set[int] = set()
    for starts in positions.values():
        if len(starts) < 2:
            continue
        for start in starts:
            covered.update(range(start, start + seed))
    return covered


def repeat_density(sequence: str, seed: int = 9) -> float:
    """Fraction of the sequence that sits inside an exact repeat of at least ``seed`` bp."""
    if not sequence:
        return 0.0
    return len(repeat_coverage(sequence, seed)) / len(sequence)


class Pair(NamedTuple):
    """An exact direct repeat: the same string at two positions."""

    length: int
    first: int
    second: int

    def motif(self, sequence: str) -> str:
        return sequence[self.first : self.first + self.length]


def exact_repeats(sequence: str, min_length: int = 9) -> list[Pair]:
    """Every maximal exact direct repeat pair of at least ``min_length`` bp, longest first.

    *Maximal* means extended as far as it goes in both directions, so a 30 bp repeat is
    reported once rather than as 22 overlapping 9-mers.
    """
    sequence = sequence.upper()
    if min_length <= 0 or len(sequence) < min_length:
        return []
    seeds: dict[str, list[int]] = {}
    for i in range(len(sequence) - min_length + 1):
        seeds.setdefault(sequence[i : i + min_length], []).append(i)
    found: set[Pair] = set()
    for starts in seeds.values():
        if len(starts) < 2:
            continue
        for a in range(len(starts)):
            for b in range(a + 1, len(starts)):
                i, j = starts[a], starts[b]
                left = 0
                while i - left > 0 and sequence[i - left - 1] == sequence[j - left - 1]:
                    left += 1
                right = min_length
                while j + right < len(sequence) and sequence[i + right] == sequence[j + right]:
                    right += 1
                found.add(Pair(right + left, i - left, j - left))
    return sorted(found, key=lambda pair: (-pair.length, pair.first, pair.second))


class Homology(NamedTuple):
    """Two regions similar enough to recombine, without being identical.

    ``offset`` is the distance between them, ``identity`` the fraction of matching bases
    over the window, and ``tm`` the estimated melting temperature of the duplex they would
    form.
    """

    first: int
    second: int
    length: int
    identity: float
    tm: float

    @property
    def offset(self) -> int:
        return self.second - self.first


def homologous_regions(
    sequence: str,
    *,
    window: int = 100,
    min_identity: float = 0.70,
    seed: int = 9,
    step: int = 3,
) -> list[Homology]:
    """Inexact direct repeats: windows of ``sequence`` that are similar but not identical.

    Two copies of one protein, spelled with different codons, share no long substring and
    are still ~79 % identical base for base. That is invisible to an exact-repeat search
    and is exactly what a vendor flags as a long repeat region, so it is searched for
    differently: every offset carrying at least one shared ``seed``-mer is scanned as an
    ungapped diagonal, and each run of windows at or above ``min_identity`` is reported
    once, as its widest extent with the identity of its best window.

    Two limitations, both deliberate. The scan is **ungapped**, which is right for two
    spellings of one protein -- they stay in register -- and wrong for two homologous genes
    that do not. And an offset is only examined when some exact ``seed``-mer is shared, so
    homology spread evenly enough to leave no seed is missed: lower ``seed`` to find it, at
    a cost in speed, since every extra offset is another pass over the sequence.
    """
    sequence = sequence.upper()
    length = len(sequence)
    if length < window or window <= 0:
        return []
    positions: dict[str, list[int]] = {}
    for i in range(length - seed + 1):
        positions.setdefault(sequence[i : i + seed], []).append(i)
    offsets: set[int] = set()
    for starts in positions.values():
        if len(starts) < 2:
            continue
        for a in range(len(starts)):
            for b in range(a + 1, len(starts)):
                offsets.add(starts[b] - starts[a])

    found: list[Homology] = []
    for offset in sorted(o for o in offsets if o >= window):
        matches = [1 if sequence[i] == sequence[i + offset] else 0
                   for i in range(length - offset)]
        if len(matches) < window:
            continue
        running = sum(matches[:window])
        best: tuple[int, int, int] | None = None  # start, end, best count
        for i in range(len(matches) - window + 1):
            if i:
                running += matches[i + window - 1] - matches[i - 1]
            above = running >= min_identity * window
            if above and best is None:
                best = (i, i + window, running)
            elif above and best is not None:
                best = (best[0], i + window, max(best[2], running))
            elif not above and best is not None:
                found.append(_homology(sequence, best, offset, window))
                best = None
        if best is not None:
            found.append(_homology(sequence, best, offset, window))
    found.sort(key=lambda h: (-h.tm, -h.length, h.first))
    # one offset can produce several runs; keep them, but drop a run contained in another
    kept: list[Homology] = []
    for hit in found:
        if any(
            other.offset == hit.offset
            and other.first <= hit.first
            and other.first + other.length >= hit.first + hit.length
            for other in kept
        ):
            continue
        kept.append(hit)
    return kept


def _homology(sequence: str, run: tuple[int, int, int], offset: int, window: int) -> Homology:
    start, end, best_count = run
    region = sequence[start:end]
    identity = best_count / window
    return Homology(
        first=start,
        second=start + offset,
        length=end - start,
        identity=identity,
        tm=duplex_tm(min(end - start, window), gc_fraction(region), identity),
    )


def gc_extremes(sequence: str, window: int, low: float, high: float) -> list[tuple[int, float]]:
    """Windows whose GC falls outside ``low``-``high``, as ``(start, gc)``."""
    sequence = sequence.upper()
    if len(sequence) < window or window <= 0:
        return []
    count = sum(1 for base in sequence[:window] if base in "GC")
    out: list[tuple[int, float]] = []
    for i in range(len(sequence) - window + 1):
        if i:
            count += (sequence[i + window - 1] in "GC") - (sequence[i - 1] in "GC")
        gc = count / window
        if gc < low or gc > high:
            out.append((i, gc))
    return out


def _merge_windows(
    windows: list[tuple[int, float]], width: int
) -> list[tuple[int, int, float]]:
    """Collapse overlapping flagged windows into regions, keeping the worst value.

    A sliding window reports the same GC-rich stretch once per offset; a reader wants the
    stretch, and a repair wants one target rather than forty.
    """
    merged: list[tuple[int, int, float]] = []
    for start, value in windows:
        if merged and start <= merged[-1][1]:
            first, last, worst = merged[-1]
            keep = worst if abs(worst - 0.5) >= abs(value - 0.5) else value
            merged[-1] = (first, max(last, start + width), keep)
        else:
            merged.append((start, start + width, value))
    return merged


def homopolymers(sequence: str, max_run: int) -> list[tuple[int, int, str]]:
    """Runs of a single base longer than ``max_run``, as ``(start, length, base)``."""
    sequence = sequence.upper()
    out: list[tuple[int, int, str]] = []
    start = 0
    for i in range(1, len(sequence) + 1):
        if i == len(sequence) or sequence[i] != sequence[start]:
            if i - start > max_run:
                out.append((start, i - start, sequence[start]))
            start = i
    return out


def complexity_report(
    sequence: str,
    rules: Rules = TWIST_GENE_FRAGMENT,
    *,
    ignore: Iterable[str] = (),
) -> list[Issue]:
    """Everything in ``sequence`` that breaks ``rules``, worst first.

    ``ignore`` drops whole issue kinds, for the case where one of them is a deliberate
    part of the design -- a vendor's length limit does not apply to a sequence that will
    be ordered in pieces, for instance.
    """
    sequence = sequence.upper()
    skip = set(ignore)
    issues: list[Issue] = []

    def add(issue: Issue) -> None:
        if issue.kind not in skip:
            issues.append(issue)

    length = len(sequence)
    if length < rules.min_length:
        add(Issue("length", f"{length} bp is below the {rules.min_length} bp minimum",
                  0, length, rules.min_length - length, 0))
    if length > rules.max_length:
        add(Issue("length", f"{length} bp is above the {rules.max_length} bp maximum",
                  0, length, length, rules.max_length))

    density = repeat_density(sequence, rules.repeat_seed)
    if density > rules.max_repeat_density:
        add(Issue(
            "repeat_density",
            f"{100 * density:.1f} % of the sequence sits in a repeat of "
            f">={rules.repeat_seed} bp (limit {100 * rules.max_repeat_density:.0f} %)",
            0, length, density, rules.max_repeat_density))

    for pair in exact_repeats(sequence, rules.max_exact_repeat + 1):
        motif = pair.motif(sequence)
        add(Issue(
            "exact_repeat",
            f"{pair.length} bp exact direct repeat, Tm "
            f"{duplex_tm(pair.length, gc_fraction(motif)):.0f} C",
            pair.first, pair.first + pair.length, pair.length, rules.max_exact_repeat,
            partner=pair.second))

    for hit in homologous_regions(
        sequence, window=rules.homology_window, min_identity=rules.max_homology,
        seed=rules.repeat_seed,
    ):
        over_tm = hit.tm >= rules.max_repeat_tm
        add(Issue(
            "homology",
            f"{hit.length} bp region {100 * hit.identity:.0f} % identical to the one "
            f"{hit.offset} bp away, Tm {hit.tm:.0f} C"
            + (" -- over the Tm limit" if over_tm else ""),
            hit.first, hit.first + hit.length,
            hit.identity, rules.max_homology, partner=hit.second))

    gc = gc_fraction(sequence)
    if not rules.min_gc <= gc <= rules.max_gc:
        add(Issue("gc_global", f"overall GC {100 * gc:.1f} %", 0, length, gc,
                  rules.max_gc if gc > rules.max_gc else rules.min_gc))
    for first, last, worst in _merge_windows(
        gc_extremes(sequence, rules.gc_window, rules.min_window_gc, rules.max_window_gc),
        rules.gc_window,
    ):
        add(Issue(
            "gc_window",
            f"GC reaches {100 * worst:.0f} % in a {rules.gc_window} bp window",
            first, last, worst,
            rules.max_window_gc if worst > rules.max_window_gc else rules.min_window_gc))

    for start, run, base in homopolymers(sequence, rules.max_homopolymer):
        add(Issue("homopolymer", f"{run} x {base}", start, start + run, run,
                  rules.max_homopolymer))

    issues.sort(key=lambda issue: (-issue.excess, issue.kind, issue.start))
    return issues


def score(sequence: str, rules: Rules = TWIST_GENE_FRAGMENT,
          *, ignore: Iterable[str] = ()) -> float:
    """One number for how badly ``sequence`` breaks ``rules``; 0 is clean.

    The sum of every issue's excess over its limit, so any repair that lowers it is an
    improvement by at least one of the vendor's own measures. Comparable only between
    sequences checked against the same rules.
    """
    return sum(max(issue.excess, 0.0) for issue in complexity_report(sequence, rules,
                                                                     ignore=ignore))


def summarise(sequence: str, rules: Rules = TWIST_GENE_FRAGMENT,
              *, ignore: Iterable[str] = ()) -> str:
    """A short human-readable verdict, for printing in a notebook."""
    issues = complexity_report(sequence, rules, ignore=ignore)
    if not issues:
        return f"{len(sequence)} bp, clean against {rules.name}"
    lines = [f"{len(sequence)} bp, {len(issues)} issue(s) against {rules.name}:"]
    lines += [f"  {issue}" for issue in issues]
    return "\n".join(lines)
