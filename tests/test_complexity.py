"""Vendor-style complexity checks.

The sequences here are planted: a known repeat at a known offset, a known GC island, a
known homopolymer. A checker that cannot find what was deliberately put in front of it is
no use on a real construct.
"""

from __future__ import annotations

import random

from pysnapgene.design import (
    TWIST_GENE_FRAGMENT,
    complexity_report,
    duplex_tm,
    exact_repeats,
    gc_fraction,
    homologous_regions,
    homopolymers,
    repeat_density,
    score,
)
from pysnapgene.design.complexity import gc_extremes, repeat_coverage


def _random_dna(length: int, seed: int, gc: float = 0.5) -> str:
    rng = random.Random(seed)
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ATGC", weights=weights, k=length))


def test_repeat_density_counts_both_copies():
    block = _random_dna(30, 1)
    sequence = block + _random_dna(40, 2) + block
    covered = repeat_coverage(sequence, 9)
    assert len(covered) == 60, "both copies of a 30 bp block are repeat"
    assert repeat_density(sequence, 9) == 60 / len(sequence)


def test_random_sequence_has_little_repeat():
    assert repeat_density(_random_dna(1500, 3), 9) < 0.05


def test_exact_repeats_are_maximal():
    block = _random_dna(25, 4)
    sequence = _random_dna(50, 5) + block + _random_dna(60, 6) + block + _random_dna(50, 7)
    hits = exact_repeats(sequence, 9)
    assert hits[0].length == 25, "a 25 bp repeat should be reported once, at full length"
    assert hits[0].motif(sequence) == block


def test_homologous_regions_finds_an_inexact_repeat():
    rng = random.Random(8)
    original = _random_dna(300, 9)
    # 85 % identical: change one base in six
    mutated = "".join(
        rng.choice("ATGC") if rng.random() < 0.15 else base for base in original
    )
    sequence = original + _random_dna(200, 10) + mutated
    hits = homologous_regions(sequence, window=100, min_identity=0.70)
    assert hits, "an 85 % identical 300 bp region should be found"
    best = hits[0]
    assert best.offset == 500
    assert 0.70 <= best.identity <= 1.0
    assert best.length >= 100


def test_homologous_regions_ignores_unrelated_sequence():
    sequence = _random_dna(800, 11)
    assert homologous_regions(sequence, window=100, min_identity=0.70) == []


def test_homology_is_invisible_to_an_exact_search():
    """The case this module exists for: no long substring, still one long repeat."""
    rng = random.Random(12)
    original = _random_dna(400, 13)
    # 80 % identity, mutated at random positions as real codon divergence is
    mutated = "".join(
        rng.choice([b for b in "ATGC" if b != base]) if rng.random() < 0.20 else base
        for base in original
    )
    sequence = original + mutated
    assert not exact_repeats(sequence, 30), "no exact repeat anywhere near the region length"
    hits = homologous_regions(sequence, window=100, min_identity=0.70)
    assert hits and hits[0].offset == 400
    assert hits[0].identity >= 0.70


def test_homology_search_needs_a_shared_seed():
    """A documented blind spot, asserted so it cannot change silently.

    Offsets are found by exact k-mer seeds. Divergence spread perfectly evenly -- every
    fifth base changed -- leaves no 9-mer in common, and the region is missed although it
    is 80 % identical. Real codon divergence is irregular and does leave seeds; a sequence
    engineered to avoid them needs a smaller ``seed``, at a cost in speed.
    """
    original = _random_dna(400, 13)
    even = "".join(
        ("A" if base != "A" else "C") if i % 5 == 0 else base
        for i, base in enumerate(original)
    )
    sequence = original + even
    assert homologous_regions(sequence, window=100, min_identity=0.70) == []
    assert homologous_regions(sequence, window=100, min_identity=0.70, seed=4)


def test_duplex_tm_rises_with_length_and_gc_and_falls_with_mismatch():
    assert duplex_tm(100, 0.5) > duplex_tm(20, 0.5)
    assert duplex_tm(100, 0.7) > duplex_tm(100, 0.4)
    assert duplex_tm(100, 0.5, identity=0.8) < duplex_tm(100, 0.5)


def test_gc_windows_are_merged_into_regions():
    sequence = _random_dna(100, 14, gc=0.5) + "GC" * 40 + _random_dna(100, 15, gc=0.5)
    flagged = gc_extremes(sequence, 50, 0.10, 0.90)
    assert flagged, "a 80 bp GC island should break a 90 % window limit"
    issues = [i for i in complexity_report(sequence) if i.kind == "gc_window"]
    assert len(issues) == 1, "overlapping windows over one island are one issue"
    assert issues[0].start <= 100 < issues[0].end


def test_homopolymers():
    assert homopolymers("AAAA" + "T" * 15 + "GG", 13) == [(4, 15, "T")]
    assert homopolymers("A" * 13, 13) == []


def test_gc_fraction():
    assert gc_fraction("GGCC") == 1.0
    assert gc_fraction("ATAT") == 0.0
    assert gc_fraction("") == 0.0


def test_a_clean_sequence_scores_zero():
    sequence = _random_dna(1000, 16)
    assert complexity_report(sequence) == []
    assert score(sequence) == 0.0


def test_report_names_the_rule_and_ranks_by_excess():
    block = _random_dna(40, 17)
    sequence = (block + _random_dna(30, 18)) * 6          # heavily repetitive
    issues = complexity_report(sequence)
    kinds = {issue.kind for issue in issues}
    assert {"repeat_density", "exact_repeat"} <= kinds
    assert [i.excess for i in issues] == sorted((i.excess for i in issues), reverse=True)
    assert score(sequence) > 0


def test_ignore_drops_a_kind():
    short = _random_dna(100, 19)
    assert any(i.kind == "length" for i in complexity_report(short))
    assert not any(i.kind == "length" for i in complexity_report(short, ignore=["length"]))


def test_rules_can_be_tightened():
    sequence = _random_dna(1000, 20)
    strict = TWIST_GENE_FRAGMENT._replace(max_repeat_density=0.0)
    assert complexity_report(sequence) == []
    assert any(i.kind == "repeat_density" for i in complexity_report(sequence, strict))


def test_issue_str_names_both_copies():
    block = _random_dna(30, 21)
    sequence = block + _random_dna(400, 22) + block
    issue = next(i for i in complexity_report(sequence) if i.kind == "exact_repeat")
    assert "and" in str(issue)
    assert issue.partner > issue.start
