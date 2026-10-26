"""Homology-based assembly: In-Fusion, Gibson, HiFi."""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from pysnapgene.cloning import (
    INFUSION_MIN_HOMOLOGY,
    ambiguous_homology,
    homology_assemble,
    infusion,
)

# Three fragments that close into a circle. Each arm is 20 bp and distinctive, so there is no
# accidental homology anywhere else.
ARM_AB = "GGATCAGTTCACGATTGCAA"
ARM_BC = "CCTTAGGACAAGTTCCGATT"
ARM_CA = "TTGCCAATGGTCAAGGTTAC"
BODY_A = "ACGTACGTACGTACGTACGTACGTAC"
BODY_B = "TTTTCCCCAAAAGGGGTTTTCCCCAA"
BODY_C = "GCGCATATGCGCGCTTAATTAAGCGC"

A = ARM_CA + BODY_A + ARM_AB
B = ARM_AB + BODY_B + ARM_BC
C = ARM_BC + BODY_C + ARM_CA


def test_two_fragments_close_into_a_circle():
    first = ARM_CA + BODY_A + ARM_AB
    second = ARM_AB + BODY_B + ARM_CA
    assembly = infusion([first, second])
    assert assembly.circular
    assert len(assembly) == len(first) + len(second) - 2 * len(ARM_AB)
    assert assembly.sequence == ARM_CA + BODY_A + ARM_AB + BODY_B


def test_each_homology_appears_once_in_the_product():
    assembly = infusion([A, B, C])
    for arm in (ARM_AB, ARM_BC, ARM_CA):
        assert assembly.sequence.count(arm) == 1, f"{arm} is duplicated in the product"


def test_product_length_is_the_sum_less_the_homology():
    assembly = infusion([A, B, C])
    total_homology = sum(len(j) for j in assembly.junctions)
    assert len(assembly) == len(A) + len(B) + len(C) - total_homology


def test_junctions_report_where_the_homology_landed():
    assembly = infusion([A, B, C])
    assert len(assembly.junctions) == 3  # two joins plus the closure
    for junction in assembly.junctions:
        start, end = junction.start, junction.start + len(junction)
        assert assembly.sequence[start:end] == junction.homology
    assert [(j.upstream, j.downstream) for j in assembly.junctions] == [(0, 1), (1, 2), (2, 0)]


def test_offsets_map_a_fragment_position_into_the_product():
    assembly = infusion([A, B, C])
    for index, fragment in enumerate((A, B, C)):
        for position in (0, 5, len(fragment) - 1):
            landed = assembly.locate(index, position)
            assert assembly.sequence[landed] == fragment[position], (
                f"fragment {index} base {position} does not land where locate() says"
            )


def test_a_linear_assembly_is_not_closed():
    assembly = homology_assemble([A, B], circular=False)
    assert not assembly.circular
    assert assembly.sequence == A + B[len(ARM_AB) :]
    assert len(assembly.junctions) == 1


def test_a_fragment_the_wrong_way_round_is_named_as_such():
    """The usual mistake, and the error should say so rather than just "no homology"."""
    flipped = str(Seq(B).reverse_complement())
    with pytest.raises(ValueError, match="wrong way round"):
        infusion([A, flipped, C])


def test_homology_shorter_than_the_minimum_is_refused():
    short = ARM_AB[-8:] + BODY_B + ARM_CA
    with pytest.raises(ValueError, match="only 8 bp"):
        infusion([ARM_CA + BODY_A + ARM_AB, short])


def test_unrelated_fragments_are_refused():
    with pytest.raises(ValueError, match="no terminal sequence"):
        infusion([A, "TTTTTTTTTTTTTTTTTTTTTTTT"])


def test_a_circle_that_does_not_close_is_refused():
    with pytest.raises(ValueError):
        infusion([ARM_CA + BODY_A + ARM_AB, ARM_AB + BODY_B + "GGGGGGGGGGGGGGGGGGGG"])


def test_longest_overlap_is_taken():
    """With 25 bp of shared sequence the join uses all of it, not just the 15 bp minimum."""
    long_arm = "GGATCAGTTCACGATTGCAACGTAA"
    assembly = infusion([ARM_CA + BODY_A + long_arm, long_arm + BODY_B + ARM_CA])
    assert len(assembly.junctions[0]) == len(long_arm)


def test_max_homology_below_the_real_overlap_refuses_the_junction():
    """``max_homology`` bounds the search, it does not truncate the join.

    A shorter suffix of the left fragment is not a prefix of the right one, so capping the
    search below the real homology finds nothing at all rather than a partial match. Worth
    knowing: the bound is there to stop a coincidental long match, not to trim a real one.
    """
    long_arm = "GGATCAGTTCACGATTGCAACGTAA"  # 25 bp
    pair = [ARM_CA + BODY_A + long_arm, long_arm + BODY_B + ARM_CA]
    assert len(infusion(pair).junctions[0]) == 25
    with pytest.raises(ValueError):
        infusion(pair, max_homology=18)


def test_repeated_homology_is_reported():
    """A homology arm that also occurs elsewhere makes the assembly unpredictable."""
    clean = infusion([A, B, C])
    assert ambiguous_homology(clean) == []
    # plant a second copy of the closing arm inside a fragment body
    planted_b = ARM_AB + BODY_B + ARM_CA + ARM_BC
    assembly = infusion([A, planted_b, C])
    assert assembly.sequence.count(ARM_CA) == 2
    flagged = {j.homology for j in ambiguous_homology(assembly)}
    assert ARM_CA in flagged, "the duplicated arm should be reported"


def test_seq_objects_are_accepted():
    assembly = infusion([Seq(A), Seq(B), Seq(C)])
    assert assembly.sequence == infusion([A, B, C]).sequence


def test_lowercase_input_is_normalised():
    assembly = infusion([A.lower(), B.lower(), C.lower()])
    assert assembly.sequence == infusion([A, B, C]).sequence


def test_no_fragments_is_an_error():
    with pytest.raises(ValueError, match="no fragments"):
        infusion([])


def test_infusion_default_matches_the_kit():
    assert INFUSION_MIN_HOMOLOGY == 15
