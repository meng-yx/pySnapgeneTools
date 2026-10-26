"""Codon optimisation and repeat removal."""

from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq

from pysnapgene.design import (
    CodonUsage,
    break_repeats,
    codon_optimize,
    longest_common_substring,
    longest_repeat,
    site_free,
    translate,
)

# A miniature usage table: four amino acids with two to four synonyms each, enough to
# exercise sampling, argmax, the fraction floor and synonymous substitution.
USAGE = CodonUsage(
    {
        "GGT": ("G", 0.34),
        "GGC": ("G", 0.40),
        "GGA": ("G", 0.11),
        "GGG": ("G", 0.15),
        "AGT": ("S", 0.15),
        "AGC": ("S", 0.28),
        "TCT": ("S", 0.15),
        "TCC": ("S", 0.15),
        "TCA": ("S", 0.12),
        "TCG": ("S", 0.15),
        "GAT": ("D", 0.63),
        "GAC": ("D", 0.37),
        "AAA": ("K", 0.74),
        "AAG": ("K", 0.26),
    }
)


def test_translate_trims_a_partial_codon():
    assert translate("ATGAAA") == "MK"
    assert translate("ATGAAAG") == "MK"
    assert translate("") == ""


def test_site_free_checks_both_strands():
    assert site_free("AAAAGGGGAAAA", ["GGTCTC"])
    assert not site_free("AAAGGTCTCAAA", ["GGTCTC"])
    # GAGACC is the reverse complement of GGTCTC, so a BsaI site on the bottom strand
    assert not site_free("AAAGAGACCAAA", ["GGTCTC"])
    assert site_free("AAAGAGACCAAA", ["GGTCTC"], both_strands=False)


def test_codons_are_ordered_by_usage():
    assert USAGE.codons("G")[0] == "GGC"
    assert USAGE.best("G") == "GGC"
    assert USAGE.fraction("GGC") == 0.40


def test_fraction_floor_never_leaves_an_amino_acid_unencodable():
    """A threshold above every synonym must still return the commonest one."""
    assert USAGE.codons("S", min_fraction=0.9) == ["AGC"]


def test_unknown_amino_acid_is_reported():
    with pytest.raises(KeyError, match="amino acid"):
        USAGE.codons("W")


def test_argmax_is_deterministic_and_encodes_the_protein():
    first = codon_optimize("GGSDK", USAGE, argmax=True)
    second = codon_optimize("GGSDK", USAGE, argmax=True)
    assert first.dna == second.dna == "GGCGGCAGCGATAAA"
    assert translate(first.dna) == "GGSDK"


def test_sampling_encodes_the_protein():
    rng = random.Random(0)
    choice = codon_optimize("GGGSGGSDGGSGGSGSGG", USAGE, rng=rng)
    assert translate(choice.dna) == "GGGSGGSDGGSGGSGSGG"


def test_avoid_keeps_the_second_spelling_away_from_the_first():
    """The reason this function exists: two copies of an element must differ."""
    protein = "GGGSGGSDGGSGGSGSGG"
    rng = random.Random(1)
    first = codon_optimize(protein, USAGE, rng=rng)
    second = codon_optimize(protein, USAGE, rng=rng, avoid=[first.dna], target=6)
    assert translate(second.dna) == protein
    assert second.shared == longest_common_substring(first.dna, second.dna)
    assert second.shared < len(protein) * 3


def test_forbidden_motifs_are_absent_from_the_result():
    rng = random.Random(2)
    # GGTGGC and its complement are reachable spellings of "GG", so this really constrains
    choice = codon_optimize("GGSGGK", USAGE, rng=rng, forbidden=["GGTGGC"])
    assert site_free(choice.dna, ["GGTGGC"])


def test_impossible_constraint_raises():
    with pytest.raises(ValueError, match="no spelling"):
        codon_optimize("KK", USAGE, rng=random.Random(3), forbidden=["AA"], tries=20)


# --- repeats ----------------------------------------------------------------


def test_longest_repeat_finds_the_longer_of_two():
    seq = "AAACGTACGTAAA" + "TTTT" + "AAACGTACGTAAA"
    hit = longest_repeat(seq, 8)
    assert hit is not None
    assert hit.motif == "AAACGTACGTAAA"
    assert seq[hit.first : hit.first + len(hit)] == hit.motif
    assert seq[hit.second : hit.second + len(hit)] == hit.motif


def test_no_repeat_returns_none():
    assert longest_repeat("ACGTACGT", 14) is None


def test_min_length_is_inclusive():
    seq = "ACGTACGT" + "TTTTTTTTTT" + "ACGTACGT"
    assert longest_repeat(seq, 8) is not None
    assert longest_repeat(seq, 9) is None or len(longest_repeat(seq, 9)) >= 9


def test_repeat_across_the_origin_is_found_only_when_circular():
    """One copy whole, the other split between the end of the plasmid and the start.

    A linear search cannot see the split copy, so it reports nothing; a plasmid really
    does carry the repeat, and a circular search finds it.
    """
    motif = "GGATCAGTTCACGA"  # 14 bp, not self-similar
    seq = motif[7:] + "GGGGCCCCGGGGCCCCAAAA" + motif + "TTTTAAAATTTT" + motif[:7]
    assert longest_repeat(seq, 14, circular=False) is None
    hit = longest_repeat(seq, 14, circular=True)
    assert hit is not None and len(hit) >= 14
    assert motif in (hit.motif, seq[hit.first : hit.first + len(hit)])


def test_a_sequence_is_not_its_own_circular_repeat():
    """Every motif reappears one turn later; that is the same occurrence, not a repeat."""
    assert longest_repeat("ACGTTGCAATCGGTACCATG", 6, circular=True) is None


def test_break_repeats_removes_the_repeat_and_keeps_the_protein():
    protein = "GGSDKGGSDKGGSDK"
    doubled = codon_optimize(protein, USAGE, argmax=True).dna
    seq = "AAAA" + doubled + "AAAA"
    assert longest_repeat(seq, 14) is not None
    result = break_repeats(seq, [(4, 4 + len(doubled))], USAGE, rng=random.Random(4), max_length=13)
    assert result.unfixable == []
    assert longest_repeat(result.sequence, 14) is None
    assert result.edits > 0
    assert translate(result.sequence[4 : 4 + len(doubled)]) == protein


def test_break_repeats_leaves_protected_spans_alone():
    protein = "GGSDKGGSDKGGSDK"
    doubled = codon_optimize(protein, USAGE, argmax=True).dna
    flank = "TTTTTTTT"
    seq = flank + doubled + flank
    result = break_repeats(
        seq, [(len(flank), len(flank) + len(doubled))], USAGE, rng=random.Random(5), max_length=13
    )
    assert result.sequence[: len(flank)] == flank
    assert result.sequence[-len(flank) :] == flank


def test_break_repeats_reports_a_repeat_it_cannot_touch():
    """Both copies outside the editable spans is a design problem, not an exception."""
    motif = codon_optimize("GGSDKGGSDK", USAGE, argmax=True).dna  # 30 bp
    seq = motif + "AAAACCCC" + motif
    result = break_repeats(seq, [], USAGE, rng=random.Random(6), max_length=13)
    assert result.sequence == seq
    assert result.unfixable and len(result.unfixable[0]) > 13


def test_break_repeats_respects_forbidden_motifs():
    protein = "GGSDKGGSDKGGSDK"
    doubled = codon_optimize(protein, USAGE, argmax=True).dna
    seq = "AAAA" + doubled + "AAAA"
    result = break_repeats(
        seq,
        [(4, 4 + len(doubled))],
        USAGE,
        rng=random.Random(7),
        max_length=13,
        forbidden=["GGTCTC"],
    )
    assert site_free(result.sequence, ["GGTCTC"])
    assert translate(result.sequence[4 : 4 + len(doubled)]) == protein


def test_codon_usage_rejects_an_empty_table():
    with pytest.raises(ValueError, match="empty"):
        CodonUsage({})


def test_synonyms_exclude_the_codon_itself():
    assert "GGC" not in USAGE.synonyms("GGC")
    assert set(USAGE.synonyms("GGC")) == {"GGT", "GGA", "GGG"}
    assert str(Seq("GGT").translate()) == "G"


def test_codon_usage_from_csv(tmp_path):
    """The usual entry point: a codon usage table downloaded as CSV."""
    path = tmp_path / "usage.csv"
    path.write_text(
        "codon,amino_acid,fraction\n"
        "GGU,G,0.34\n"  # RNA spelling, as several published tables use
        "GGC,G,0.40\n"
        "GAU,D,0.63\n"
        "GAC,D,0.37\n"
    )
    usage = CodonUsage.from_csv(path)
    assert usage.amino_acids() == ["D", "G"]
    assert usage.best("G") == "GGC", "U must be read as T"
    assert usage.codons("D") == ["GAT", "GAC"]
    assert usage.fraction("GGT") == 0.34
    assert codon_optimize("GD", usage, argmax=True).dna == "GGCGAT"
