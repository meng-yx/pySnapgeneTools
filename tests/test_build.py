"""Protein to DNA: the deterministic builder.

What has to hold for the output to be orderable: the protein is untouched, frozen parts
are untouched, forbidden motifs are absent from what was encoded, two copies of one element
are spelled differently, and the same inputs give the same sequence every time.
"""

from __future__ import annotations

import random

import pytest
from Bio.Seq import Seq

from pysnapgene.design import (
    TWIST_GENE_FRAGMENT,
    CodonUsage,
    build_dna,
    complexity_report,
    exact_repeats,
    fix_complexity,
    frozen,
    identity,
    identity_floor,
    protein,
    repeat_density,
    score,
    translate,
)

# Real E. coli-ish usage for the amino acids the test proteins use, with enough synonyms
# to leave the builder a choice.
USAGE = CodonUsage(
    {
        "GGT": ("G", 0.34), "GGC": ("G", 0.40), "GGA": ("G", 0.11), "GGG": ("G", 0.15),
        "AGT": ("S", 0.15), "AGC": ("S", 0.28), "TCT": ("S", 0.15), "TCC": ("S", 0.15),
        "TCA": ("S", 0.12), "TCG": ("S", 0.15),
        "GAT": ("D", 0.63), "GAC": ("D", 0.37),
        "AAA": ("K", 0.74), "AAG": ("K", 0.26),
        "CTG": ("L", 0.50), "CTT": ("L", 0.10), "CTC": ("L", 0.10), "CTA": ("L", 0.04),
        "TTA": ("L", 0.13), "TTG": ("L", 0.13),
        "GAA": ("E", 0.68), "GAG": ("E", 0.32),
        "GTG": ("V", 0.37), "GTT": ("V", 0.26), "GTC": ("V", 0.22), "GTA": ("V", 0.15),
        "ACC": ("T", 0.44), "ACT": ("T", 0.17), "ACA": ("T", 0.13), "ACG": ("T", 0.27),
        "ATG": ("M", 1.00),
        "TGG": ("W", 1.00),
        "CAT": ("H", 0.57), "CAC": ("H", 0.43),
        "CCG": ("P", 0.52), "CCA": ("P", 0.19), "CCT": ("P", 0.16), "CCC": ("P", 0.12),
        "TTT": ("F", 0.57), "TTC": ("F", 0.43),
        "CGT": ("R", 0.38), "CGC": ("R", 0.40), "CGG": ("R", 0.10), "CGA": ("R", 0.06),
        "AAT": ("N", 0.45), "AAC": ("N", 0.55),
        "CAA": ("Q", 0.35), "CAG": ("Q", 0.65),
        "ATT": ("I", 0.51), "ATC": ("I", 0.42), "ATA": ("I", 0.07),
        "GCG": ("A", 0.36), "GCC": ("A", 0.27), "GCA": ("A", 0.21), "GCT": ("A", 0.16),
        "TAT": ("Y", 0.57), "TAC": ("Y", 0.43),
        "TGT": ("C", 0.45), "TGC": ("C", 0.55),
    }
)

PEPTIDE = (
    "MKSTVLGDWEQTAAYNLDQVLEQGGVSSLLQNLAVSVTPIQRIVRSGENALKIDIHVIIPYEGLSADQMAQIEEVFK"
    "VVYPVDDHHFKVILPYGTLVIDGVTPNMLNYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLITPDGSMLFRVTIN"
)
SHORT = "MKSTVLGDWEQTAAYNLDQ"


def _random_dna(length: int, seed: int) -> str:
    return "".join(random.Random(seed).choices("ATGC", k=length))


def test_the_protein_is_what_was_asked_for():
    design = build_dna([protein("p", PEPTIDE)], USAGE)
    assert translate(design.sequence) == PEPTIDE


def test_the_same_inputs_give_the_same_sequence():
    plan = [protein("a", PEPTIDE), frozen("mid", "ATGCATGC"), protein("b", PEPTIDE,
                                                                     unlike=["a"])]
    first = build_dna(plan, USAGE)
    second = build_dna(plan, USAGE)
    assert first.sequence == second.sequence
    # and no hidden global state: a different random stream must not matter
    random.seed(1234)
    assert build_dna(plan, USAGE).sequence == first.sequence


def test_frozen_parts_are_kept_verbatim():
    cassette = _random_dna(90, 1)
    design = build_dna(
        [protein("head", SHORT), frozen("cassette", cassette), protein("tail", SHORT)],
        USAGE,
    )
    assert design.part("cassette") == cassette
    assert cassette in design.sequence


def test_a_frozen_part_may_contain_a_forbidden_motif():
    """A ccdB cassette carries its own BsaI sites; the parts beside it still get encoded."""
    cassette = "AAAA" + "GGTCTC" + _random_dna(60, 2) + "GAGACC" + "AAAA"
    design = build_dna(
        [protein("head", SHORT), frozen("cassette", cassette), protein("tail", SHORT)],
        USAGE,
        forbidden=["GGTCTC"],
    )
    assert design.part("cassette") == cassette
    assert design.sequence.count("GGTCTC") == 1, "no new site was created"
    assert design.sequence.count("GAGACC") == 1


def test_forbidden_motifs_are_kept_out_of_encoded_parts():
    design = build_dna([protein("p", PEPTIDE * 2)], USAGE, forbidden=["GGTCTC", "GAATTC"])
    for motif in ("GGTCTC", "GAATTC"):
        assert motif not in design.sequence
        assert str(Seq(motif).reverse_complement()) not in design.sequence


def test_unlike_makes_the_second_copy_different():
    plain = build_dna([protein("a", PEPTIDE), protein("b", PEPTIDE)], USAGE)
    apart = build_dna(
        [protein("a", PEPTIDE), protein("b", PEPTIDE, unlike=["a"])], USAGE
    )
    same = identity(plain.part("a"), plain.part("b"))
    different = identity(apart.part("a"), apart.part("b"))
    # even without unlike the encoder avoids repeating what it has already written, so the
    # second copy is not verbatim; unlike is what drives the identity down to the floor
    assert 0.80 < same < 0.95
    assert different < 0.75
    assert different < same - 0.10
    assert translate(apart.part("b")) == PEPTIDE


def test_unlike_lowers_the_repeat_density_too():
    plain = build_dna([protein("a", PEPTIDE), protein("b", PEPTIDE)], USAGE, repair=False)
    apart = build_dna([protein("a", PEPTIDE), protein("b", PEPTIDE, unlike=["a"])], USAGE,
                      repair=False)
    assert repeat_density(apart.sequence) < repeat_density(plain.sequence)
    assert max((p.length for p in exact_repeats(apart.sequence, 9)), default=0) < 20


def test_unlike_must_name_an_earlier_part():
    with pytest.raises(ValueError, match="not built before it"):
        build_dna([protein("b", SHORT, unlike=["a"]), protein("a", SHORT)], USAGE)


def test_duplicate_names_are_refused():
    with pytest.raises(ValueError, match="share a name"):
        build_dna([protein("a", SHORT), protein("a", SHORT)], USAGE)


def test_an_editable_dna_part_must_be_in_frame():
    with pytest.raises(ValueError, match="whole number of codons"):
        frozen("x", "ATGCA", editable=True)
    assert frozen("x", "ATGCA").dna == "ATGCA", "a frozen part needs no frame"


def test_an_empty_or_broken_protein_is_refused():
    with pytest.raises(ValueError, match="empty"):
        protein("x", "")
    with pytest.raises(ValueError, match="internal stop"):
        protein("x", "MK*MK")


def test_spans_locate_every_part():
    design = build_dna(
        [protein("a", SHORT), frozen("mid", "ATGCATGCAT"), protein("b", SHORT)], USAGE
    )
    for name, (start, end) in design.spans.items():
        assert design.sequence[start:end] == design.part(name)
    assert design.spans["a"][1] == design.spans["mid"][0]
    assert design.spans["mid"][1] == design.spans["b"][0]


def test_repair_lowers_the_score_without_changing_the_protein():
    repetitive = build_dna([protein("p", PEPTIDE), protein("q", PEPTIDE)], USAGE,
                           repair=False)
    before = score(repetitive.sequence)
    assert before > 0
    fixed = fix_complexity(
        repetitive.sequence, [(0, len(repetitive.sequence))], USAGE
    )
    assert score(fixed.sequence) < before
    assert translate(fixed.sequence) == translate(repetitive.sequence)
    assert fixed.edits > 0


def test_repair_leaves_everything_outside_the_editable_spans_alone():
    cassette = _random_dna(99, 3)
    start = build_dna([protein("p", PEPTIDE), frozen("c", cassette),
                       protein("q", PEPTIDE)], USAGE, repair=False)
    at = start.spans["c"]
    editable = [(0, at[0]), (at[1], len(start.sequence))]
    fixed = fix_complexity(start.sequence, editable, USAGE)
    assert fixed.sequence[at[0]:at[1]] == cassette


def test_an_unfixable_issue_is_reported_not_raised():
    """Two identical frozen copies cannot be repaired; the report says so."""
    block = _random_dna(120, 4)
    sequence = block + _random_dna(300, 5) + block
    fixed = fix_complexity(sequence, [], USAGE)      # nothing editable at all
    assert fixed.sequence == sequence
    assert fixed.edits == 0
    assert any(issue.kind == "exact_repeat" for issue in fixed.issues)
    assert not fixed.clean


def test_a_build_comes_back_clean_and_says_what_it_changed():
    design = build_dna([protein("p", PEPTIDE * 2)], USAGE)
    assert design.clean, design.report()
    # a tandem duplication is the hard case; the encoder gets most of the way there and
    # the repair finishes it off with a handful of codons
    assert design.edits < 0.02 * len(design.sequence)
    assert design.rounds == design.edits or design.rounds >= 1


def test_report_mentions_what_is_left():
    design = build_dna([protein("p", SHORT)], USAGE)       # 57 bp: under the length floor
    assert "length" in design.report()
    assert not design.clean


def test_ignore_passes_through_to_the_rules():
    design = build_dna([protein("p", SHORT)], USAGE, ignore=["length"])
    assert design.clean


def test_rules_are_honoured():
    strict = TWIST_GENE_FRAGMENT._replace(max_repeat_density=0.0, min_length=0)
    design = build_dna([protein("p", PEPTIDE)], USAGE, rules=strict, rounds=50)
    assert repeat_density(design.sequence) <= repeat_density(
        build_dna([protein("p", PEPTIDE)], USAGE, rules=strict, repair=False).sequence
    )


def test_identity_measures_what_exact_repeats_cannot():
    a = _random_dna(300, 6)
    b = "".join(base if i % 10 else "A" for i, base in enumerate(a))
    assert 0.85 < identity(a, b) < 1.0
    assert identity(a, a) == 1.0
    assert identity("", "x") == 0.0


def test_a_real_design_clears_the_vendor_rules():
    """The shape that was quoted as complex: two copies of one protein round a cassette."""
    cassette = _random_dna(469, 7)
    design = build_dna(
        [
            frozen("arm", "AATTGGTCTCAAGGA"),
            protein("copy_a", PEPTIDE),
            frozen("ccdB", cassette),
            protein("copy_b", PEPTIDE, unlike=["copy_a"]),
        ],
        USAGE,
        forbidden=["GGTCTC"],
    )
    assert design.clean, design.report()
    assert identity(design.part("copy_a"), design.part("copy_b")) < 0.70
    assert complexity_report(design.sequence) == []


def test_identity_floor_is_a_floor():
    floor = identity_floor(PEPTIDE, USAGE)
    design = build_dna([protein("a", PEPTIDE), protein("b", PEPTIDE, unlike=["a"])], USAGE)
    measured = identity(design.part("a"), design.part("b"))
    assert 0.5 < floor < 0.75
    assert measured >= floor - 1e-9, "no spelling can beat the floor"
    assert measured < floor + 0.10, "unlike should land close to it"


def test_identity_floor_is_one_for_a_protein_with_no_choice():
    assert identity_floor("MW", USAGE) == 1.0
    assert identity_floor("", USAGE) == 0.0
