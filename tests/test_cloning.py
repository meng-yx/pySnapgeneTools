"""Type IIS digestion and Golden Gate assembly."""

from __future__ import annotations

import pytest
from Bio.Restriction import BsaI, EcoRI, SapI
from Bio.Seq import Seq

from pysnapgene.cloning import digest, excise, golden_gate, translate_orf

# A toy destination vector shaped like a real one: ATG, a short N-terminal tag, then the
# Golden Gate cassette, then a tail and a stop. Both BsaI sites point *inwards*, so the
# digest carries them away with the stuffer and the assembled product has none left --
# that orientation is the whole trick of a Golden Gate destination vector.
#
#   ATG tag  HEAD | AGGA n GAGACC  stuffer  GGTCTC n | TTCC  tail  stop  filler
#                  <------------- released ------->
TAG = "ATGGGTGGTAGC"  # M G G S
HEAD = "GGTTC"  # puts the AGGA junction two bases into a codon, as pET29b-NanoBiT does
STUFFER = "GGATGCATGCATGCATGCATGCATGCAAA"
TAILFRAME = "GGTGGTTCCGG"
FILLER = "CCCGGGTTTAAACCCGGGTTTAAA"
VECTOR = (
    TAG
    + HEAD
    + "AGGA"
    + "G"
    + "GAGACC"
    + STUFFER
    + "GGTCTC"
    + "A"
    + "TTCC"
    + TAILFRAME
    + "TGA"
    + FILLER
)

# An insert with the matching overhangs and its own sites pointing outwards, so they are
# discarded instead of ending up in the product
BODY = "AGCATGAAAGGCTGGGATCCGAAA"
INSERT = "CCCC" + "GGTCTC" + "A" + "AGGA" + BODY + "TTCC" + "G" + "GAGACC" + "GGGG"

#: the fusion the assembly is meant to produce
EXPECTED_PROTEIN = "MGGSGSGSMKGWDPKFRWFR"


def test_digest_reports_the_overhang_at_each_end():
    fragments = digest(VECTOR, BsaI, circular=True)
    assert len(fragments) == 2
    assert {(f.left, f.right) for f in fragments} == {("AGGA", "TTCC"), ("TTCC", "AGGA")}


def test_circular_digest_accounts_for_every_base_once():
    fragments = digest(VECTOR, BsaI, circular=True)
    assert sum(len(f) for f in fragments) == len(VECTOR)
    # the fragments are the vector, rotated to the first cut
    rejoined = "".join(f.seq for f in fragments)
    assert rejoined in VECTOR + VECTOR
    assert len(rejoined) == len(VECTOR)


def test_linear_digest_has_blunt_outer_ends():
    fragments = digest(INSERT, BsaI, circular=False)
    assert len(fragments) == 3
    assert fragments[0].left == "" and fragments[-1].right == ""
    assert "".join(f.seq for f in fragments) == INSERT


def test_excise_separates_backbone_from_stuffer():
    backbone, stuffer = excise(VECTOR, BsaI)
    assert len(backbone) > len(stuffer)
    assert STUFFER in stuffer.seq
    assert (stuffer.left, stuffer.right) == ("AGGA", "TTCC")
    assert (backbone.left, backbone.right) == ("TTCC", "AGGA")
    # both recognition sites leave with the stuffer
    assert len(BsaI.search(Seq(backbone.seq), linear=True)) == 0


def test_excise_rejects_a_vector_not_cut_twice():
    with pytest.raises(ValueError, match="cut twice"):
        excise(TAG + "GGTCTC" + "A" + "AGGA" + STUFFER + FILLER, BsaI)


def test_golden_gate_replaces_the_stuffer_with_the_insert():
    product = golden_gate(VECTOR, INSERT, BsaI)
    assert BODY in product
    assert STUFFER not in product
    assert len(product) == len(VECTOR) - 47 + 28  # stuffer out, insert in
    assert len(BsaI.search(Seq(product), linear=False)) == 0, "the sites must be consumed"


def test_the_product_reads_as_one_fusion_in_frame():
    product = golden_gate(VECTOR, INSERT, BsaI)
    start = product.find(TAG)
    assert start >= 0
    assert translate_orf(product, start) == EXPECTED_PROTEIN


def test_the_junctions_contribute_residues_of_their_own():
    """Why simulating beats counting: the overhangs add residues at both ends.

    Between the tag (``MGGS``) and the insert's own ``M`` the product reads ``GSGS`` -- four
    residues that exist in neither the vector's tag nor the insert's gene.
    """
    product = golden_gate(VECTOR, INSERT, BsaI)
    protein = translate_orf(product, product.find(TAG))
    tag, gene = "MGGS", "MKGWDPK"
    between = protein[len(tag) : protein.index(gene)]
    assert between == "GSGS"
    assert len(between) == 4


def test_overhang_mismatch_is_reported_rather_than_silently_assembled():
    """An insert ordered against the wrong overhangs looks fine until it will not go in."""
    wrong = "CCCC" + "GGTCTC" + "A" + "TGCA" + BODY + "ACGT" + "G" + "GAGACC" + "GGGG"
    with pytest.raises(ValueError, match="overhangs"):
        golden_gate(VECTOR, wrong, BsaI)


def test_blunt_enzyme_is_refused():
    with pytest.raises(ValueError, match="blunt|does not cut"):
        digest(VECTOR, EcoRI, circular=True)


def test_uncut_sequence_is_reported():
    with pytest.raises(ValueError, match="does not cut"):
        digest("AAAACCCCGGGGTTTT", BsaI, circular=True)


def test_another_type_iis_enzyme_uses_its_own_overhang_length():
    """SapI leaves three bases, not four, so the length must come from the enzyme."""
    vector = (
        TAG
        + HEAD
        + "GGA"
        + "G"
        + "GAAGAGC"
        + STUFFER
        + "GCTCTTC"
        + "A"
        + "TCC"
        + TAILFRAME
        + "TGA"
        + FILLER
    )
    backbone, stuffer = excise(vector, SapI)
    assert (stuffer.left, stuffer.right) == ("GGA", "TCC")
    assert (backbone.left, backbone.right) == ("TCC", "GGA")
    assert len(stuffer.left) == 3


def test_translate_orf_wraps_the_origin():
    """An ORF in a plasmid frequently starts near the end and finishes after position 0."""
    orf = "ATG" + "AAAGGC" * 4 + "TAA"
    plasmid = orf[12:] + "CCCCGGGGCCCCGGGG" + orf[:12]
    start = len(plasmid) - 12
    assert translate_orf(plasmid, start, circular=True) == "MKGKGKGKG"


def test_translate_orf_without_a_stop_raises():
    with pytest.raises(ValueError, match="no stop codon"):
        translate_orf("ATG" + "AAA" * 20, 0, circular=False)


def test_fragment_start_locates_the_fragment_in_the_sequence():
    """``start`` must be the fragment's own offset, origin-wrapping included.

    Biopython returns a circular digest rotated so the first fragment begins at the *last*
    cut, so a start taken straight from the cut list in order is off by one position.
    """
    for fragment in digest(VECTOR, BsaI, circular=True):
        doubled = VECTOR + VECTOR
        assert doubled[fragment.start : fragment.start + len(fragment)] == fragment.seq
        assert 0 <= fragment.start < len(VECTOR)


def test_fragment_start_in_a_linear_digest():
    offset = 0
    for fragment in digest(INSERT, BsaI, circular=False):
        assert fragment.start == offset
        assert INSERT[fragment.start : fragment.start + len(fragment)] == fragment.seq
        offset += len(fragment)


def test_fragment_start_with_three_cut_sites():
    """More than two sites, to catch an off-by-one that two sites can hide."""
    triple = "GGTCTCA" + "AAAACCCC" * 3 + "GGTCTCA" + "GGGGTTTT" * 3 + "GGTCTCA" + "CCCCAAAA" * 4
    fragments = digest(triple, BsaI, circular=True)
    assert len(fragments) == 3
    doubled = triple + triple
    for fragment in fragments:
        assert doubled[fragment.start : fragment.start + len(fragment)] == fragment.seq
