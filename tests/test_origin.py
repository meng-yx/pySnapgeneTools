"""Moving a circular map's origin.

The thing worth testing is not the sequence rotation -- that is one slice -- but that
every feature still covers the same bases afterwards, including the one that ends up
straddling the new origin, and that the result survives a trip through the writer.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from pysnapgene import Q_COLOR, feature_start, read_snapgene, set_origin, write_snapgene


def _bases(record, feature):
    """The bases a feature covers, read off the record it belongs to."""
    return str(feature.extract(record.seq))


def _by_label(record, label):
    return next(
        f for f in record.features if str(f.qualifiers.get("label", [""])[0]) == label
    )


def test_sequence_is_rotated(record):
    original = str(record.seq)
    set_origin(record, 40)
    assert str(record.seq) == original[40:] + original[:40]
    assert len(record.seq) == len(original)


def test_every_feature_still_covers_the_same_bases(record):
    before = {
        str(f.qualifiers["label"][0]): _bases(record, f)
        for f in record.features
        if "label" in f.qualifiers
    }
    set_origin(record, 40)
    after = {
        str(f.qualifiers["label"][0]): _bases(record, f)
        for f in record.features
        if "label" in f.qualifiers
    }
    assert after == before


def test_a_feature_straddling_the_new_origin_becomes_two_parts(record):
    # GFP_frag is 29-68; an origin at 40 cuts it in two
    set_origin(record, 40)
    gfp = _by_label(record, "GFP_frag")
    assert isinstance(gfp.location, CompoundLocation)
    assert [(int(p.start), int(p.end)) for p in gfp.location.parts] == [
        (len(record.seq) - 11, len(record.seq)),
        (0, 28),
    ]
    assert len(gfp) == 39


def test_a_reverse_feature_straddling_the_origin_keeps_its_bases(record):
    before = _bases(record, _by_label(record, "spacer_rev"))
    set_origin(record, 80)  # spacer_rev is 68-98 on the minus strand
    moved = _by_label(record, "spacer_rev")
    assert isinstance(moved.location, CompoundLocation)
    assert moved.location.strand == -1
    assert _bases(record, moved) == before


def test_rotating_back_restores_the_original(record):
    original = str(record.seq)
    spans = {
        str(f.qualifiers["label"][0]): (int(f.location.start), int(f.location.end))
        for f in record.features
        if "label" in f.qualifiers
    }
    set_origin(record, 40)
    set_origin(record, len(original) - 40)
    assert str(record.seq) == original
    assert spans == {
        str(f.qualifiers["label"][0]): (int(f.location.start), int(f.location.end))
        for f in record.features
        if "label" in f.qualifiers
    }


def test_an_already_wrapping_feature_can_become_contiguous(record):
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [
                    SimpleLocation(len(record.seq) - 10, len(record.seq), strand=1),
                    SimpleLocation(0, 6, strand=1),
                ]
            ),
            type="misc_feature",
            qualifiers={"label": ["wraps"], Q_COLOR: ["#ff0000"]},
        )
    )
    before = _bases(record, _by_label(record, "wraps"))
    set_origin(record, len(record.seq) - 10)
    moved = _by_label(record, "wraps")
    assert isinstance(moved.location, SimpleLocation)
    assert (int(moved.location.start), int(moved.location.end)) == (0, 16)
    assert _bases(record, moved) == before


def test_origin_can_be_named_by_a_feature_label(record):
    set_origin(record, label="His6")
    assert str(record.seq).startswith("ATGCATCACCATCAC")
    assert int(_by_label(record, "His6").location.start) == 0


def test_feature_start_refuses_an_ambiguous_label(record):
    record.features.append(
        SeqFeature(
            SimpleLocation(0, 10, strand=1),
            type="misc_feature",
            qualifiers={"label": ["His6"]},
        )
    )
    with pytest.raises(KeyError, match="2 features"):
        feature_start(record, "His6")
    with pytest.raises(KeyError, match="no feature labelled"):
        feature_start(record, "absent")


def test_a_linear_molecule_is_refused(record):
    record.annotations["topology"] = "linear"
    with pytest.raises(ValueError, match="circular"):
        set_origin(record, 40)


def test_the_arguments_are_exclusive(record):
    with pytest.raises(TypeError):
        set_origin(record)
    with pytest.raises(TypeError):
        set_origin(record, 40, label="His6")
    with pytest.raises(ValueError, match="outside"):
        set_origin(record, len(record.seq))


def test_zero_is_a_no_op(record):
    before = [(int(f.location.start), int(f.location.end)) for f in record.features]
    set_origin(record, 0)
    assert [(int(f.location.start), int(f.location.end)) for f in record.features] == before


def test_the_rotated_map_round_trips_through_the_writer(record, tmp_path):
    set_origin(record, 40)
    bases = {
        str(f.qualifiers["label"][0]): _bases(record, f)
        for f in record.features
        if "label" in f.qualifiers
    }
    path = tmp_path / "rotated.dna"
    write_snapgene(record, path)
    back = read_snapgene(path)
    assert str(back.seq) == str(record.seq)
    for label, expected in bases.items():
        assert _bases(back, _by_label(back, label)) == expected, label


def test_a_full_length_feature_survives(record):
    length = len(record.seq)
    record.features.append(
        SeqFeature(
            SimpleLocation(0, length, strand=1),
            type="source",
            qualifiers={"label": ["whole"]},
        )
    )
    set_origin(record, 40)
    whole = _by_label(record, "whole")
    assert len(whole) == length


def test_a_stale_primer_binding_site_does_not_travel(record):
    primer = _by_label(record, "seq_fwd")
    primer.qualifiers["snapgene_site_range"] = ["123-145"]
    primer.qualifiers["snapgene_site_xml"] = ["<BindingSite .../>"]
    set_origin(record, 40)
    moved = _by_label(record, "seq_fwd")
    assert "snapgene_site_range" not in moved.qualifiers
    assert "snapgene_site_xml" not in moved.qualifiers


def test_rotation_is_what_the_restriction_map_says(record):
    """A landmark found by sequence search lands where set_origin says it does."""
    motif = "TTAATTAA"
    original = str(record.seq)
    at = original.find(motif)
    assert at > 0
    rec = SeqRecord(
        Seq(original),
        id="pTest001",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    set_origin(rec, at)
    assert str(rec.seq).startswith(motif)
    assert str(rec.seq).count(motif) == original.count(motif)
