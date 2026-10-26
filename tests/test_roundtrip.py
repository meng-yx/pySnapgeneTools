"""Round-trip behaviour: what goes in must come back out."""

from __future__ import annotations

import pytest
from Bio.SeqFeature import SeqFeature, SimpleLocation

from pysnapgene import (
    DEFAULT_COLOR,
    FLAG_CIRCULAR,
    FLAG_DOUBLE_STRANDED,
    Q_COLOR,
    Q_TRANSLATED,
    SnapGeneError,
    read_snapgene,
    to_bytes,
    write_snapgene,
)


def _reload(record, tmp_path, **kwargs):
    path = tmp_path / "out.dna"
    write_snapgene(record, path, **kwargs)
    return read_snapgene(path)


def test_sequence_and_topology_survive(record, tmp_path):
    back = _reload(record, tmp_path)
    assert str(back.seq) == str(record.seq)
    assert back.annotations["topology"] == "circular"


def test_molecule_is_double_stranded_by_default(record, tmp_path):
    back = _reload(record, tmp_path)
    flags = back.annotations["snapgene_flags"]
    assert flags & FLAG_DOUBLE_STRANDED, "SnapGene would show this as single stranded"
    assert flags & FLAG_CIRCULAR
    assert back.annotations["double_stranded"] is True


def test_linear_topology(record, tmp_path):
    back = _reload(record, tmp_path, circular=False)
    assert back.annotations["topology"] == "linear"
    assert not back.annotations["snapgene_flags"] & FLAG_CIRCULAR
    # still double stranded
    assert back.annotations["snapgene_flags"] & FLAG_DOUBLE_STRANDED


def test_feature_positions_and_strands(record, tmp_path):
    back = _reload(record, tmp_path)
    expected = {
        (int(f.location.start), int(f.location.end), f.location.strand) for f in record.features
    }
    got = {(int(f.location.start), int(f.location.end), f.location.strand) for f in back.features}
    assert got == expected


def test_colours_are_preserved(record, tmp_path):
    back = _reload(record, tmp_path)
    by_label = {f.qualifiers["label"][0]: f for f in back.features}
    assert by_label["promoter"].qualifiers[Q_COLOR] == ["#ffffff"]
    assert by_label["GFP_frag"].qualifiers[Q_COLOR] == ["#75c6a9"]
    assert by_label["His6"].qualifiers[Q_COLOR] == ["#c7b0e3"]
    assert by_label["spacer_rev"].qualifiers[Q_COLOR] == ["#d6d6d6"]


def test_uncoloured_feature_gets_the_default(record, tmp_path):
    record.features.append(
        SeqFeature(
            SimpleLocation(0, 12, strand=1),
            type="misc_feature",
            qualifiers={"label": ["no_colour"]},
        )
    )
    back = _reload(record, tmp_path)
    plain = [f for f in back.features if f.qualifiers["label"][0] == "no_colour"][0]
    assert plain.qualifiers[Q_COLOR] == [DEFAULT_COLOR]


def test_custom_default_colour(record, tmp_path):
    record.features.append(
        SeqFeature(
            SimpleLocation(0, 12, strand=1),
            type="misc_feature",
            qualifiers={"label": ["no_colour"]},
        )
    )
    back = _reload(record, tmp_path, default_color="#123456")
    plain = [f for f in back.features if f.qualifiers["label"][0] == "no_colour"][0]
    assert plain.qualifiers[Q_COLOR] == ["#123456"]


def test_clean_orf_is_marked_translated(record, tmp_path):
    back = _reload(record, tmp_path)
    by_label = {f.qualifiers["label"][0]: f for f in back.features}
    assert by_label["GFP_frag"].qualifiers.get(Q_TRANSLATED) == ["1"]
    assert by_label["His6"].qualifiers.get(Q_TRANSLATED) == ["1"]
    assert "translation" in by_label["GFP_frag"].qualifiers


def test_translation_matches_the_sequence(record, tmp_path):
    back = _reload(record, tmp_path)
    gfp = [f for f in back.features if f.qualifiers["label"][0] == "GFP_frag"][0]
    expected = str(record.seq[29:68].translate())
    assert gfp.qualifiers["translation"] == [expected]


def test_non_coding_feature_is_not_translated(record, tmp_path):
    back = _reload(record, tmp_path)
    spacer = [f for f in back.features if f.qualifiers["label"][0] == "spacer_rev"][0]
    assert Q_TRANSLATED not in spacer.qualifiers


def test_cds_that_is_not_a_clean_orf_is_demoted(record, tmp_path):
    """A span that is not a whole number of codons must not claim a protein."""
    record.features.append(
        SeqFeature(
            SimpleLocation(0, 20, strand=1),  # 20 nt: not divisible by three
            type="CDS",
            qualifiers={"label": ["ragged"], Q_COLOR: ["#ff0000"]},
        )
    )
    back = _reload(record, tmp_path)
    ragged = [f for f in back.features if f.qualifiers["label"][0] == "ragged"][0]
    assert Q_TRANSLATED not in ragged.qualifiers
    assert "translation" not in ragged.qualifiers
    assert ragged.qualifiers[Q_COLOR] == ["#ff0000"]  # colour is still kept


def test_cds_with_internal_stop_is_demoted(record, tmp_path):
    record.features.append(
        SeqFeature(
            SimpleLocation(29, 71, strand=1),  # runs past the ORF's stop codon
            type="CDS",
            qualifiers={"label": ["readthrough"]},
        )
    )
    back = _reload(record, tmp_path)
    bad = [f for f in back.features if f.qualifiers["label"][0] == "readthrough"][0]
    assert Q_TRANSLATED not in bad.qualifiers


def test_primers_land_in_the_primer_packet(record, tmp_path):
    back = _reload(record, tmp_path)
    primers = [f for f in back.features if f.type == "primer_bind"]
    assert len(primers) == 1
    assert primers[0].qualifiers["label"] == ["seq_fwd"]
    assert int(primers[0].location.start) == 122
    assert int(primers[0].location.end) == 145


def test_short_primers_are_dropped(record, tmp_path):
    record.features.append(
        SeqFeature(
            SimpleLocation(10, 18, strand=1),  # 8 nt, below SnapGene's threshold
            type="primer_bind",
            qualifiers={"label": ["too_short"]},
        )
    )
    back = _reload(record, tmp_path)
    labels = {f.qualifiers["label"][0] for f in back.features}
    assert "too_short" not in labels
    assert "seq_fwd" in labels


def test_duplicate_features_are_collapsed(record, tmp_path):
    record.features.append(record.features[0])
    back = _reload(record, tmp_path)
    promoters = [f for f in back.features if f.qualifiers["label"][0] == "promoter"]
    assert len(promoters) == 1


def test_reverse_strand_cds_translates_from_the_right_strand(tmp_path, record):
    # a reverse-strand ORF: reverse complement of ATG...TAA
    from Bio.Seq import Seq as _Seq

    orf = "ATGAAACCCGGGTTTTAA"
    record.seq = record.seq + _Seq(orf).reverse_complement()
    start = len(record.seq) - len(orf)
    record.features.append(
        SeqFeature(
            SimpleLocation(start, start + len(orf), strand=-1),
            type="CDS",
            qualifiers={"label": ["rev_orf"]},
        )
    )
    back = _reload(record, tmp_path)
    rev = [f for f in back.features if f.qualifiers["label"][0] == "rev_orf"][0]
    assert rev.qualifiers.get(Q_TRANSLATED) == ["1"]
    assert rev.qualifiers["translation"] == [str(_Seq(orf).translate())]


def test_notes_carry_identity(record, tmp_path):
    back = _reload(record, tmp_path)
    assert back.id == "pTest001"
    assert "synthetic test plasmid" in back.description


def test_unknown_packets_are_preserved(record, tmp_path):
    from pysnapgene import Packet

    extra = Packet(0x03, b"pretend enzyme set")
    path = tmp_path / "with_extra.dna"
    write_snapgene(record, path, preserved_packets=[extra])
    back = read_snapgene(path)
    kept = back.annotations["snapgene_unknown_packets"]
    assert any(p.type == 0x03 and p.data == b"pretend enzyme set" for p in kept)


def test_unknown_packets_survive_a_second_round_trip(record, tmp_path):
    from pysnapgene import Packet

    first = tmp_path / "a.dna"
    write_snapgene(record, first, preserved_packets=[Packet(0x03, b"xyz")])
    middle = read_snapgene(first)
    second = tmp_path / "b.dna"
    write_snapgene(middle, second)
    back = read_snapgene(second)
    assert any(
        p.type == 0x03 and p.data == b"xyz" for p in back.annotations["snapgene_unknown_packets"]
    )


def test_write_to_file_object(record, tmp_path):
    path = tmp_path / "handle.dna"
    with path.open("wb") as handle:
        assert write_snapgene(record, handle) is None
    assert read_snapgene(path).seq == record.seq


def test_read_from_bytes(record):
    payload = to_bytes(record)
    assert str(read_snapgene(payload).seq) == str(record.seq)


def test_empty_sequence_rejected(record):
    record.seq = record.seq[:0]
    with pytest.raises(ValueError):
        to_bytes(record)


def test_not_a_snapgene_file():
    with pytest.raises(SnapGeneError):
        read_snapgene(b"\x09\x00\x00\x00\x0eNOTSNAPG\x00\x01\x00\x0f\x00\x14")


def test_missing_cookie_rejected(record):
    payload = to_bytes(record)
    # strip the leading cookie packet (5 byte header + 14 byte payload)
    with pytest.raises(SnapGeneError):
        read_snapgene(payload[19:])


def test_real_file_round_trip(real_dna_path, tmp_path):
    """Optional: only runs when PYSNAPGENE_TEST_DNA points at a real file."""
    original = read_snapgene(real_dna_path)
    path = tmp_path / "real.dna"
    write_snapgene(original, path)
    back = read_snapgene(path)

    # an unedited map comes back byte for byte
    assert path.read_bytes() == real_dna_path.read_bytes()

    assert str(back.seq) == str(original.seq)
    assert back.annotations["topology"] == original.annotations["topology"]

    def key(feature):
        return (
            int(feature.location.start),
            int(feature.location.end),
            feature.location.strand,
            feature.qualifiers.get("label", ["?"])[0],
            feature.qualifiers.get("snapgene_color", ["?"])[0],
        )

    # short primer fragments are dropped on purpose, so compare what survives
    kept = {
        key(f)
        for f in original.features
        if not (f.type == "primer_bind" and int(f.location.end) - int(f.location.start) < 15)
    }
    assert kept == {key(f) for f in back.features}

    def shown(record):
        return {
            key(f): (
                f.type,
                f.qualifiers.get("translation"),
                f.qualifiers.get("snapgene_segment_names"),
            )
            for f in record.features
        }

    # same types, same proteins, same segment names
    after = shown(back)
    assert {k: v for k, v in shown(original).items() if k in after} == after
