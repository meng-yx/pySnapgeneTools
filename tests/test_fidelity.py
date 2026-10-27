"""Details of real SnapGene files that a read/write cycle must not change.

Each test builds a small file by hand, shaped like what SnapGene writes, so the
suite still needs no real plasmid.
"""

from __future__ import annotations

import re
import struct

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from pysnapgene import Q_TRANSLATED, iter_packets, pack_packet, read_snapgene, to_bytes
from pysnapgene.writer.constants import Q_SEGMENT_NAMES

ORF = "ATG" + "AAA" * 5 + "GGC" * 5 + "TAA"  # 36 nt, a clean ORF
SEQ = "ccccc" * 2 + str(Seq(ORF).reverse_complement()).lower() + "ttttt" * 10  # 96 nt

NOTES = (
    "<Notes><UUID>abc</UUID><Type>Synthetic</Type>"
    "<Description>a test vector</Description>"
    '<Created UTC="1:2:3">2024.1.2</Created>'
    '<LastModified UTC="1:2:3">2024.1.2</LastModified>'
    '<References><Reference title="A paper" pubMedID="1"/></References>'
    "</Notes>"
)


def _file(features_xml: str, primers_xml: str = "", notes: str = NOTES) -> bytes:
    out = pack_packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 15, 20))
    out += pack_packet(0x00, b"\x03" + SEQ.encode())
    out += pack_packet(0x0A, f'<?xml version="1.0"?><Features>{features_xml}</Features>'.encode())
    if primers_xml:
        out += pack_packet(0x05, f'<?xml version="1.0"?><Primers>{primers_xml}</Primers>'.encode())
    out += pack_packet(0x06, notes.encode())
    return out


def _packet(data: bytes, kind: int) -> str:
    return next(p.data for p in iter_packets(data) if p.type == kind).decode()


def _segments(data: bytes) -> list[str]:
    return re.findall(r"<Segment [^>]*/>", _packet(data, 0x0A))


def test_split_minus_strand_cds_stays_translated():
    # the ORF sits reverse-complemented at 11-46, split in two like a signal peptide
    feature = (
        '<Feature name="split" type="CDS" directionality="2">'
        '<Segment range="11-31" color="#ccffcc" translated="1"/>'
        '<Segment name="signal" range="32-46" color="#ccffcc" translated="1"/>'
        "</Feature>"
    )
    record = read_snapgene(_file(feature))
    split = record.features[0]
    assert str(split.location.extract(Seq(SEQ.upper()))) == ORF

    back = read_snapgene(to_bytes(record))
    split = back.features[0]
    assert split.type == "CDS"
    assert split.qualifiers.get(Q_TRANSLATED) == ["1"]
    assert split.qualifiers["translation"] == [str(Seq(ORF).translate())]


def test_segments_keep_their_order_and_names():
    feature = (
        '<Feature name="promoter" type="promoter" directionality="2">'
        '<Segment name="-10" range="11-16" color="#ffffff"/>'
        '<Segment range="17-30" color="#ff0000"/>'
        '<Segment name="-35" range="31-36" color="#00ff00"/>'
        "</Feature>"
    )
    data = _file(feature)
    record = read_snapgene(data)
    assert record.features[0].qualifiers[Q_SEGMENT_NAMES] == ["-10", "", "-35"]

    out = to_bytes(record)
    assert _segments(out) == [
        '<Segment name="-10" range="11-16" color="#ffffff" type="standard"/>',
        '<Segment range="17-30" color="#ff0000" type="standard"/>',
        '<Segment name="-35" range="31-36" color="#00ff00" type="standard"/>',
    ]


def test_feature_across_the_origin_is_one_segment():
    for direction in ("1", "2"):
        feature = (
            f'<Feature name="ori" type="rep_origin" directionality="{direction}">'
            '<Segment range="91-6" color="#ffff00"/></Feature>'
        )
        record = read_snapgene(_file(feature))
        assert len(record.features[0].location) == 12
        out = to_bytes(record)
        assert _segments(out) == ['<Segment range="91-6" color="#ffff00" type="standard"/>'], (
            direction
        )


def test_minus_strand_origin_feature_reads_in_the_right_order():
    feature = (
        '<Feature name="wrap" type="misc_feature" directionality="2">'
        '<Segment range="91-6" color="#ffff00"/></Feature>'
    )
    record = read_snapgene(_file(feature))
    expected = str(Seq(SEQ[90:] + SEQ[:6]).reverse_complement())
    assert str(record.features[0].location.extract(Seq(SEQ))) == expected


def test_primer_bind_feature_stays_a_feature():
    feature = (
        '<Feature name="M13 rev" type="primer_bind" directionality="1">'
        '<Segment range="50-66" color="#a020f0"/></Feature>'
    )
    out = to_bytes(read_snapgene(_file(feature)))
    assert 'name="M13 rev"' in _packet(out, 0x0A)
    assert not any(p.type == 0x05 for p in iter_packets(out))


def test_primer_description_survives():
    primer = (
        '<Primer name="fwd" sequence="TTTTTTTTTTTTTTTTTT" description="a forward primer">'
        '<BindingSite location="60-77" boundStrand="0"/></Primer>'
    )
    back = read_snapgene(to_bytes(read_snapgene(_file("", primer))))
    fwd = [f for f in back.features if f.type == "primer_bind"][0]
    assert fwd.qualifiers["note"] == ["a forward primer"]


def test_notes_are_kept_whole():
    record = read_snapgene(_file(""))
    assert record.description == "a test vector"
    notes = _packet(to_bytes(record), 0x06)
    assert "<References>" in notes and 'title="A paper"' in notes
    assert "<Description>a test vector</Description>" in notes
    assert "<UUID>abc</UUID>" in notes
    assert "AccessionNumber" not in notes  # nothing was renamed, so nothing invented
    assert "2024.1.2</LastModified>" not in notes  # but the edit date moves on


def test_redescribing_the_record_updates_the_notes():
    record = read_snapgene(_file(""))
    record.id = "pNew"
    record.description = "a modified vector"
    notes = _packet(to_bytes(record), 0x06)
    assert "<Description>a modified vector</Description>" in notes
    assert "<References>" in notes
    # AccessionNumber means a GenBank accession to SnapGene; never invent one
    assert "AccessionNumber" not in notes


def test_record_without_a_name_gets_the_file_name(tmp_path):
    path = tmp_path / "pThing.dna"
    path.write_bytes(_file(""))
    assert read_snapgene(path).id == "pThing"


def test_sequence_case_is_kept():
    out = to_bytes(read_snapgene(_file("")))
    dna = next(p.data for p in iter_packets(out) if p.type == 0x00)
    assert dna[1:].decode() == SEQ


def test_fresh_record_writes_no_placeholder_identity():
    record = SeqRecord(Seq("ATGC" * 10))
    notes = _packet(to_bytes(record), 0x06)
    assert "unknown" not in notes
    assert "AccessionNumber" not in notes


def test_new_primer_is_added_alongside_existing_features():
    feature = (
        '<Feature name="M13 rev" type="primer_bind" directionality="1">'
        '<Segment range="50-66" color="#a020f0"/></Feature>'
    )
    record = read_snapgene(_file(feature))
    record.features.append(
        SeqFeature(
            SimpleLocation(60, 78, strand=1), type="primer_bind", qualifiers={"label": ["new_fwd"]}
        )
    )
    back = read_snapgene(to_bytes(record))
    labels = sorted(f.qualifiers["label"][0] for f in back.features)
    assert labels == ["M13 rev", "new_fwd"]


# --- byte-for-byte fidelity ---------------------------------------------------

SNAPGENE_STYLE = (
    '<Feature recentID="0" name="box" type="misc_feature" allowSegmentOverlaps="0" '
    'consecutiveTranslationNumbering="1"><Segment range="50-60" color="#99ccff" '
    'type="standard"/><Q name="gene"><V text="&lt;html&gt;x&lt;/html&gt;"/></Q>'
    '<Q name="note"><V text="a &quot;quoted&quot; note"/></Q></Feature>'
    '<Feature recentID="3" name="arrow" directionality="2" type="promoter" '
    'swappedSegmentNumbering="1" allowSegmentOverlaps="0" '
    'consecutiveTranslationNumbering="1"><Segment range="11-20" color="#ffffff" '
    'type="standard"/><Q name="codon_start"><V int="1"/></Q>'
    '<Q name="direction"><V predef="LEFT"/></Q></Feature>'
)
SNAPGENE_PRIMER = (
    '<Primer recentID="0" name="fwd" sequence="TTTTTTTTTTTTTTTTTT" description="d">'
    '<BindingSite location="60-77" boundStrand="0" annealedBases="TTTTTTTTTTTTTTTTTT" '
    'meltingTemperature="40"><Component hybridizedRange="60-77" '
    'bases="TTTTTTTTTTTTTTTTTT"/></BindingSite><BindingSite simplified="1" '
    'location="60-77" boundStrand="0" annealedBases="TTTTTTTTTTTTTTTTTT" '
    'meltingTemperature="40"><Component hybridizedRange="60-77" '
    'bases="TTTTTTTTTTTTTTTTTT"/></BindingSite></Primer>'
)


def _snapgene_file() -> bytes:
    """A file laid out exactly as SnapGene writes one."""
    cookie = pack_packet(0x09, struct.pack(">8sHHH", b"SnapGene", 1, 15, 19))
    dna = pack_packet(0x00, b"\x03" + SEQ.encode())
    extra = pack_packet(0x11, b'<AlignableSequences trimStringency="Medium"/>')
    features = pack_packet(
        0x0A,
        f'<?xml version="1.0"?><Features nextValidID="4">{SNAPGENE_STYLE}</Features>\n'.encode(),
    )
    params = (
        '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1" '
        'minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" '
        'minimumFivePrimeAnnealing="15"/>'
    )
    primers = pack_packet(
        0x05,
        (
            f'<?xml version="1.0"?><Primers nextValidID="1">{params}{SNAPGENE_PRIMER}</Primers>\n'
        ).encode(),
    )
    notes = pack_packet(0x06, NOTES.encode())
    tail = pack_packet(0x1C, b'<?xml version="1.0"?><EnzymeVisibilities vals=""/>\n')
    return cookie + dna + extra + features + primers + notes + tail


def test_unedited_file_is_written_back_byte_for_byte():
    data = _snapgene_file()
    assert to_bytes(read_snapgene(data)) == data


def test_feature_without_direction_stays_without_one():
    record = read_snapgene(_snapgene_file())
    box = [f for f in record.features if f.qualifiers["label"] == ["box"]][0]
    assert box.location.strand is None
    assert 'name="box" type=' in _packet(to_bytes(record), 0x0A)  # no directionality


def test_editing_keeps_everything_else_in_place():
    data = _snapgene_file()
    record = read_snapgene(data)
    record.features.append(
        SeqFeature(SimpleLocation(30, 40), type="misc_feature", qualifiers={"label": ["new"]})
    )
    out = to_bytes(record)
    assert [p.type for p in iter_packets(out)] == [p.type for p in iter_packets(data)]
    features = _packet(out, 0x0A)
    assert SNAPGENE_STYLE in features  # existing features untouched
    assert 'recentID="4" name="new" type="misc_feature"' in features
    assert 'nextValidID="5"' in features
    assert _packet(out, 0x05) == _packet(data, 0x05)  # primers untouched
    assert "2024.1.2</LastModified>" not in _packet(out, 0x06)  # edit date moves on


def test_new_primer_is_written_the_way_snapgene_writes_one():
    record = read_snapgene(_snapgene_file())
    record.features.append(
        SeqFeature(
            SimpleLocation(10, 28, strand=-1),
            type="primer_bind",
            qualifiers={"label": ["rev"], "note": ["a reverse primer"]},
        )
    )
    primers = _packet(to_bytes(record), 0x05)
    assert SNAPGENE_PRIMER in primers  # the existing one is untouched
    new = re.search(r'<Primer recentID="1" name="rev".*?</Primer>', primers).group(0)
    bases = str(Seq(SEQ[10:28]).reverse_complement()).upper()
    assert f'sequence="{bases}" description="a reverse primer"' in new
    assert new.count("<BindingSite") == 2 and 'simplified="1"' in new
    assert f'<Component hybridizedRange="10-27" bases="{bases}"/>' in new
    assert re.search(r'meltingTemperature="\d+"', new)
    assert 'nextValidID="2"' in primers


def test_fresh_record_gets_snapgenes_standard_packets():
    out = to_bytes(SeqRecord(Seq("ATGC" * 10)))
    kinds = [p.type for p in iter_packets(out)]
    assert kinds == [0x09, 0x00, 0x11, 0x08, 0x0A, 0x06, 0x0E, 0x1C]


# --- details that cost real map content when they were wrong --------------------


def test_explicit_label_qualifier_is_not_dropped():
    """SnapGene usually leaves the label in the name attribute, but not always.

    A file that carries an explicit ``<Q name="label">`` must get it back: dropping it
    rewrites the feature XML and loses the byte-for-byte round trip.
    """
    feature = (
        '<Feature recentID="0" name="insert" type="misc_feature" directionality="1">'
        '<Segment range="1-30" color="#a6acb3" type="standard"/>'
        '<Q name="label"><V text="insert"/></Q></Feature>'
    )
    written = to_bytes(read_snapgene(_file(feature)))
    assert '<Q name="label"><V text="insert"/></Q>' in _packet(written, 0x0A)


def test_label_is_not_invented_when_the_file_has_none():
    feature = (
        '<Feature recentID="0" name="insert" type="misc_feature" directionality="1">'
        '<Segment range="1-30" color="#a6acb3" type="standard"/></Feature>'
    )
    written = to_bytes(read_snapgene(_file(feature)))
    assert '<Q name="label"' not in _packet(written, 0x0A)


def test_same_name_and_span_but_different_qualifiers_are_both_kept():
    """Two annotations alike in name, span and type are still two annotations."""
    feature = (
        '<Feature recentID="0" name="insert" type="misc_feature" directionality="1">'
        '<Segment range="1-30" color="#a6acb3" type="standard"/></Feature>'
        '<Feature recentID="1" name="insert" type="misc_feature" directionality="1">'
        '<Segment range="1-30" color="#a6acb3" type="standard"/>'
        '<Q name="note"><V text="the other one"/></Q></Feature>'
    )
    record = read_snapgene(_file(feature))
    assert len(record.features) == 2
    assert _packet(to_bytes(record), 0x0A).count("<Feature ") == 2


def test_short_binding_site_from_a_file_is_kept():
    """A sub-15 bp site SnapGene wrote is an extra 5' match, not a primer of its own.

    The length floor exists for primers built in Python, which SnapGene erases on open.
    """
    primer = (
        '<Primer recentID="0" name="p1" sequence="CCCCCCCCCC">'
        '<BindingSite location="10-19" boundStrand="0" annealedBases="CCCCCCCCCC"'
        ' meltingTemperature="28"/></Primer>'
    )
    record = read_snapgene(_file("", primer))
    assert len(record.features) == 1
    assert 'location="10-19"' in _packet(to_bytes(record), 0x05)


def test_binding_site_across_the_origin_round_trips():
    """A site that wraps the origin used to raise ValueError on read."""
    primer = (
        '<Primer recentID="0" name="wrap" sequence="TTTTTCCCCC">'
        '<BindingSite location="90-4" boundStrand="0" annealedBases="TTTTTCCCCC"'
        ' meltingTemperature="28"/></Primer>'
    )
    record = read_snapgene(_file("", primer))
    site = record.features[0]
    assert len(site.location.parts) == 2, "a wrapping site is two segments"
    assert (int(site.location.parts[0].start), int(site.location.parts[1].end)) == (90, 5)
    assert 'location="90-4"' in _packet(to_bytes(record), 0x05)


def test_translated_feature_with_a_partial_codon_does_not_warn():
    """A clipped feature keeps SnapGene's translate setting, without a partial codon.

    Biopython warns when asked to translate a span that is not a whole number of codons,
    and is due to make that an error, so the writer trims before translating.
    """
    import warnings

    feature = (
        '<Feature recentID="0" name="clipped" type="CDS" directionality="1">'
        '<Segment range="30-51" color="#75c6a9" type="standard" translated="1"/>'
        "</Feature>"
    )  # 22 bases: deliberately not a multiple of three
    record = read_snapgene(_file(feature))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        written = to_bytes(record)
    xml = _packet(written, 0x0A)
    assert 'translated="1"' in xml, "SnapGene's own setting must survive"
    assert '<Q name="translation"' in xml


def test_out_of_range_feature_is_refused_not_written():
    """A location outside the sequence must raise, not produce an unreadable file.

    Writing it emits a SnapGene range like ``-2-4``, which the reader cannot parse, so the
    failure would otherwise surface much later as a corrupt map.
    """
    import pytest

    record = SeqRecord(Seq(SEQ), id="pTest", description="pTest")
    record.annotations["topology"] = "circular"
    record.features.append(
        SeqFeature(
            SimpleLocation(-3, 3, strand=1),  # a site straddling the origin, written wrongly
            type="protein_bind",
            qualifiers={"label": ["XhoI"]},
        )
    )
    with pytest.raises(ValueError, match="outside the"):
        to_bytes(record)


def test_origin_crossing_feature_in_a_compound_location_is_fine():
    """The right way to say the same thing, and it round-trips."""
    from Bio.SeqFeature import CompoundLocation

    record = SeqRecord(Seq(SEQ), id="pTest", description="pTest")
    record.annotations["topology"] = "circular"
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [SimpleLocation(len(SEQ) - 3, len(SEQ), strand=1), SimpleLocation(0, 3, strand=1)]
            ),
            type="protein_bind",
            qualifiers={"label": ["XhoI"]},
        )
    )
    back = read_snapgene(to_bytes(record))
    site = [f for f in back.features if f.qualifiers.get("label", [""])[0] == "XhoI"][0]
    assert len(site.location.parts) == 2
    assert int(site.location.parts[0].start) == len(SEQ) - 3
    assert int(site.location.parts[1].end) == 3
