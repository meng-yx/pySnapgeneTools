"""Read a SnapGene ``.dna`` file into a Biopython :class:`SeqRecord`.

Biopython can already parse SnapGene, but it discards everything that is not part
of the generic sequence model: segment colours, the native feature type spelling,
the translation flag, and every packet it does not recognise. This reader keeps all
of it, so a file can be read, edited and written back without losing the map.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import Any, BinaryIO
from xml.dom.minidom import Document, parseString

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .constants import (
    DEFAULT_COLOR,
    FLAG_CIRCULAR,
    FLAG_DOUBLE_STRANDED,
    Q_ATTR_PREFIX,
    Q_COLOR,
    Q_NATIVE_TYPE,
    Q_SEGMENT_COLORS,
    Q_TRANSLATED,
)
from .packets import COOKIE, DNA, FEATURES, NOTES, PRIMERS, Packet, iter_packets

__all__ = ["read_snapgene", "SnapGeneError"]


class SnapGeneError(ValueError):
    """Raised when a file is not valid SnapGene."""


def _text(node: Any, attribute: str, default: str | None = None) -> str | None:
    return node.getAttribute(attribute) if node.hasAttribute(attribute) else default


def _child_text(doc: Document, tag: str) -> str | None:
    nodes = doc.getElementsByTagName(tag)
    if nodes and nodes[0].childNodes:
        first = nodes[0].firstChild
        if first is not None and first.nodeType == first.TEXT_NODE:
            return first.data
    return None


def _parse_range(rangespec: str, strand: int, length: int) -> SimpleLocation:
    """SnapGene ranges are 1-based inclusive and may wrap the origin."""
    start_text, _, end_text = rangespec.partition("-")
    start = int(start_text) - 1
    end = int(end_text)
    if start > end:  # wraps the origin
        return CompoundLocation(
            [
                SimpleLocation(start, length, strand=strand),
                SimpleLocation(0, end, strand=strand),
            ]
        )
    return SimpleLocation(start, end, strand=strand)


def _parse_dna(packet: Packet, record: SeqRecord) -> None:
    if not packet.data:
        raise SnapGeneError("empty DNA packet")
    flags = packet.data[0]
    record.seq = Seq(packet.data[1:].decode("ascii"))
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular" if flags & FLAG_CIRCULAR else "linear"
    record.annotations["snapgene_flags"] = flags
    record.annotations["double_stranded"] = bool(flags & FLAG_DOUBLE_STRANDED)


def _parse_features(packet: Packet, record: SeqRecord) -> None:
    xml = parseString(packet.data.decode("utf-8"))
    length = len(record.seq)
    for node in xml.getElementsByTagName("Feature"):
        strand = -1 if _text(node, "directionality", "1") == "2" else 1

        location = None
        colors: list[str] = []
        translated = False
        for segment in node.getElementsByTagName("Segment"):
            if _text(segment, "type", "standard") == "gap":
                continue
            rangespec = _text(segment, "range")
            if not rangespec:
                continue
            colors.append(_text(segment, "color", DEFAULT_COLOR) or DEFAULT_COLOR)
            if _text(segment, "translated") == "1":
                translated = True
            piece = _parse_range(rangespec, strand, length)
            if location is None:
                location = piece
            elif strand == -1:
                location = piece + location
            else:
                location = location + piece
        if location is None:
            continue  # a feature with no usable segment is not placeable

        qualifiers: dict[str, list[str]] = {}
        for qualifier in node.getElementsByTagName("Q"):
            name = _text(qualifier, "name")
            if not name:
                continue
            values: list[Any] = []
            for value in qualifier.getElementsByTagName("V"):
                for attribute in ("text", "predef", "int"):
                    if value.hasAttribute(attribute):
                        raw = value.attributes[attribute].value
                        values.append(int(raw) if attribute == "int" else raw)
                        break
            qualifiers[name] = values

        name = _text(node, "name")
        if name and "label" not in qualifiers:
            qualifiers["label"] = [name]

        native_type = _text(node, "type", "misc_feature") or "misc_feature"
        qualifiers[Q_NATIVE_TYPE] = [native_type]
        qualifiers[Q_COLOR] = [colors[0]]
        if len(set(colors)) > 1:
            qualifiers[Q_SEGMENT_COLORS] = colors
        if translated:
            qualifiers[Q_TRANSLATED] = ["1"]
        # keep every other Feature attribute so the writer can put it back
        for key in node.attributes.keys():
            if key in ("name", "type", "directionality"):
                continue
            qualifiers[f"{Q_ATTR_PREFIX}{key}"] = [node.attributes[key].value]

        record.features.append(
            SeqFeature(location, type=native_type.rstrip("_"), qualifiers=qualifiers)
        )


def _parse_primers(packet: Packet, record: SeqRecord) -> None:
    xml = parseString(packet.data.decode("utf-8"))
    length = len(record.seq)
    for node in xml.getElementsByTagName("Primer"):
        qualifiers: dict[str, list[str]] = {}
        name = _text(node, "name")
        if name:
            qualifiers["label"] = [name]
        sequence = _text(node, "sequence")
        if sequence:
            qualifiers["snapgene_primer_sequence"] = [sequence]

        placed: set[tuple[int, int, int]] = set()
        for site in node.getElementsByTagName("BindingSite"):
            rangespec = _text(site, "location")
            if not rangespec:
                continue
            strand = -1 if _text(site, "boundStrand", "0") == "1" else 1
            # primer coordinates sit one lower than a feature's at both ends
            start_text, _, end_text = rangespec.partition("-")
            location = SimpleLocation(int(start_text), int(end_text) + 1, strand=strand)
            key = (int(location.start), int(location.end), strand)
            if key in placed:
                continue  # SnapGene repeats each site as a "simplified" duplicate
            placed.add(key)
            if int(location.end) > length:
                continue
            record.features.append(
                SeqFeature(location, type="primer_bind", qualifiers=dict(qualifiers))
            )


def _parse_notes(packet: Packet, record: SeqRecord) -> None:
    xml = parseString(packet.data.decode("utf-8"))
    notes: dict[str, str] = {}
    for tag in (
        "UUID",
        "Type",
        "Created",
        "LastModified",
        "AccessionNumber",
        "Comments",
        "Description",
    ):
        value = _child_text(xml, tag)
        if value is not None:
            notes[tag] = value
    record.annotations["snapgene_notes"] = notes
    if notes.get("Type") == "Synthetic":
        record.annotations["data_file_division"] = "SYN"
    accession = notes.get("AccessionNumber")
    comments = notes.get("Comments")
    if accession:
        record.id = accession
    if comments:
        record.description = comments
        if not accession:
            record.id = comments.split(" ", 1)[0]
        record.name = (accession or comments.split(" ", 1)[0])[:20]


def read_snapgene(
    source: str | Path | bytes | BinaryIO,
    *,
    preserve_unknown_packets: bool = True,
) -> SeqRecord:
    """Read a SnapGene ``.dna`` file.

    Parameters
    ----------
    source:
        Path, raw bytes, or an open binary file object.
    preserve_unknown_packets:
        Keep packets this library does not model (enzyme sets, alignments, edit
        history) in ``record.annotations["snapgene_unknown_packets"]`` so that
        :func:`pysnapgene.write_snapgene` can put them back.

    Returns
    -------
    SeqRecord
        Features carry their SnapGene colour in the ``snapgene_color`` qualifier,
        the native type spelling in ``snapgene_native_type``, and ``snapgene_translated``
        when SnapGene shows a translation for them.
    """
    if isinstance(source, (str, Path)):
        data = Path(source).read_bytes()
    elif isinstance(source, (bytes, bytearray)):
        data = bytes(source)
    else:
        data = source.read()

    record = SeqRecord(Seq(""), id="", name="", description="")
    unknown: list[Packet] = []
    saw_cookie = False
    saw_dna = False

    for packet in iter_packets(data):
        if packet.type == COOKIE:
            if len(packet.data) < 14:
                raise SnapGeneError("short cookie packet")
            magic, seq_type, export_version, import_version = struct.unpack(
                ">8sHHH", packet.data[:14]
            )
            if magic != b"SnapGene":
                raise SnapGeneError("not a SnapGene file: bad magic cookie")
            record.annotations["snapgene_cookie"] = (
                seq_type,
                export_version,
                import_version,
            )
            saw_cookie = True
        elif packet.type == DNA:
            if saw_dna:
                raise SnapGeneError("more than one DNA packet")
            _parse_dna(packet, record)
            saw_dna = True
        elif packet.type == FEATURES:
            _parse_features(packet, record)
        elif packet.type == PRIMERS:
            _parse_primers(packet, record)
        elif packet.type == NOTES:
            _parse_notes(packet, record)
        elif preserve_unknown_packets:
            unknown.append(packet)

    if not saw_cookie:
        raise SnapGeneError("not a SnapGene file: no cookie packet")
    if not saw_dna:
        raise SnapGeneError("file contains no DNA packet")

    if preserve_unknown_packets:
        record.annotations["snapgene_unknown_packets"] = unknown
    return record
