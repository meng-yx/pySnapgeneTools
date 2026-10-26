"""Read a SnapGene ``.dna`` file into a Biopython :class:`SeqRecord`.

Biopython can already parse SnapGene, but it discards everything that is not part
of the generic sequence model: segment colours, the native feature type spelling,
the translation flag, and every packet it does not recognise. This reader keeps all
of it, plus enough bookkeeping (attribute order, IDs, packet order) that an unedited
map is written back byte for byte.
"""

from __future__ import annotations

import re
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
    Q_ATTR_ORDER,
    Q_ATTR_PREFIX,
    Q_COLOR,
    Q_NATIVE_TYPE,
    Q_PREDEF,
    Q_PRIMER_SEQUENCE,
    Q_SEGMENT_COLORS,
    Q_SEGMENT_NAMES,
    Q_SITE_XML,
    Q_TRANSLATED,
)
from .packets import COOKIE, DNA, FEATURES, NOTES, PRIMERS, Packet, iter_packets

__all__ = ["read_snapgene", "SnapGeneError", "STRAND_FROM_DIRECTIONALITY"]

#: SnapGene ``directionality`` -> Biopython strand. No attribute means no
#: direction, which SnapGene draws as a box rather than an arrow.
STRAND_FROM_DIRECTIONALITY = {"1": 1, "2": -1, "3": 0}


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


def _parse_range(rangespec: str, strand: int | None, length: int) -> SimpleLocation:
    """SnapGene ranges are 1-based inclusive and may wrap the origin."""
    start_text, _, end_text = rangespec.partition("-")
    start = int(start_text) - 1
    end = int(end_text)
    if start > end:  # wraps the origin
        before = SimpleLocation(start, length, strand=strand)
        after = SimpleLocation(0, end, strand=strand)
        # Biopython lists parts in the order they are read, 5' to 3'
        return CompoundLocation([after, before] if strand == -1 else [before, after])
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


def _next_valid_id(xml: Document, tag: str) -> int | None:
    root = xml.getElementsByTagName(tag)
    if root and root[0].hasAttribute("nextValidID"):
        try:
            return int(root[0].getAttribute("nextValidID"))
        except ValueError:
            return None
    return None


def _parse_features(packet: Packet, record: SeqRecord) -> None:
    xml = parseString(packet.data.decode("utf-8"))
    record.annotations["snapgene_features_next_id"] = _next_valid_id(xml, "Features")
    length = len(record.seq)
    for node in xml.getElementsByTagName("Feature"):
        strand = STRAND_FROM_DIRECTIONALITY.get(_text(node, "directionality") or "")

        location = None
        colors: list[str] = []
        segment_names: list[str] = []
        translated = False
        for segment in node.getElementsByTagName("Segment"):
            if _text(segment, "type", "standard") == "gap":
                continue
            rangespec = _text(segment, "range")
            if not rangespec:
                continue
            colors.append(_text(segment, "color", DEFAULT_COLOR) or DEFAULT_COLOR)
            segment_names.append(_text(segment, "name") or "")
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

        qualifiers: dict[str, list[Any]] = {}
        predef: list[str] = []
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
                        if attribute == "predef" and name not in predef:
                            predef.append(name)
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
        if any(segment_names):
            qualifiers[Q_SEGMENT_NAMES] = segment_names
        if translated:
            qualifiers[Q_TRANSLATED] = ["1"]
        if predef:
            qualifiers[Q_PREDEF] = predef
        # keep every other Feature attribute, and their order, for the writer
        qualifiers[Q_ATTR_ORDER] = list(node.attributes.keys())
        for key in node.attributes.keys():
            if key in ("name", "type", "directionality"):
                continue
            qualifiers[f"{Q_ATTR_PREFIX}{key}"] = [node.attributes[key].value]

        record.features.append(
            SeqFeature(location, type=native_type.rstrip("_"), qualifiers=qualifiers)
        )


_PRIMER_RE = re.compile(r"<Primer\b.*?</Primer>|<Primer\b[^>]*/>", re.S)
_SITE_RE = re.compile(r"<BindingSite\b.*?</BindingSite>|<BindingSite\b[^>]*/>", re.S)
_SITE_KEY_RE = re.compile(r'\blocation="(\d+)-(\d+)"[^>]*?\bboundStrand="(\d)"')


def _parse_primers(packet: Packet, record: SeqRecord) -> None:
    text = packet.data.decode("utf-8")
    xml = parseString(text)
    record.annotations["snapgene_primers_next_id"] = _next_valid_id(xml, "Primers")
    params = re.search(r"<HybridizationParams\b[^>]*/>", text)
    if params:
        record.annotations["snapgene_hybridization_xml"] = params.group(0)

    length = len(record.seq)
    raw_primers = _PRIMER_RE.findall(text)
    nodes = xml.getElementsByTagName("Primer")
    for index, node in enumerate(nodes):
        qualifiers: dict[str, list[Any]] = {}
        name = _text(node, "name")
        if name:
            qualifiers["label"] = [name]
        sequence = _text(node, "sequence")
        if sequence:
            qualifiers[Q_PRIMER_SEQUENCE] = [sequence]
        description = _text(node, "description")
        if description:
            qualifiers["note"] = [description]
        qualifiers[Q_ATTR_ORDER] = list(node.attributes.keys())
        for key in node.attributes.keys():
            if key in ("name", "sequence", "description"):
                continue
            qualifiers[f"{Q_ATTR_PREFIX}{key}"] = [node.attributes[key].value]

        # the binding sites' own XML, grouped by place, so unedited ones go back as-is
        raw_sites: dict[tuple[int, int, str], list[str]] = {}
        if index < len(raw_primers):
            for site in _SITE_RE.findall(raw_primers[index]):
                match = _SITE_KEY_RE.search(site)
                if match:
                    key = (int(match.group(1)), int(match.group(2)), match.group(3))
                    raw_sites.setdefault(key, []).append(site)

        placed: set[tuple[int, int, int]] = set()
        for site in node.getElementsByTagName("BindingSite"):
            rangespec = _text(site, "location")
            if not rangespec:
                continue
            bound = _text(site, "boundStrand", "0") or "0"
            strand = -1 if bound == "1" else 1
            # primer coordinates sit one lower than a feature's at both ends
            start_text, _, end_text = rangespec.partition("-")
            location = SimpleLocation(int(start_text), int(end_text) + 1, strand=strand)
            key = (int(location.start), int(location.end), strand)
            if key in placed:
                continue  # SnapGene repeats each site as a "simplified" duplicate
            placed.add(key)
            if int(location.end) > length:
                continue
            site_qualifiers = dict(qualifiers)
            raw = raw_sites.get((int(start_text), int(end_text), bound))
            if raw:
                site_qualifiers[Q_SITE_XML] = raw
            record.features.append(
                SeqFeature(location, type="primer_bind", qualifiers=site_qualifiers)
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
    # the whole packet too, so references and other fields survive a rewrite
    record.annotations["snapgene_notes_xml"] = packet.data.decode("utf-8")
    if notes.get("Type") == "Synthetic":
        record.annotations["data_file_division"] = "SYN"
    accession = notes.get("AccessionNumber")
    description = notes.get("Description") or notes.get("Comments")
    if accession:
        record.id = accession
        record.name = accession[:20]
    if description:
        record.description = description


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
        :func:`pysnapgene_writer.write_snapgene` can put them back.

    Returns
    -------
    SeqRecord
        Features carry their SnapGene colour in the ``snapgene_color`` qualifier,
        the native type spelling in ``snapgene_native_type``, and ``snapgene_translated``
        when SnapGene shows a translation for them. A feature SnapGene draws without
        an arrow has ``strand=None``.
    """
    if isinstance(source, (str, Path)):
        data = Path(source).read_bytes()
    elif isinstance(source, (bytes, bytearray)):
        data = bytes(source)
    else:
        data = source.read()

    record = SeqRecord(Seq(""), id="", name="", description="")
    unknown: list[Packet] = []
    order: list[int] = []
    originals: dict[int, bytes] = {}
    saw_cookie = False
    saw_dna = False

    for packet in iter_packets(data):
        order.append(packet.type)
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
        if packet.type in (DNA, FEATURES, PRIMERS, NOTES):
            originals[packet.type] = packet.data

    if not saw_cookie:
        raise SnapGeneError("not a SnapGene file: no cookie packet")
    if not saw_dna:
        raise SnapGeneError("file contains no DNA packet")

    if not record.id and isinstance(source, (str, Path)):
        record.id = record.name = Path(source).stem
    # lets the writer tell whether the caller renamed or redescribed the record
    record.annotations["snapgene_identity"] = (record.id, record.description)
    # lets the writer keep SnapGene's packet order and leave the notes untouched
    # when nothing else changed
    record.annotations["snapgene_packet_order"] = order
    record.annotations["snapgene_original_packets"] = originals

    if preserve_unknown_packets:
        record.annotations["snapgene_unknown_packets"] = unknown
    return record
