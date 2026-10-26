"""Write a Biopython :class:`SeqRecord` out as a SnapGene ``.dna`` file.

Biopython has no SnapGene writer -- ``snapgene`` does not appear in
``Bio.SeqIO._FormatToWriter`` -- so this module supplies one.

The XML is laid out the way SnapGene lays it out (attribute order, qualifiers
sorted by name, each primer binding site with its melting temperature and a
"simplified" twin, packets in SnapGene's order), so a map that was read and not
edited is written back byte for byte, and an edited one looks like SnapGene saved
it.

Three details decide whether a map built from scratch looks right:

* the **double-stranded bit** in the DNA packet flag byte. Leave it clear and the
  molecule is shown as single stranded;
* a **colour on every segment**. SnapGene does not invent one, so a file written
  without the attribute is uniformly grey;
* ``translated="1"`` plus the usual CDS qualifiers, which is what ticks
  *Translate this feature in Sequence view*.
"""

from __future__ import annotations

import re
import struct
import uuid
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, BinaryIO

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .constants import (
    DEFAULT_COLOR,
    DEFAULT_COOKIE,
    DEFAULT_DNA_FLAGS,
    FLAG_CIRCULAR,
    MIN_PRIMER_LENGTH,
    Q_ATTR_ORDER,
    Q_ATTR_PREFIX,
    Q_COLOR,
    Q_LABEL_IS_QUALIFIER,
    Q_NATIVE_TYPE,
    Q_PREDEF,
    Q_PRIMER_SEQUENCE,
    Q_SEGMENT_COLORS,
    Q_SEGMENT_NAMES,
    Q_SITE_RANGE,
    Q_SITE_XML,
    Q_TRANSLATED,
    RESERVED_QUALIFIERS,
)
from .packets import COOKIE, DNA, FEATURES, NOTES, PRIMERS, Packet, pack_packet

__all__ = ["write_snapgene", "to_bytes"]

#: Feature types SnapGene will offer a translation for.
_TRANSLATABLE = {"CDS", "gene", "mat_peptide", "sig_peptide", "transit_peptide"}

#: Biopython strand -> SnapGene ``directionality``; ``None`` writes no attribute.
_DIRECTIONALITY = {1: "1", -1: "2", 0: "3"}

#: The order SnapGene writes attributes in, for attributes a feature did not
#: arrive with.
_FEATURE_ATTR_ORDER = [
    "recentID",
    "name",
    "directionality",
    "translationMW",
    "type",
    "swappedSegmentNumbering",
    "allowSegmentOverlaps",
    "cleavageArrows",
    "readingFrame",
    "consecutiveTranslationNumbering",
    "hitsStopCodon",
    "detectionMode",
]
_PRIMER_ATTR_ORDER = ["recentID", "name", "sequence", "description"]

#: SnapGene's packet order. Packets missing from a file are slotted in by it.
_PACKET_ORDER = [0x09, 0x00, 0x03, 0x0B, 0x07, 0x11, 0x08, 0x0A, 0x05, 0x06, 0x0D, 0x0E, 0x1C]

#: Small packets SnapGene writes into every file, with their default content.
#: A file built from scratch gets them too, so it looks like SnapGene wrote it.
_DEFAULT_PACKETS = {
    0x11: b'<AlignableSequences trimStringency="Medium"/>',
    0x08: (
        b"<AdditionalSequenceProperties><UpstreamStickiness>0</UpstreamStickiness>"
        b"<DownstreamStickiness>0</DownstreamStickiness>"
        b"<UpstreamModification>FivePrimePhosphorylated</UpstreamModification>"
        b"<DownstreamModification>FivePrimePhosphorylated</DownstreamModification>"
        b"</AdditionalSequenceProperties>"
    ),
    0x0E: b'<?xml version="1.0"?><CustomEnzymeSets/>\n',
    0x1C: b'<?xml version="1.0"?><EnzymeVisibilities vals=""/>\n',
}

_HYBRIDIZATION_PARAMS = (
    '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1" '
    'minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" '
    'minimumFivePrimeAnnealing="15"/>'
)


def _attr(value: object) -> str:
    """Escape an attribute value the way SnapGene does, in double quotes."""
    text = str(value)
    for raw, escaped in (
        ("&", "&amp;"),
        ("<", "&lt;"),
        (">", "&gt;"),
        ('"', "&quot;"),
        ("\n", "&#10;"),
        ("\r", "&#13;"),
        ("\t", "&#9;"),
    ):
        text = text.replace(raw, escaped)
    return f'"{text}"'


def _text_node(value: object) -> str:
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _render_attributes(attributes: dict[str, str], stored: list[str] | None, canon: list[str]):
    order = [k for k in (stored or []) if k in attributes]
    order += [k for k in canon if k in attributes and k not in order]
    order += [k for k in attributes if k not in order]
    return " ".join(f"{k}={_attr(attributes[k])}" for k in order)


def _label(feature: SeqFeature) -> str:
    values = feature.qualifiers.get("label")
    if values:
        return str(values[0])
    return feature.type or "feature"


def _first(feature: SeqFeature, key: str) -> str | None:
    values = feature.qualifiers.get(key)
    return str(values[0]) if values else None


def _stored_attributes(feature: SeqFeature) -> dict[str, str]:
    return {
        key[len(Q_ATTR_PREFIX) :]: str(values[0])
        for key, values in feature.qualifiers.items()
        if key.startswith(Q_ATTR_PREFIX) and values
    }


def _segment_colors(feature: SeqFeature, n_parts: int, fallback: str) -> list[str]:
    explicit = feature.qualifiers.get(Q_SEGMENT_COLORS)
    if explicit and len(explicit) == n_parts:
        return [str(c) for c in explicit]
    single = _first(feature, Q_COLOR) or fallback
    return [single] * n_parts


def _segment_names(feature: SeqFeature, n_parts: int) -> list[str]:
    names = feature.qualifiers.get(Q_SEGMENT_NAMES)
    if names and len(names) == n_parts:
        return [str(n) for n in names]
    return [""] * n_parts


def _file_segments(feature: SeqFeature, length: int) -> list[tuple[int, int]]:
    """Segment ranges as SnapGene stores them: 1-based inclusive, ascending.

    Biopython lists a minus-strand feature's parts 5' to 3', i.e. descending, and
    splits a feature that crosses the origin in two; SnapGene wants ascending
    order and a single range such as ``2470-372`` for the origin-crossing part.
    """
    parts = (
        feature.location.parts
        if isinstance(feature.location, CompoundLocation)
        else [feature.location]
    )
    if feature.location.strand == -1:
        parts = list(reversed(parts))
    ranges = [(int(p.start), int(p.end)) for p in parts]
    merged: list[tuple[int, int]] = []
    for start, end in ranges:
        if merged and merged[-1][1] == length and start == 0:
            merged[-1] = (merged[-1][0], end)  # rejoin across the origin
        else:
            merged.append((start, end))
    return [(start + 1, end) for start, end in merged]


def _is_clean_orf(nucleotides: str) -> tuple[bool, str]:
    """A span is only annotated as a translated CDS if it really is one.

    Feature boundaries carried over from another map are often not in frame, and
    marking those as CDS makes SnapGene display a nonsense protein.
    """
    if len(nucleotides) < 3 or len(nucleotides) % 3:
        return False, ""
    protein = str(Seq(nucleotides).translate())
    if "*" in protein[:-1]:
        return False, ""
    return True, protein


def _feature_nt(sequence: str, feature: SeqFeature) -> str:
    # extract() reverse-complements each part, which is right for multi-part
    # minus-strand features such as a CDS split around a signal peptide
    return str(feature.location.extract(Seq(sequence)))


def _value_xml(name: str, value: Any, predef: set[str]) -> str:
    if isinstance(value, int) and not isinstance(value, bool):
        return f'<V int="{value}"/>'
    kind = "predef" if name in predef else "text"
    return f"<V {kind}={_attr(value)}/>"


def _features_xml(
    sequence: str,
    features: Sequence[SeqFeature],
    default_color: str,
    translate: bool,
    next_id: int | None,
) -> str:
    used = [
        int(v)
        for f in features
        if (v := _first(f, f"{Q_ATTR_PREFIX}recentID")) is not None and v.isdigit()
    ]
    fresh_id = max([next_id or 0, *(i + 1 for i in used)])

    chunks: list[str] = []
    for feature in features:
        from_snapgene = Q_NATIVE_TYPE in feature.qualifiers
        label = _label(feature)
        native_type = _first(feature, Q_NATIVE_TYPE) or feature.type or "misc_feature"
        ranges = _file_segments(feature, len(sequence))
        colors = _segment_colors(feature, len(ranges), default_color)
        names = _segment_names(feature, len(ranges))
        stored = _stored_attributes(feature)

        # A feature SnapGene wrote keeps its own translate setting. A new one is
        # translated when it is a translatable type and a clean reading frame.
        as_cds = False
        protein = ""
        if from_snapgene:
            if _first(feature, Q_TRANSLATED) == "1":
                nucleotides = _feature_nt(sequence, feature)
                ok, protein = _is_clean_orf(nucleotides)
                as_cds = True
                if not ok:
                    # SnapGene's own setting is honoured even where the span is not a clean
                    # reading frame -- clipped features and alternative starts both occur --
                    # but the partial codon is trimmed, because Biopython warns on it now and
                    # is due to make it an error.
                    whole = nucleotides[: len(nucleotides) // 3 * 3]
                    protein = str(Seq(whole).translate()) if whole else ""
        elif (translate and native_type.rstrip("_") in _TRANSLATABLE) or _first(
            feature, Q_TRANSLATED
        ) == "1":
            as_cds, protein = _is_clean_orf(_feature_nt(sequence, feature))
            if not as_cds and native_type.rstrip("_") in _TRANSLATABLE:
                native_type = "misc_feature"

        segments = "".join(
            '<Segment{name} range="{start}-{end}" color={color} type="standard"{tr}/>'.format(
                name=f" name={_attr(seg_name)}" if seg_name else "",
                start=start,
                end=end,
                color=_attr(color),
                tr=' translated="1"' if as_cds else "",
            )
            for (start, end), color, seg_name in zip(ranges, colors, names)
        )

        qualifiers: dict[str, list[Any]] = {}
        for key, values in feature.qualifiers.items():
            if key in RESERVED_QUALIFIERS or key.startswith(Q_ATTR_PREFIX):
                continue
            if (
                key == "label"
                and values
                and str(values[0]) == label
                and Q_LABEL_IS_QUALIFIER not in feature.qualifiers
            ):
                continue  # the reader synthesised it from the name attribute
            qualifiers[key] = list(values)
        if as_cds:
            stored_translation = _first(feature, "translation")
            # keep SnapGene's own text when it says the same thing: it marks
            # cleavage sites with commas and shows Met for alternative starts
            same = stored_translation is not None and (
                stored_translation.replace(",", "")[1:] == protein[1:]
            )
            qualifiers["codon_start"] = [1]
            qualifiers["transl_table"] = [1]
            qualifiers["translation"] = [stored_translation if same else protein]
        predef = set(feature.qualifiers.get(Q_PREDEF, []))
        qualifier_xml = "".join(
            f"<Q name={_attr(key)}>"
            + "".join(_value_xml(key, v, predef) for v in qualifiers[key])
            + "</Q>"
            for key in sorted(qualifiers)
        )

        attributes: dict[str, str] = {}
        recent = stored.get("recentID")
        if recent is None:
            recent = str(fresh_id)
            fresh_id += 1
        attributes["recentID"] = recent
        attributes["name"] = label
        direction = _DIRECTIONALITY.get(feature.location.strand)
        if direction:
            attributes["directionality"] = direction
        attributes["type"] = native_type
        for key, value in stored.items():
            if key not in ("recentID", "name", "type", "directionality", "translated"):
                attributes[key] = value
        attributes.setdefault("allowSegmentOverlaps", "0")
        attributes.setdefault("consecutiveTranslationNumbering", "1")
        if as_cds and "translationMW" not in attributes:
            try:
                from Bio.SeqUtils import molecular_weight

                mass = molecular_weight(protein.replace("*", ""), seq_type="protein")
                attributes["translationMW"] = f"{mass:.2f}"
            except Exception:  # pragma: no cover - MW is cosmetic
                pass
        elif not as_cds:
            for key in ("translationMW", "readingFrame", "hitsStopCodon"):
                attributes.pop(key, None)

        rendered = _render_attributes(
            attributes, feature.qualifiers.get(Q_ATTR_ORDER), _FEATURE_ATTR_ORDER
        )
        chunks.append(f"<Feature {rendered}>{segments}{qualifier_xml}</Feature>")
    return f'<?xml version="1.0"?><Features nextValidID="{fresh_id}">{"".join(chunks)}</Features>\n'


def _melting_temperature(bases: str) -> int:
    """Close to SnapGene's figure (within a few degrees); SnapGene recalculates."""
    from Bio.SeqUtils import MeltingTemp

    tm = MeltingTemp.Tm_NN(
        bases, nn_table=MeltingTemp.DNA_NN3, Na=50, dnac1=500, dnac2=0, saltcorr=5
    )
    return int(round(tm))


def _site_xml(feature: SeqFeature, sequence: str, start: int, end: int, bound: int) -> str:
    # A site read from a file keeps its own range text, which is the only way an
    # origin-crossing site survives: a Biopython location reports the span's extremes.
    rangespec = _first(feature, Q_SITE_RANGE) or f"{start}-{end - 1}"
    location = f'location="{rangespec}" boundStrand="{bound}"'
    stored = feature.qualifiers.get(Q_SITE_XML)
    if stored and all(location in site for site in stored):
        return "".join(stored)
    annealed = sequence[start:end]
    if bound == 1:
        annealed = str(Seq(annealed).reverse_complement())
    tm = _melting_temperature(annealed)
    body = (
        f'annealedBases={_attr(annealed)} meltingTemperature="{tm}">'
        f'<Component hybridizedRange="{start}-{end - 1}" bases={_attr(annealed)}/>'
        "</BindingSite>"
    )
    return f'<BindingSite {location} {body}<BindingSite simplified="1" {location} {body}'


def _primers_xml(sequence: str, primers: Sequence[SeqFeature], record: SeqRecord) -> str:
    grouped: dict[tuple[str, str], list[SeqFeature]] = {}
    for feature in primers:
        start, end = int(feature.location.start), int(feature.location.end)
        bases = _first(feature, Q_PRIMER_SEQUENCE)
        if not bases:
            bases = sequence[start:end]
            if feature.location.strand == -1:
                bases = str(Seq(bases).reverse_complement())
        grouped.setdefault((_label(feature), bases), []).append(feature)

    used = [
        int(v)
        for f in primers
        if (v := _first(f, f"{Q_ATTR_PREFIX}recentID")) is not None and v.isdigit()
    ]
    next_id = record.annotations.get("snapgene_primers_next_id") or 0
    fresh_id = max([next_id, *(i + 1 for i in used)])

    chunks: list[str] = []
    for (label, bases), sites in grouped.items():
        head = sites[0]
        stored = _stored_attributes(head)
        recent = stored.get("recentID")
        if recent is None:
            recent = str(fresh_id)
            fresh_id += 1
        attributes = {"recentID": recent, "name": label, "sequence": bases}
        description = _first(head, "note")
        if description:
            attributes["description"] = description
        for key, value in stored.items():
            if key not in attributes:
                attributes[key] = value
        rendered = _render_attributes(
            attributes, head.qualifiers.get(Q_ATTR_ORDER), _PRIMER_ATTR_ORDER
        )
        binding = "".join(
            _site_xml(
                site,
                sequence,
                int(site.location.start),
                int(site.location.end),
                1 if site.location.strand == -1 else 0,
            )
            for site in sites
        )
        chunks.append(f"<Primer {rendered}>{binding}</Primer>")

    params = record.annotations.get("snapgene_hybridization_xml") or _HYBRIDIZATION_PARAMS
    return (
        '<?xml version="1.0"?>'
        f'<Primers nextValidID="{fresh_id}">{params}{"".join(chunks)}</Primers>\n'
    )


#: Biopython's defaults for a record nobody named.
_PLACEHOLDERS = {"", "<unknown id>", "<unknown name>", "<unknown description>"}


def _meaningful(value: object) -> str:
    text = str(value or "")
    return "" if text in _PLACEHOLDERS else text


def _replace_element(xml: str, tag: str, text: str, attributes: str = "") -> str:
    """Replace ``<tag ...>...</tag>`` in place, or add it just before ``</Notes>``."""
    element = f"<{tag}{attributes}>{_text_node(text)}</{tag}>"
    pattern = re.compile(rf"<{tag}\b[^>]*>.*?</{tag}>|<{tag}\b[^>]*/>", re.S)
    if pattern.search(xml):
        return pattern.sub(lambda _: element, xml, count=1)
    return xml.replace("</Notes>", f"{element}\n</Notes>", 1)


def _notes_xml(record: SeqRecord, edited: bool) -> str:
    now = datetime.now(timezone.utc)
    today = f"{now.year}.{now.month}.{now.day}"
    clock = f' UTC="{now.hour}:{now.minute}:{now.second}"'
    accession = _meaningful(record.id)
    description = _meaningful(record.description)

    raw = record.annotations.get("snapgene_notes_xml")
    if raw:
        # Start from the notes that were read, so references, sequence class and
        # everything else this library does not model are kept as they were.
        identity_changed = (record.id, record.description) != record.annotations.get(
            "snapgene_identity"
        )
        if not (edited or identity_changed):
            return raw
        notes = _replace_element(raw, "LastModified", today, clock)
        if identity_changed:
            if description:
                notes = _replace_element(notes, "Description", description)
            if accession and "<AccessionNumber" in notes:
                notes = _replace_element(notes, "AccessionNumber", accession)
        return notes

    stored = record.annotations.get("snapgene_notes", {}) or {}
    fields: list[tuple[str, str, str]] = [
        ("UUID", "", stored.get("UUID") or str(uuid.uuid4())),
        ("Type", "", stored.get("Type", "Synthetic")),
        ("ConfirmedExperimentally", "", "0"),
        ("Description", "", description),
        ("AccessionNumber", "", accession),
        ("Created", clock, stored.get("Created", today)),
        ("LastModified", clock, today),
        ("SequenceClass", "", "UNA"),
        ("TransformedInto", "", "unspecified"),
    ]
    body = "".join(f"<{k}{a}>{_text_node(v)}</{k}>\n" for k, a, v in fields if v)
    return f"<Notes>\n{body}</Notes>\n"


def _ordered(packets: dict[int, bytes], unknown: Sequence[Packet], order: list[int] | None):
    """Lay packets out in the file's own order, else SnapGene's."""
    queue = list(unknown)
    out: list[tuple[int, bytes]] = []
    if order:
        placed: set[int] = set()
        for kind in order:
            if kind in packets and kind not in placed:
                out.append((kind, packets[kind]))
                placed.add(kind)
            else:
                for i, packet in enumerate(queue):
                    if packet.type == kind:
                        out.append((packet.type, packet.data))
                        del queue[i]
                        break
        missing = [(k, v) for k, v in packets.items() if k not in placed]
    else:
        missing = list(packets.items())
    missing += [(p.type, p.data) for p in queue]

    def rank(kind: int) -> int:
        return _PACKET_ORDER.index(kind) if kind in _PACKET_ORDER else len(_PACKET_ORDER)

    for kind, data in missing:
        at = next((i for i, (k, _) in enumerate(out) if rank(k) > rank(kind)), len(out))
        out.insert(at, (kind, data))
    return out


def to_bytes(
    record: SeqRecord,
    *,
    circular: bool | None = None,
    flags: int | None = None,
    default_color: str = DEFAULT_COLOR,
    translate_cds: bool = True,
    min_primer_length: int = MIN_PRIMER_LENGTH,
    preserved_packets: Iterable[Packet] | None = None,
) -> bytes:
    """Serialise ``record`` as SnapGene file bytes. See :func:`write_snapgene`."""
    sequence = str(record.seq).upper()
    if not sequence:
        raise ValueError("record has no sequence")

    if circular is None:
        circular = record.annotations.get("topology") == "circular"
    if flags is None:
        flags = record.annotations.get("snapgene_flags")
        if flags is None:
            flags = DEFAULT_DNA_FLAGS
        flags = int(flags) & ~FLAG_CIRCULAR
    flags = int(flags) | (FLAG_CIRCULAR if circular else 0)

    features: list[SeqFeature] = []
    primers: list[SeqFeature] = []
    seen: set[tuple[object, ...]] = set()
    for feature in record.features:
        if feature.type == "source":
            continue
        # Only an annotation that is identical in every respect is a duplicate. A map can
        # legitimately carry two features with the same name, span and type that differ in
        # their qualifiers, and collapsing those loses part of the map.
        key = (
            _label(feature),
            int(feature.location.start),
            int(feature.location.end),
            feature.location.strand,
            feature.type,
            repr(sorted((k, [str(x) for x in v]) for k, v in feature.qualifiers.items())),
        )
        if key in seen:
            continue  # the same annotation carried in twice
        seen.add(key)
        # A primer_bind that came from SnapGene's Features list (it carries a
        # native type) stays a feature; any other primer_bind becomes a primer.
        if feature.type == "primer_bind" and not _first(feature, Q_NATIVE_TYPE):
            span = int(feature.location.end) - int(feature.location.start)
            # A site SnapGene itself wrote is kept whatever its length: a short one is
            # an additional 5' match of a longer primer, not a primer of its own, and
            # dropping it both loses data and breaks the byte-for-byte round trip. The
            # floor exists for primers built in Python, which SnapGene would erase.
            if span < min_primer_length and Q_ATTR_ORDER not in feature.qualifiers:
                continue  # SnapGene would drop it on open
            primers.append(feature)
        else:
            features.append(feature)

    originals: dict[int, bytes] = record.annotations.get("snapgene_original_packets") or {}
    cookie = record.annotations.get("snapgene_cookie") or DEFAULT_COOKIE
    packets: dict[int, bytes] = {
        COOKIE: struct.pack(">8sHHH", b"SnapGene", *cookie),
        # written in the record's own case; SnapGene itself stores lower case
        DNA: struct.pack(">B", flags) + str(record.seq).encode("ascii"),
        FEATURES: _features_xml(
            sequence,
            features,
            default_color,
            translate_cds,
            record.annotations.get("snapgene_features_next_id"),
        ).encode("utf-8"),
    }
    if primers or PRIMERS in originals:
        packets[PRIMERS] = _primers_xml(sequence, primers, record).encode("utf-8")
    edited = any(packets.get(k) != originals.get(k) for k in (DNA, FEATURES, PRIMERS))
    packets[NOTES] = _notes_xml(record, edited).encode("utf-8")

    order = record.annotations.get("snapgene_packet_order")
    if preserved_packets is None:
        unknown = list(record.annotations.get("snapgene_unknown_packets") or [])
    else:
        unknown = list(preserved_packets)
    if not order:
        # a map built from scratch gets the small packets SnapGene always writes
        present = {p.type for p in unknown}
        for kind, data in _DEFAULT_PACKETS.items():
            if kind not in present:
                unknown.append(Packet(kind, data))

    return b"".join(pack_packet(kind, data) for kind, data in _ordered(packets, unknown, order))


def write_snapgene(
    record: SeqRecord,
    target: str | Path | BinaryIO,
    **kwargs: object,
) -> Path | None:
    """Write ``record`` to ``target`` as a SnapGene ``.dna`` file.

    Parameters
    ----------
    record:
        Any :class:`SeqRecord`. Feature colours are taken from the
        ``snapgene_color`` qualifier when present, otherwise ``default_color``.
        A record read with :func:`pysnapgene.read_snapgene` and not edited
        is written back byte for byte.
    target:
        Path or open binary file object.
    circular:
        Override the topology. Defaults to ``record.annotations["topology"]``.
    flags:
        Override the whole DNA flag byte. The circular bit is always set from
        ``circular``.
    default_color:
        Colour for features that carry none.
    translate_cds:
        Mark new translatable features so SnapGene shows their protein. A feature
        is only marked when its span really is a clean reading frame. Features
        read from SnapGene keep their own setting.
    min_primer_length:
        Primer binding sites shorter than this are dropped, matching SnapGene's
        own behaviour on open.
    preserved_packets:
        Packets to write verbatim, by default those kept by
        :func:`pysnapgene.read_snapgene`.

    Returns
    -------
    pathlib.Path or None
        The path written, or ``None`` when a file object was supplied.
    """
    payload = to_bytes(record, **kwargs)  # type: ignore[arg-type]
    if isinstance(target, (str, Path)):
        path = Path(target)
        path.write_bytes(payload)
        return path
    target.write(payload)
    return None
