"""Write a Biopython :class:`SeqRecord` out as a SnapGene ``.dna`` file.

Biopython has no SnapGene writer -- ``snapgene`` does not appear in
``Bio.SeqIO._FormatToWriter`` -- so this module supplies one.

Three details decide whether the resulting map looks right when SnapGene opens it:

* the **double-stranded bit** in the DNA packet flag byte. Leave it clear and the
  molecule is shown as single stranded;
* a **colour on every segment**. SnapGene does not invent one, so a file written
  without the attribute is uniformly grey;
* ``translated="1"`` plus the usual CDS qualifiers, which is what ticks
  *Translate this feature in Sequence view*.
"""

from __future__ import annotations

import struct
import uuid
from collections.abc import Iterable, Sequence
from datetime import date, datetime, timezone
from pathlib import Path
from typing import BinaryIO
from xml.sax.saxutils import escape, quoteattr

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .constants import (
    DEFAULT_COLOR,
    DEFAULT_COOKIE,
    DEFAULT_DNA_FLAGS,
    FLAG_CIRCULAR,
    MIN_PRIMER_LENGTH,
    Q_ATTR_PREFIX,
    Q_COLOR,
    Q_NATIVE_TYPE,
    Q_SEGMENT_COLORS,
    Q_TRANSLATED,
    RESERVED_QUALIFIERS,
)
from .packets import COOKIE, DNA, FEATURES, NOTES, PRIMERS, Packet, pack_packet

__all__ = ["write_snapgene", "to_bytes"]

#: Feature types SnapGene will offer a translation for.
_TRANSLATABLE = {"CDS", "gene", "mat_peptide", "sig_peptide", "transit_peptide"}


def _label(feature: SeqFeature) -> str:
    values = feature.qualifiers.get("label")
    if values:
        return str(values[0])
    return feature.type or "feature"


def _first(feature: SeqFeature, key: str) -> str | None:
    values = feature.qualifiers.get(key)
    return str(values[0]) if values else None


def _segment_colors(feature: SeqFeature, n_parts: int, fallback: str) -> list[str]:
    explicit = feature.qualifiers.get(Q_SEGMENT_COLORS)
    if explicit and len(explicit) == n_parts:
        return [str(c) for c in explicit]
    single = _first(feature, Q_COLOR) or fallback
    return [single] * n_parts


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
    parts = (
        feature.location.parts
        if isinstance(feature.location, CompoundLocation)
        else [feature.location]
    )
    nucleotides = "".join(sequence[int(p.start) : int(p.end)] for p in parts)
    if feature.location.strand == -1:
        nucleotides = str(Seq(nucleotides).reverse_complement())
    return nucleotides


def _features_xml(
    sequence: str, features: Sequence[SeqFeature], default_color: str, translate: bool
) -> str:
    chunks: list[str] = []
    for index, feature in enumerate(features):
        label = _label(feature)
        native_type = _first(feature, Q_NATIVE_TYPE) or feature.type or "misc_feature"
        parts = (
            feature.location.parts
            if isinstance(feature.location, CompoundLocation)
            else [feature.location]
        )
        colors = _segment_colors(feature, len(parts), default_color)

        as_cds = False
        protein = ""
        wants_translation = (
            translate
            and native_type.rstrip("_") in _TRANSLATABLE
            or _first(feature, Q_TRANSLATED) == "1"
        )
        if wants_translation:
            as_cds, protein = _is_clean_orf(_feature_nt(sequence, feature))
            if not as_cds and native_type.rstrip("_") in _TRANSLATABLE:
                native_type = "misc_feature"

        segments = "".join(
            '<Segment range="{start}-{end}" color="{color}" type="standard"{tr}/>'.format(
                start=int(part.start) + 1,
                end=int(part.end),
                color=color,
                tr=' translated="1"' if as_cds else "",
            )
            for part, color in zip(parts, colors)
        )

        qualifier_xml: list[str] = []
        for key, values in feature.qualifiers.items():
            if key in RESERVED_QUALIFIERS or key.startswith(Q_ATTR_PREFIX):
                continue
            if key in ("codon_start", "transl_table", "translation") and as_cds:
                continue  # rewritten below from the real sequence
            rendered = "".join(f"<V text={quoteattr(str(v))}/>" for v in values)
            qualifier_xml.append(f"<Q name={quoteattr(str(key))}>{rendered}</Q>")
        if as_cds:
            qualifier_xml.append('<Q name="codon_start"><V int="1"/></Q>')
            qualifier_xml.append('<Q name="transl_table"><V int="1"/></Q>')
            qualifier_xml.append(f'<Q name="translation"><V text={quoteattr(protein)}/></Q>')

        # Build the attribute set as a mapping so a preserved attribute can never
        # collide with one we emit ourselves and produce invalid XML.
        attributes: dict[str, str] = {
            "recentID": str(index),
            "name": label,
            "type": native_type,
            "directionality": "2" if feature.location.strand == -1 else "1",
            "allowSegmentOverlaps": "0",
            "consecutiveTranslationNumbering": "1",
        }
        for key, values in feature.qualifiers.items():
            if not key.startswith(Q_ATTR_PREFIX) or not values:
                continue
            attribute = key[len(Q_ATTR_PREFIX) :]
            if attribute in ("recentID", "name", "type", "directionality", "translated"):
                continue
            attributes[attribute] = str(values[0])
        if as_cds and "translationMW" not in attributes:
            try:
                from Bio.SeqUtils import molecular_weight

                mass = molecular_weight(protein.replace("*", ""), seq_type="protein")
                attributes["translationMW"] = f"{mass:.2f}"
            except Exception:  # pragma: no cover - MW is cosmetic
                pass
        elif not as_cds:
            attributes.pop("translationMW", None)

        rendered_attributes = " ".join(
            f"{name}={quoteattr(value)}" for name, value in attributes.items()
        )
        chunks.append(
            f"<Feature {rendered_attributes}>{segments}{''.join(qualifier_xml)}</Feature>"
        )
    return (
        f'<?xml version="1.0"?><Features nextValidID="{len(chunks)}">{"".join(chunks)}</Features>'
    )


def _primers_xml(sequence: str, primers: Sequence[SeqFeature]) -> str:
    grouped: dict[tuple[str, str], list[tuple[int, int, int]]] = {}
    for feature in primers:
        start, end = int(feature.location.start), int(feature.location.end)
        bases = _first(feature, "snapgene_primer_sequence")
        if not bases:
            bases = sequence[start:end]
            if feature.location.strand == -1:
                bases = str(Seq(bases).reverse_complement())
        bound = 1 if feature.location.strand == -1 else 0
        grouped.setdefault((_label(feature), bases), []).append((start, end, bound))

    chunks: list[str] = []
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    for index, ((label, bases), sites) in enumerate(grouped.items()):
        binding = "".join(
            f'<BindingSite location="{start}-{end - 1}" boundStrand="{bound}" '
            f"annealedBases={quoteattr(bases)}/>"
            for start, end, bound in sorted(sites)
        )
        chunks.append(
            f"<Primer recentID={quoteattr(str(index))} name={quoteattr(label)} "
            f'sequence={quoteattr(bases)} dateAdded="{stamp}">{binding}</Primer>'
        )
    params = (
        '<HybridizationParams minContinuousMatchLen="10" allowMismatch="1" '
        'minMeltingTemperature="40" showAdditionalFivePrimeMatches="1" '
        'minimumFivePrimeAnnealing="15"/>'
    )
    return (
        '<?xml version="1.0"?>'
        f'<Primers nextValidID="{len(chunks)}">{params}{"".join(chunks)}</Primers>'
    )


def _notes_xml(record: SeqRecord) -> str:
    stored = record.annotations.get("snapgene_notes", {}) or {}
    today = date.today().strftime("%Y.%m.%d")
    fields = {
        "UUID": stored.get("UUID") or str(uuid.uuid4()),
        "Type": stored.get("Type", "Synthetic"),
        "Created": stored.get("Created", today),
        "LastModified": today,
        "AccessionNumber": record.id or stored.get("AccessionNumber", "") or "sequence",
        "Comments": record.description or stored.get("Comments", "") or record.id or "",
    }
    body = "".join(f"<{k}>{escape(str(v))}</{k}>" for k, v in fields.items())
    return f"<Notes>{body}</Notes>"


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
    seen: set[tuple[str, int, int, int | None]] = set()
    for feature in record.features:
        if feature.type == "source":
            continue
        key = (
            _label(feature),
            int(feature.location.start),
            int(feature.location.end),
            feature.location.strand,
        )
        if key in seen:
            continue  # the same annotation carried in twice
        seen.add(key)
        if feature.type == "primer_bind":
            span = int(feature.location.end) - int(feature.location.start)
            if span < min_primer_length:
                continue  # SnapGene would drop it on open
            primers.append(feature)
        else:
            features.append(feature)

    cookie = record.annotations.get("snapgene_cookie") or DEFAULT_COOKIE
    out = bytearray()
    out += pack_packet(COOKIE, struct.pack(">8sHHH", b"SnapGene", *cookie))
    out += pack_packet(DNA, struct.pack(">B", flags) + sequence.encode("ascii"))
    out += pack_packet(
        FEATURES,
        _features_xml(sequence, features, default_color, translate_cds).encode("utf-8"),
    )
    if primers:
        out += pack_packet(PRIMERS, _primers_xml(sequence, primers).encode("utf-8"))
    out += pack_packet(NOTES, _notes_xml(record).encode("utf-8"))

    if preserved_packets is None:
        preserved_packets = record.annotations.get("snapgene_unknown_packets") or []
    for packet in preserved_packets:
        out += pack_packet(packet.type, packet.data)
    return bytes(out)


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
        Mark translatable features so SnapGene shows their protein. A feature is
        only marked when its span really is a clean reading frame.
    min_primer_length:
        Primer binding sites shorter than this are dropped, matching SnapGene's
        own behaviour on open.
    preserved_packets:
        Extra packets to append verbatim, by default those kept by
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
