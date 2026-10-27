"""Move a circular map's origin, the way SnapGene's *Set Origin* does.

A plasmid has no natural first base, so a map's origin is a presentation choice. A
simulated product comes out wherever the assembly left it -- Golden Gate rotates to the
backbone's start, In-Fusion to the first fragment's -- which is rarely where anyone wants
to read the map from. Numbering a series of vectors from the same landmark, usually the
origin of replication, makes their coordinates comparable.

Rotating means rewriting every feature, and a feature that straddles the new origin has to
become a two-part location: SnapGene stores that as one wrapping range, but a negative or
past-the-end coordinate cannot be written at all.
"""

from __future__ import annotations

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from .constants import Q_SEGMENT_COLORS, Q_SEGMENT_NAMES, Q_SITE_RANGE, Q_SITE_XML

__all__ = ["set_origin", "feature_start"]


def feature_start(record: SeqRecord, label: str) -> int:
    """The 0-based start of the one feature labelled ``label``.

    Raises if no feature carries the label, or if several do -- an ambiguous landmark
    would silently pick one of them.
    """
    found = [
        f
        for f in record.features
        if str(f.qualifiers.get("label", [""])[0]) == label
    ]
    if not found:
        raise KeyError(f"no feature labelled {label!r}")
    if len(found) > 1:
        starts = ", ".join(str(int(f.location.start)) for f in found)
        raise KeyError(f"{len(found)} features labelled {label!r}, starting at {starts}")
    return int(found[0].location.start)


def _logical_parts(feature: SeqFeature, length: int) -> list[tuple[int, int]]:
    """A feature's spans in 5'->3' order, an origin-crossing pair rejoined.

    The rejoined span runs past ``length``, i.e. in the same unwrapped coordinates
    SnapGene's wrapping ranges imply, so rotating it is plain arithmetic.
    """
    parts = (
        list(feature.location.parts)
        if isinstance(feature.location, CompoundLocation)
        else [feature.location]
    )
    if feature.location.strand == -1:
        parts = list(reversed(parts))  # Biopython lists these 5'->3', i.e. descending
    spans = [(int(p.start), int(p.end)) for p in parts]
    rejoined: list[tuple[int, int]] = []
    for start, end in spans:
        if rejoined and rejoined[-1][1] == length and start == 0:
            rejoined[-1] = (rejoined[-1][0], rejoined[-1][1] + end)
        else:
            rejoined.append((start, end))
    return rejoined


def _rotate_span(start: int, end: int, position: int, length: int) -> list[tuple[int, int]]:
    """One span in the rotated frame, split in two if it now crosses the origin."""
    new_start = (start - position) % length
    new_end = new_start + (end - start)
    if new_end <= length:
        return [(new_start, new_end)]
    return [(new_start, length), (0, new_end - length)]


def set_origin(
    record: SeqRecord,
    position: int | None = None,
    *,
    label: str | None = None,
) -> SeqRecord:
    """Return ``record`` rotated so that one base becomes base 1.

    Parameters
    ----------
    record:
        A circular record. Rotating a linear molecule would change the molecule, so
        that is refused.
    position:
        0-based index of the base to put first.
    label:
        Alternatively, the label of the feature whose start to put first -- the usual
        case being the origin of replication, e.g. ``label="colE1_high_copy"``.

    The record is modified in place and returned. Features keep their qualifiers, types,
    strands and lengths; one that straddles the new origin comes back as a two-part
    location, which the writer stores as a single wrapping range.

    A primer's stored binding site (``snapgene_site_range``, ``snapgene_site_xml``) is
    dropped, because those coordinates describe the old frame. The writer recomputes the
    site from the primer sequence, so nothing is lost but the original XML; the primer
    itself, including its own location, is rotated like any other feature.
    """
    if (position is None) == (label is None):
        raise TypeError("pass either position or label, not both")
    if record.annotations.get("topology") != "circular":
        raise ValueError(
            "set_origin only applies to a circular molecule; rotating a linear one "
            "would join its two ends"
        )
    length = len(record.seq)
    if label is not None:
        position = feature_start(record, label)
    assert position is not None
    if not 0 <= position < length:
        raise ValueError(f"origin {position} is outside the {length} bp sequence")
    if position == 0:
        return record

    rotated = []
    for feature in record.features:
        spans: list[tuple[int, int]] = []
        for start, end in _logical_parts(feature, length):
            spans.extend(_rotate_span(start, end, position, length))
        strand = feature.location.strand
        parts = [SimpleLocation(a, b, strand=strand) for a, b in spans]
        if len(parts) == 1:
            location = parts[0]
        else:
            location = CompoundLocation(list(reversed(parts)) if strand == -1 else parts)
        qualifiers = {k: list(v) for k, v in feature.qualifiers.items()}
        if len(parts) != len(feature.location.parts):
            # per-segment lists no longer line up with the segments
            qualifiers.pop(Q_SEGMENT_COLORS, None)
            qualifiers.pop(Q_SEGMENT_NAMES, None)
        qualifiers.pop(Q_SITE_RANGE, None)
        qualifiers.pop(Q_SITE_XML, None)
        rotated.append(
            SeqFeature(location, type=feature.type, id=feature.id, qualifiers=qualifiers)
        )

    record.seq = record.seq[position:] + record.seq[:position]
    record.features = rotated
    return record
