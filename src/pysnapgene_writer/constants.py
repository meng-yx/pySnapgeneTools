"""Flags, qualifier keys and defaults shared by the reader and the writer."""

from __future__ import annotations

__all__ = [
    "FLAG_CIRCULAR",
    "FLAG_DOUBLE_STRANDED",
    "FLAG_DAM_METHYLATED",
    "FLAG_DCM_METHYLATED",
    "FLAG_ECOKI_METHYLATED",
    "DEFAULT_DNA_FLAGS",
    "DEFAULT_COOKIE",
    "Q_COLOR",
    "Q_SEGMENT_COLORS",
    "Q_SEGMENT_NAMES",
    "Q_NATIVE_TYPE",
    "Q_TRANSLATED",
    "Q_ATTR_PREFIX",
    "Q_PRIMER_SEQUENCE",
    "Q_ATTR_ORDER",
    "Q_PREDEF",
    "Q_SITE_XML",
    "RESERVED_QUALIFIERS",
    "DEFAULT_COLOR",
    "MIN_PRIMER_LENGTH",
]

# --- DNA packet flag byte ---------------------------------------------------
# Bit 1 is the one that catches people out: leave it clear and SnapGene shows the
# molecule as single stranded.
FLAG_CIRCULAR = 0x01
FLAG_DOUBLE_STRANDED = 0x02
FLAG_DAM_METHYLATED = 0x04
FLAG_DCM_METHYLATED = 0x08
FLAG_ECOKI_METHYLATED = 0x10

#: Double stranded and methylated, which is what SnapGene itself writes for a
#: plasmid propagated in *E. coli*. The circular bit is added by the writer.
DEFAULT_DNA_FLAGS = (
    FLAG_DOUBLE_STRANDED | FLAG_DAM_METHYLATED | FLAG_DCM_METHYLATED | FLAG_ECOKI_METHYLATED
)

#: ``(sequence type, export version, import version)`` as written by SnapGene 8.
DEFAULT_COOKIE = (1, 15, 20)

# --- qualifier keys used to round-trip SnapGene-only display properties -----
# Per-segment lists (colours, names) are in file order: ascending position, as
# SnapGene stores them, whatever the strand.
Q_COLOR = "snapgene_color"
Q_SEGMENT_COLORS = "snapgene_segment_colors"
Q_SEGMENT_NAMES = "snapgene_segment_names"
Q_NATIVE_TYPE = "snapgene_native_type"
Q_TRANSLATED = "snapgene_translated"
Q_ATTR_PREFIX = "snapgene_attr_"
Q_PRIMER_SEQUENCE = "snapgene_primer_sequence"

# bookkeeping that lets an unedited map be written back byte for byte
Q_ATTR_ORDER = "snapgene_attribute_order"  # attribute order of the <Feature>/<Primer>
Q_PREDEF = "snapgene_predef"  # qualifier names whose values were <V predef=...>
Q_SITE_XML = "snapgene_site_xml"  # a primer binding site's original XML

#: Keys the writer consumes itself rather than emitting as ``<Q>`` elements.
RESERVED_QUALIFIERS = frozenset(
    {
        Q_COLOR,
        Q_SEGMENT_COLORS,
        Q_SEGMENT_NAMES,
        Q_NATIVE_TYPE,
        Q_TRANSLATED,
        Q_PRIMER_SEQUENCE,
        Q_ATTR_ORDER,
        Q_PREDEF,
        Q_SITE_XML,
    }
)

#: SnapGene's own default feature grey.
DEFAULT_COLOR = "#a6acb3"

#: SnapGene discards a primer whose binding site is shorter than this, so the
#: writer does too rather than producing a file that loses primers on open.
MIN_PRIMER_LENGTH = 15
