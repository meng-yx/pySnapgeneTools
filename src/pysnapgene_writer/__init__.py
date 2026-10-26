"""pysnapgene-writer -- a minimalist writer (and reader) for SnapGene ``.dna`` files.

Biopython parses SnapGene but cannot write it, and its parser drops the
SnapGene-only display properties: segment colours, the native feature type
spelling, the translation flag, and every packet it does not model.

This package reads those properties into ordinary ``SeqRecord`` qualifiers and
writes them back, so a map survives a Python round trip.

    >>> from pysnapgene_writer import read_snapgene, write_snapgene
    >>> record = read_snapgene("plasmid.dna")
    >>> record.features[0].qualifiers["snapgene_color"]
    ['#993366']
    >>> write_snapgene(record, "copy.dna")
"""

from __future__ import annotations

from .constants import (
    DEFAULT_COLOR,
    DEFAULT_DNA_FLAGS,
    FLAG_CIRCULAR,
    FLAG_DAM_METHYLATED,
    FLAG_DCM_METHYLATED,
    FLAG_DOUBLE_STRANDED,
    FLAG_ECOKI_METHYLATED,
    Q_COLOR,
    Q_NATIVE_TYPE,
    Q_TRANSLATED,
)
from .packets import Packet, iter_packets, pack_packet
from .reader import SnapGeneError, read_snapgene
from .writer import to_bytes, write_snapgene

__version__ = "0.1.0"

__all__ = [
    "read_snapgene",
    "write_snapgene",
    "to_bytes",
    "SnapGeneError",
    "Packet",
    "iter_packets",
    "pack_packet",
    "Q_COLOR",
    "Q_NATIVE_TYPE",
    "Q_TRANSLATED",
    "DEFAULT_COLOR",
    "DEFAULT_DNA_FLAGS",
    "FLAG_CIRCULAR",
    "FLAG_DOUBLE_STRANDED",
    "FLAG_DAM_METHYLATED",
    "FLAG_DCM_METHYLATED",
    "FLAG_ECOKI_METHYLATED",
    "__version__",
]
