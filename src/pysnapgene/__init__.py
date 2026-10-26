"""pysnapgene-tools -- SnapGene files, sequence design and cloning simulation.

Three subpackages, each usable on its own:

``pysnapgene.writer``
    Read and write SnapGene ``.dna`` files. Biopython parses SnapGene but cannot
    write it, and its parser drops the SnapGene-only display properties: segment
    colours, the native feature type spelling, the translation flag, feature
    direction, and every packet it does not model. This one keeps all of it.

``pysnapgene.design``
    Codon optimisation against a usage table, and the repeat analysis that keeps a
    designed plasmid from recombining with itself: find the longest repeat, then
    remove it with minimal synonymous edits.

``pysnapgene.cloning``
    Type IIS digestion and Golden Gate assembly in silico, so a designed vector can
    be checked against the construct it is meant to produce.

    >>> from pysnapgene import read_snapgene, write_snapgene
    >>> record = read_snapgene("plasmid.dna")
    >>> record.features[0].qualifiers["snapgene_color"]
    ['#993366']
    >>> write_snapgene(record, "copy.dna")
"""

from __future__ import annotations

from .writer import (
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
    Packet,
    SnapGeneError,
    iter_packets,
    pack_packet,
    read_snapgene,
    to_bytes,
    write_snapgene,
)

__version__ = "0.2.0"

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
