"""Shared fixtures.

The tests build their own SnapGene files rather than shipping a plasmid, so the
suite runs anywhere. Set ``PYSNAPGENE_TEST_DNA`` to a real ``.dna`` path to also
exercise the parser against a file SnapGene itself wrote.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from pysnapgene_writer import Q_COLOR, Q_TRANSLATED

# a tiny circular "plasmid": a promoter, a clean ORF, a reverse CDS and a primer
SEQUENCE = (
    "TTGACAATTAATCATCGGCTCGTATAATG"  # 29 nt, promoter-ish
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGTAA"  # 39 nt ORF, ends with a stop
    "GGCGCGCCTTAATTAAGCGGCCGCACGCGT"  # 30 nt spacer
    "ATGCATCACCATCACCATCACTAA"  # 24 nt His tag ORF
    "CCTGCAGGTCGACTCTAGAGGATCCCCGGG"  # 30 nt
)


def _feature(start, end, strand, label, ftype, colour, translated=False):
    quals = {"label": [label], Q_COLOR: [colour]}
    if translated:
        quals[Q_TRANSLATED] = ["1"]
    return SeqFeature(SimpleLocation(start, end, strand=strand), type=ftype, qualifiers=quals)


@pytest.fixture
def record() -> SeqRecord:
    rec = SeqRecord(
        Seq(SEQUENCE),
        id="pTest001",
        name="pTest001",
        description="pTest001 synthetic test plasmid",
    )
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    rec.features = [
        _feature(0, 29, 1, "promoter", "promoter", "#ffffff"),
        _feature(29, 68, 1, "GFP_frag", "CDS", "#75c6a9", translated=True),
        _feature(98, 122, 1, "His6", "CDS", "#c7b0e3", translated=True),
        _feature(68, 98, -1, "spacer_rev", "misc_feature", "#d6d6d6"),
        SeqFeature(
            SimpleLocation(122, 145, strand=1),
            type="primer_bind",
            qualifiers={"label": ["seq_fwd"]},
        ),
    ]
    return rec


@pytest.fixture
def real_dna_path() -> Path:
    path = os.environ.get("PYSNAPGENE_TEST_DNA")
    if not path or not Path(path).exists():
        pytest.skip("set PYSNAPGENE_TEST_DNA to a real SnapGene .dna file")
    return Path(path)
