# pysnapgene-writer

A minimalist Python package to **write** SnapGene `.dna` files. It reads them too, so a
map can be edited in Python and saved back with nothing lost. Its only dependency is
Biopython.

[![CI](https://github.com/meng-yx/pysnapgene-writer/actions/workflows/ci.yml/badge.svg)](https://github.com/meng-yx/pysnapgene-writer/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Written for SnapGene 8.2** (tested with SnapGene 8.2.2 on macOS). Files from other
SnapGene versions are read and written the same way, but have not been tested.

## Why

Biopython can parse SnapGene, but:

* **it cannot write it.** `snapgene` does not appear in `Bio.SeqIO._FormatToWriter`, so a
  plasmid designed in Python has to be exported as GenBank and re-saved by hand.
* **its parser is lossy.** Segment colours, the native feature type spelling, the
  translation flag, and every packet it does not model (enzyme sets, alignments, edit
  history) are dropped. Read a map and write it back and you get a grey, stripped file.

`pysnapgene-writer` writes SnapGene files and round-trips all of that.

```python
from pysnapgene_writer import read_snapgene, write_snapgene

record = read_snapgene("plasmid.dna")          # a normal Biopython SeqRecord
record.features[3].qualifiers["snapgene_color"]
# ['#993366']

write_snapgene(record, "copy.dna")             # opens in SnapGene, colours intact
```

## Install

```bash
pip install pysnapgene-writer
```

From source, with the development environment:

```bash
git clone https://github.com/meng-yx/pysnapgene-writer
cd pysnapgene-writer
conda env create -f environment.yml
conda activate pysnapgene-writer
pip install -e ".[dev]"
pytest
```

## What survives a round trip

A map that is read and written back without edits comes out **byte for byte identical**,
so SnapGene treats the copy exactly like the original. After an edit, only the edited
parts change. Everything else stays as SnapGene wrote it: attribute order, feature and
primer IDs, packet order, and the notes, apart from the *last modified* date.

| | preserved |
|---|---|
| sequence and topology (circular / linear) | yes |
| double-stranded / methylation flags | yes |
| feature positions, strands, multi-segment features | yes |
| features with **no direction** (drawn as boxes) stay without one (`strand=None`) | yes |
| **feature colours**, per segment | yes |
| native SnapGene feature type (`misc_feature_` and friends) | yes |
| feature qualifiers (`<Q>` elements) | yes |
| **translation flag** and the translated protein | yes |
| **segment names** (`-10`, `-35`, signal sequence) and features that cross the origin | yes |
| primers, as primers, with their binding sites and descriptions | yes |
| `primer_bind` *features* (e.g. coloured M13 sites) stay features | yes |
| notes, whole: description, references, UUID, dates | yes |
| unrecognised packets (enzyme sets, alignments, history) | yes, verbatim |

SnapGene-only display properties are stored as ordinary qualifiers, so a `SeqRecord`
stays a `SeqRecord` and still works with the rest of Biopython:

```python
feature.qualifiers["snapgene_color"]        # ['#75c6a9']
feature.qualifiers["snapgene_native_type"]  # ['misc_feature_']
feature.qualifiers["snapgene_translated"]   # ['1'] when SnapGene shows a protein
feature.qualifiers["snapgene_segment_names"]  # ['-10', '', '-35'], ascending
feature.qualifiers["snapgene_attr_<name>"]  # any other <Feature> attribute
```

## Building a map from scratch

```python
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pysnapgene_writer import write_snapgene

record = SeqRecord(Seq("ATG" + "GTGAGCAAGGGCGAGGAG" * 20 + "TAA"),
                   id="pDemo", description="pDemo demo construct")
record.annotations["topology"] = "circular"
record.features.append(
    SeqFeature(SimpleLocation(0, 363, strand=1), type="CDS",
               qualifiers={"label": ["my_orf"], "snapgene_color": ["#75c6a9"]})
)
write_snapgene(record, "pDemo.dna")
```

`examples/quickstart.ipynb` reads pUC19, adds an annotation and a sequencing primer, and
writes `pUC19_modified.dna`. The map is not shipped here, because SnapGene's maps are
copyrighted. Download it from the
[SnapGene pUC19 page](https://www.snapgene.com/plasmids/basic_cloning_vectors/pUC19)
into `examples/` first.

A `primer_bind` feature you add becomes a primer in SnapGene's Primers list, written the
way SnapGene writes one: binding site, melting temperature (an estimate within a few
degrees of SnapGene's; SnapGene recalculates it), and the "simplified" copy of each site.
SnapGene does not draw primers on the map until *Show Primers* is on. To make a
`primer_bind` an ordinary feature instead, set `snapgene_native_type` to `primer_bind`.

A feature's arrow comes from its strand: `1` or `-1` draws an arrow, `None` draws a box.

## Three details that decide whether the map looks right

These are the ones that cost time to work out, so they are worth stating plainly.

**1. The double-stranded bit.** The DNA packet begins with a flag byte:

```
0x01 circular | 0x02 double-stranded | 0x04 Dam | 0x08 Dcm | 0x10 EcoKI
```

Leave `0x02` clear and SnapGene opens the file as **single-stranded** DNA — no
complementary strand, no translation. `pysnapgene-writer` sets `0x1e` plus the circular bit by
default, which is what SnapGene writes for a plasmid propagated in *E. coli*.

**2. Colour is per segment and is never inferred.** A file written without a `color`
attribute renders entirely grey. Pass `default_color=` or set `snapgene_color` per feature.

**3. Translation is `translated="1"` on the segment**, together with `codon_start`,
`transl_table` and `translation` qualifiers. That is what ticks *Translate this feature in
Sequence view*.

On that last point, `pysnapgene-writer` will not mark a new feature as a translated CDS
unless the span really is one — a whole number of codons with no internal stop. Feature boundaries
carried over from another map frequently are not in frame, and marking those as CDS makes
SnapGene display a nonsense protein. Such features keep their colour and position but are
written as `misc_feature`. Features read from a SnapGene file keep their own setting.

## API

```python
read_snapgene(source, *, preserve_unknown_packets=True) -> SeqRecord
```
`source` is a path, raw `bytes`, or an open binary file.

```python
write_snapgene(record, target, *, circular=None, flags=None,
               default_color="#a6acb3", translate_cds=True,
               min_primer_length=15, preserved_packets=None) -> Path | None
to_bytes(record, **same_kwargs) -> bytes
```

* `circular` — defaults to `record.annotations["topology"]`.
* `flags` — override the whole flag byte; the circular bit is always set from `circular`.
* `min_primer_length` — binding sites shorter than this are dropped, because SnapGene
  discards them on open and would otherwise warn that primers had been erased.

The packet layer is public too, for anyone wanting to inspect a file:

```python
from pysnapgene_writer import iter_packets
for packet in iter_packets(open("plasmid.dna", "rb")):
    print(f"0x{packet.type:02x} {packet.name:28} {len(packet)} bytes")
```

## File format, briefly

A `.dna` file is a flat run of TLV packets: one type byte, a big-endian `uint32` length,
then the payload.

| packet | contents |
|---|---|
| `0x09` | cookie: `SnapGene` magic + sequence type + export/import versions |
| `0x00` | flag byte, then the sequence |
| `0x05` | primers, as XML |
| `0x06` | notes, as XML |
| `0x0A` | features, as XML; ranges are **1-based inclusive** |
| others | enzyme sets, alignments, history — preserved verbatim, not parsed |

Primer binding-site coordinates sit one lower than a feature's at both ends; the reader
and writer both account for that.

## Limitations

* Files written from scratch are **minimal**: cookie, DNA, features, primers and notes,
  plus the small default packets SnapGene always writes. They carry no enzyme sets or
  edit history unless those packets came from a file you read (in which case they are
  preserved).
* Editing a map that was read keeps its SnapGene edit history as it was. SnapGene's
  history view will not show the edits made in Python.
* Protein and RNA SnapGene files are not supported; DNA only.
* The format is reverse-engineered from public descriptions and from inspecting files.
  It is validated by round-tripping real SnapGene output, not against a published spec.
* `.prot`, `.pdna` and SnapGene project files are out of scope.

## Testing

```bash
pytest                       # synthetic fixtures, runs anywhere
PYSNAPGENE_TEST_DNA=/path/to/real.dna pytest   # also round-trips a real file
```

The suite builds its own SnapGene files rather than shipping a plasmid. Point
`PYSNAPGENE_TEST_DNA` at any file SnapGene wrote and an extra test verifies that the
file is written back byte for byte, and that the sequence and every feature position,
strand, colour and translation survive a read/write/read cycle.

## Contributing

Issues and pull requests welcome. Please run `ruff check src tests`,
`ruff format src tests` and `pytest` before opening one.

## Licence

MIT — see [LICENSE](LICENSE).

SnapGene is a trademark of GSL Biotech LLC. This project is not affiliated with or
endorsed by GSL Biotech LLC.
