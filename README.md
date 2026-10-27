# pysnapgene-tools

Python tools for plasmid design: read and **write** SnapGene `.dna` files, codon-optimise
and de-repeat a designed sequence, and simulate Golden Gate assembly. Its only dependency
is Biopython.

[![CI](https://github.com/meng-yx/pysnapgene-writer/actions/workflows/ci.yml/badge.svg)](https://github.com/meng-yx/pysnapgene-writer/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Written for SnapGene 8.2** (tested with SnapGene 8.2.2 on macOS). Files from other
SnapGene versions are read and written the same way, but have not been tested.

| subpackage | what it does |
|---|---|
| [`pysnapgene.writer`](#pysnapgenewriter--snapgene-dna-files) | read and write SnapGene `.dna`, keeping colours, feature types, direction and translation |
| [`pysnapgene.design`](#pysnapgenedesign--codon-choice-and-repeats) | codon-optimise against a usage table; find and remove direct repeats |
| [`pysnapgene.cloning`](#pysnapgenecloning--cloning-in-silico) | type IIS digestion, Golden Gate, and homology assembly (In-Fusion / Gibson / HiFi) |

```python
from pysnapgene import read_snapgene, write_snapgene

record = read_snapgene("plasmid.dna")          # a normal Biopython SeqRecord
record.features[3].qualifiers["snapgene_color"]
# ['#993366']

write_snapgene(record, "copy.dna")             # opens in SnapGene, colours intact
```

## Why

Biopython can parse SnapGene, but:

* **it cannot write it.** `snapgene` does not appear in `Bio.SeqIO._FormatToWriter`, so a
  plasmid designed in Python has to be exported as GenBank and re-saved by hand.
* **its parser is lossy.** Segment colours, the native feature type spelling, the
  translation flag, feature direction, and every packet it does not model (enzyme sets,
  alignments, edit history) are dropped. Read a map and write it back and you get a grey,
  stripped file.

The `design` and `cloning` subpackages exist for the other half of the job. A designed
plasmid is only as good as the checks run on it, and two checks keep coming up: that no
element is spelled the same way twice (or one recombination event deletes what lies
between the copies), and that the construct the vector actually produces is the one
intended — which, with Golden Gate, means assembling it rather than counting bases, since
the overhangs contribute residues of their own at both junctions.

## Install

**Not published yet** — there is no `pysnapgene-tools` on PyPI, so install from source:

```bash
git clone https://github.com/meng-yx/pysnapgene-writer    # repo still has the old name
cd pysnapgene-writer
conda env create -f environment.yml
conda activate pysnapgene-tools
pip install -e ".[dev]"
pytest
```

Or straight into an environment you already have, which is how this is used day to day:

```bash
pip install -e /path/to/the/checkout --no-deps
```

Packaging is ready for release — `pyproject.toml` builds a wheel and `.github/workflows/ci.yml`
publishes on a `v*` tag — but nothing has been uploaded and the name is unclaimed. Publishing
needs, in order: a PyPI account holding the `pysnapgene-tools` name, a trusted publisher (or the
`PYPI_API_TOKEN` secret) on the `release` environment, and a tag. Until then `pip install
pysnapgene-tools` installs nothing, or worse installs whatever someone else registers under that
name, so do not put the command in anyone's instructions.

## `pysnapgene.writer` — SnapGene `.dna` files

### What survives a round trip

A map that is read and written back without edits is intended to come out **byte for byte
identical**, so SnapGene treats the copy exactly like the original; after an edit, only the
edited parts change. Everything else stays as SnapGene wrote it: attribute order, feature
and primer IDs, packet order, and the notes, apart from the *last modified* date.

Measured over a collection of 524 real `.dna` files, byte-for-byte identity currently holds
for 43 of them. The usual remaining difference is the order in which SnapGene groups the
"simplified" copies of a primer's binding sites — it groups them per primer, this package
emits them per site. Everything in the table below is preserved either way.

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

### Building a map from scratch

```python
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pysnapgene import write_snapgene

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

### Three details that decide whether the map looks right

These are the ones that cost time to work out, so they are worth stating plainly.

**1. The double-stranded bit.** The DNA packet begins with a flag byte:

```
0x01 circular | 0x02 double-stranded | 0x04 Dam | 0x08 Dcm | 0x10 EcoKI
```

Leave `0x02` clear and SnapGene opens the file as **single-stranded** DNA — no
complementary strand, no translation. `pysnapgene.writer` sets `0x1e` plus the circular bit by
default, which is what SnapGene writes for a plasmid propagated in *E. coli*.

**2. Colour is per segment and is never inferred.** A file written without a `color`
attribute renders entirely grey. Pass `default_color=` or set `snapgene_color` per feature.

**3. Translation is `translated="1"` on the segment**, together with `codon_start`,
`transl_table` and `translation` qualifiers. That is what ticks *Translate this feature in
Sequence view*.

On that last point, the writer will not mark a new feature as a translated CDS
unless the span really is one — a whole number of codons with no internal stop. Feature boundaries
carried over from another map frequently are not in frame, and marking those as CDS makes
SnapGene display a nonsense protein. Such features keep their colour and position but are
written as `misc_feature`. Features read from a SnapGene file keep their own setting.

### API

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
from pysnapgene import iter_packets
for packet in iter_packets(open("plasmid.dna", "rb")):
    print(f"0x{packet.type:02x} {packet.name:28} {len(packet)} bytes")
```

## `pysnapgene.design` — codon choice and repeats

```python
from pysnapgene.design import CodonUsage, codon_optimize, longest_repeat, break_repeats

usage = CodonUsage.from_csv("ecoli_codon_usage.csv")   # codon, amino_acid, fraction
```

**Codon optimisation.** `codon_optimize(protein, usage, ...)` samples codons in proportion
to their usage and keeps the best of `tries` attempts, or takes the commonest codon
everywhere with `argmax=True`. Two arguments matter more than the rest:

* `forbidden` — motifs to keep out, checked on **both strands**, so passing `["GGTCTC"]`
  excludes `GAGACC` as well. This is how an internal BsaI site is kept out of a Golden Gate
  part.
* `avoid` — sequences the result should resemble as little as possible. *Best* means the
  shortest substring shared with anything in `avoid`, which is what you want when the same
  element appears twice in one construct:

```python
first = codon_optimize("LEVLFQGP", usage, forbidden=["GGTCTC"])
second = codon_optimize("LEVLFQGP", usage, avoid=[first.dna], forbidden=["GGTCTC"])
second.shared        # longest run the two spellings share, in bp
```

**Repeats.** `longest_repeat(sequence, min_length)` returns the longest substring occurring
at least twice, as `(motif, first, second)`, or `None`. Pass `circular=True` for a plasmid,
where a repeat can straddle the origin — a linear search cannot see those.

`break_repeats(sequence, editable, usage, max_length=13)` then removes them with minimal
synonymous edits: one codon at a time, working outwards from the middle of the repeat so it
is split rather than trimmed. `editable` is a list of half-open `(start, end)` spans that
may be changed, each starting on a codon boundary; everything else — a preserved cassette,
a backbone — is left alone. A repeat with no editable position in either copy is returned
in `.unfixable` rather than raising, because the fix is then a design decision.

```python
result = break_repeats(plasmid, editable_spans, usage, max_length=13,
                       forbidden=["GGTCTC"])
result.sequence, result.edits, result.unfixable
```

Why 13: that is roughly where a repeat stops being expected by chance. A ~6 kb plasmid
contains on the order of 7.6 coincidental 10-mers but only ~0.1 coincidental 13-mers, so
13 bp is about the shortest limit reachable by design. `DEFAULT_REPEAT_LIMIT` holds it.

Only **direct** repeats on the same strand are found; inverted repeats are not modelled.

## `pysnapgene.cloning` — cloning in silico

```python
from Bio.Restriction import BsaI
from pysnapgene.cloning import digest, excise, golden_gate, translate_orf

backbone, stuffer = excise(vector, BsaI)        # a vector cut exactly twice
product = golden_gate(vector, insert, BsaI)     # the stuffer replaced by the insert
protein = translate_orf(product, atg_position)  # the fusion, reading frame wrapped
```

`digest` splits on the **top-strand** cut and reports the overhang at each end, so for a 5'
overhang enzyme each fragment's sequence begins with the overhang it shares with the
fragment before it. Overhang length comes from the enzyme, so SapI's three bases work as
well as BsaI's four.

### Homology assembly — In-Fusion, Gibson, HiFi

```python
from pysnapgene.cloning import infusion, ambiguous_homology

assembly = infusion([linearised_backbone, synthesised_fragment])   # 15 bp minimum
assembly.sequence                 # the circular product
assembly.junctions                # (upstream, downstream, homology, start) per join
assembly.locate(1, 240)           # where base 240 of fragment 1 landed in the product
ambiguous_homology(assembly)      # junctions whose homology occurs twice -- unpredictable
```

`homology_assemble` is the same function with the homology length as an argument (Gibson and
HiFi usually want 20–40 bp). Fragments are taken in the order and orientation given — nothing
is reordered or flipped — and a junction that *would* work reverse-complemented says so,
because that is the usual mistake:

```
ValueError: no homology into fragment 1, but its reverse complement shares 20 bp
            -- the fragment is the wrong way round
```

`max_homology` bounds the search; it does not truncate a join. Setting it below the real
overlap makes the junction fail rather than shorten, since a shorter suffix of one fragment is
not a prefix of the next.

`offsets` and `locate` are what make the product annotatable: every feature of every input
fragment can be moved onto the product, which is how an accurate map of an assembly gets
built rather than guessed at.

### Golden Gate

`golden_gate` picks the insert fragment whose overhangs match the vector's stuffer and
raises if none does — which is the check worth having, because an insert ordered against
the wrong overhangs looks perfectly fine until it refuses to assemble. The product is
returned rotated to the backbone's start, and `translate_orf` reads through the origin,
which an ORF in a plasmid frequently does.

The reason to assemble rather than measure: between a tag in the vector and the gene in the
insert, the junctions contribute residues belonging to neither, and a different number at
each end. Only the product tells you how long the linker really is.

3' overhang enzymes and blunt cutters are refused rather than guessed at.

## The .dna file format, briefly

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
* `design` and `cloning` work on plain strings, not `SeqRecord`s, so a map read
  with `read_snapgene` is passed in as `str(record.seq)`.
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
