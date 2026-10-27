# pysnapgene-tools

Python tools for plasmid design: read and **write** SnapGene `.dna` files, codon-optimise and
de-repeat a designed sequence, and simulate the cloning step that produces it. Its only
dependency is Biopython.

[![CI](https://github.com/meng-yx/pysnapgene-writer/actions/workflows/ci.yml/badge.svg)](https://github.com/meng-yx/pysnapgene-writer/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Written for SnapGene 8.2** (tested with SnapGene 8.2.2 on macOS). Files from other SnapGene
versions are read and written the same way, but have not been tested.

```python
from pysnapgene import read_snapgene, write_snapgene

record = read_snapgene("plasmid.dna")          # a normal Biopython SeqRecord
record.features[3].qualifiers["snapgene_color"]
# ['#993366']

write_snapgene(record, "copy.dna")             # opens in SnapGene, colours intact
```

## Why

Biopython can parse SnapGene, but **it cannot write it** — `snapgene` does not appear in
`Bio.SeqIO._FormatToWriter`, so a plasmid designed in Python has to be exported as GenBank and
re-saved by hand — and **its parser is lossy**: segment colours, the native feature type spelling,
the translation flag, feature direction, and every packet it does not model are dropped. Read a
map and write it back and you get a grey, stripped file.

The other two subpackages exist for the rest of the job. A designed plasmid is only as good as
the checks run on it, and two keep coming up: that no element is spelled the same way twice (or
one recombination deletes what lies between the copies), and that the construct the vector
actually produces is the one intended — which, for any assembly, means building it rather than
counting bases, since the junctions contribute sequence of their own.

## The three subpackages

| | what it does | detail |
|---|---|---|
| **`pysnapgene.writer`** | Reads and writes `.dna`, keeping feature colours, native types, direction, translation flags, primers, notes and unrecognised packets. Re-exported at the top level, so `from pysnapgene import read_snapgene` works. An unedited read/write is intended to be byte-for-byte identical; it currently is for 43 of 524 real files tested, and everything listed as preserved survives either way. `set_origin` rotates a circular map to a new first base, splitting whatever now straddles the origin. | [docs/writer.md](docs/writer.md) |
| **`pysnapgene.design`** | `CodonUsage` + `codon_optimize` to encode a protein against a usage table, keeping chosen motifs out on both strands and keeping two copies of one element from sharing a long run. `longest_repeat`, `longest_shared_substring` and `break_repeats` to find direct repeats — origin-crossing ones included — and remove them with minimal synonymous edits inside spans you mark editable. | [docs/design.md](docs/design.md) |
| **`pysnapgene.cloning`** | `digest` and `excise` for type IIS enzymes, taking the overhang length from the enzyme. `golden_gate` to assemble an insert into a destination vector, refusing a mismatched overhang. `infusion` / `homology_assemble` for In-Fusion, Gibson and HiFi, reporting each junction and where every input fragment landed, so the product can be annotated rather than guessed at. `translate_orf` reads a fusion through the origin. | [docs/cloning.md](docs/cloning.md) |

```python
from pysnapgene.design import CodonUsage, codon_optimize, break_repeats
from pysnapgene.cloning import digest, excise, golden_gate, infusion, translate_orf
```

### What it does not do

* Protein and RNA SnapGene files, `.prot`, `.pdna` and project files are out of scope — DNA only.
* Inverted repeats are not modelled, only direct ones.
* 3' overhang enzymes and blunt cutters are refused rather than guessed at.
* Assembly efficiency is not modelled: a four-piece assembly simulates as readily as a two-piece
  one, and only the second is pleasant at the bench.
* `design` and `cloning` work on plain strings, not `SeqRecord`s, so a map read with
  `read_snapgene` is passed in as `str(record.seq)`.
* The `.dna` format is reverse-engineered from public descriptions and from inspecting files,
  validated by round-tripping real SnapGene output rather than against a published spec.

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

## Testing

```bash
pytest                                          # synthetic fixtures, runs anywhere
PYSNAPGENE_TEST_DNA=/path/to/real.dna pytest    # also round-trips a real file
```

The suite builds its own SnapGene files rather than shipping a plasmid, because SnapGene's maps
are copyrighted. Point `PYSNAPGENE_TEST_DNA` at any file SnapGene wrote and an extra test
verifies that the file is written back byte for byte, and that the sequence and every feature
position, strand, colour and translation survive a read/write/read cycle.

`examples/quickstart.ipynb` reads pUC19, adds an annotation and a sequencing primer, and writes a
copy. Download pUC19 into `examples/` first — see [docs/writer.md](docs/writer.md).

## Contributing

Issues and pull requests welcome. Please run `ruff check src tests`, `ruff format src tests` and
`pytest` before opening one.

## Licence

MIT — see [LICENSE](LICENSE).

SnapGene is a trademark of GSL Biotech LLC. This project is not affiliated with or endorsed by
GSL Biotech LLC.
