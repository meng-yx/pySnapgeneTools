# `pysnapgene.design` — protein to orderable DNA

Turn a parts list into one DNA sequence, deterministically; score it the way a synthesis vendor
does; and find and remove the repeats that let a designed plasmid recombine with itself.

[← back to the README](../README.md)

```python
from pysnapgene.design import CodonUsage, build_dna, protein, frozen

usage = CodonUsage.from_csv("ecoli_codon_usage.csv")   # codon, amino_acid, fraction
```

## `build_dna`: the whole step, in one call

```python
design = build_dna(
    [
        frozen("overhang", "AATTGGTCTCAAGGA"),      # kept verbatim
        protein("tag_a", "VTGYRLFEEIL"),
        frozen("ccdB", cassette),                   # its own BsaI sites are left alone
        protein("tag_b", "VTGYRLFEEIL", unlike=["tag_a"]),
    ],
    usage,
    forbidden=["GGTCTC"],
)
design.sequence, design.part("tag_b"), design.spans, design.issues
print(design.report())
```

**Nothing is sampled.** Each codon is the first in a total order — fewest bases shared with the
copies named in `unlike`, then fewest bases repeated from what is already built, then commonest
codon, then alphabetical — so the same parts give the same DNA on any machine. There is no seed
to record and no run to repeat before ordering.

Three arguments carry the design decisions:

* `protein(name, peptide, unlike=[...])` — `unlike` names earlier parts this one should share as
  little as possible with. Use it for the second copy of anything appearing twice. It is the
  difference between a construct a vendor quotes as standard and one it quotes as complex: see
  *Why two copies of one protein are the hard case* below.
* `frozen(name, dna)` — a cassette, overhang, promoter or backbone kept exactly as given. Frozen
  parts are **scanned** like everything else, so a repeat against a backbone is still found, but
  no edit may touch them. `forbidden` is not applied to them either: a ccdB cassette flanked by
  its own BsaI sites is the normal case, and rewriting it would be wrong.
* `rules` / `ignore` — which vendor rules to hold the result to, and which to skip. After
  encoding, `build_dna` calls `fix_complexity` until the rules pass or nothing more can be done.

What it will not do: change the protein, reorder parts, or insert or delete a base. A construct
whose complexity is architectural comes back still flagged, with `issues` saying which rule — the
remedy for that is a design decision, not a codon.

## Checking a sequence the way a vendor does

```python
from pysnapgene.design import complexity_report, summarise, TWIST_GENE_FRAGMENT

for issue in complexity_report(sequence, TWIST_GENE_FRAGMENT):
    print(issue)       # homology at 2-515 and 1005-1518: 514 bp region 84 % identical ...
print(summarise(sequence))
```

`Rules` is a `NamedTuple` of thresholds, so one knob is `TWIST_GENE_FRAGMENT._replace(
max_repeat_density=0.12)`. The defaults model the warnings Twist returns on gene fragments:
length 300–5000 bp, no more than 15 % of the sequence inside a repeat of 9 bp or longer, no
inexact repeat over 70 % identical across 100 bp or with a duplex Tm at or above 68 °C, overall
GC 25–65 % (10–90 % in a 50 bp window), and homopolymers of 13 or less. They are an
approximation: the vendor scores with a model it does not publish, so a sequence that passes here
can still be quoted as complex. What they are reliable for is **comparison** — a sequence that
scores better here is better by the measures the vendor names.

`score(sequence, rules)` collapses the report to one number (0 is clean) by summing how far each
issue is past its limit, which is what `fix_complexity` minimises. The individual measures are
public too: `repeat_density`, `exact_repeats`, `homologous_regions`, `duplex_tm`, `gc_fraction`,
`homopolymers`.

## Why two copies of one protein are the hard case

Two independently optimised spellings of one 157 aa protein share no run longer than 12 bp and
are still **79 % identical base for base**. An exact-repeat search sees nothing; a vendor sees one
long repeat region, and a cell sees a recombination substrate.

`identity(a, b)` measures it, and `homologous_regions` finds such regions without an alignment, by
scanning every offset that carries a shared k-mer seed. `unlike` drives it down — but only so far:

```python
identity_floor("VFTLEDFVGDWEQ...", usage)     # 0.61
```

`identity_floor` is the lowest identity two spellings of that protein could possibly have. Met and
Trp have one codon each, and in most families the first two bases are fixed by the amino acid, so
a typical floor is around 60 %, not 0. Compute it before trying to drive two copies further apart:
if the measured identity is already at the floor, no codon choice will improve it.

## Repairing a sequence you already have

```python
fixed = fix_complexity(sequence, editable_spans, usage, forbidden=["GGTCTC"])
fixed.sequence, fixed.edits, fixed.issues
```

Each round takes the worst remaining issue, walks the codons inside it and inside its partner
region, and takes the **first** swap that lowers the score — first in a fixed order: by position,
then by codon usage, then alphabetically. First rather than best keeps the edit count down and
the run time proportional to the number of edits, and because the order is total the result is
still reproducible. An issue no single codon can improve is reported, not raised.

`editable` lists half-open `(start, end)` spans that may change, each starting on a codon
boundary; everything outside them is protected.

## Codon optimisation, one element at a time

`codon_optimize` predates `build_dna` and is still the direct way to spell a single element.
`codon_optimize(protein, usage, ...)` samples codons in proportion to their usage and keeps the
best of `tries` attempts, or takes the commonest codon everywhere with `argmax=True`. Two
arguments matter more than the rest:

* `forbidden` — motifs to keep out, checked on **both strands**, so passing `["GGTCTC"]` excludes
  `GAGACC` as well. This is how an internal BsaI site is kept out of a Golden Gate part.
* `avoid` — sequences the result should resemble as little as possible. *Best* means the shortest
  substring shared with anything in `avoid`, which is what you want when the same element appears
  twice in one construct:

```python
first = codon_optimize("LEVLFQGP", usage, forbidden=["GGTCTC"])
second = codon_optimize("LEVLFQGP", usage, avoid=[first.dna], forbidden=["GGTCTC"])
second.shared        # longest run the two spellings share, in bp
```

`CodonUsage` also answers the ordinary questions directly: `codons(aa, min_fraction=…)` commonest
first, `best(aa)`, `fraction(codon)`, and `synonyms(codon)`. A `min_fraction` that would exclude
every synonym still returns the commonest one, because a protein has to be encodable.

## Repeats

Why this matters: a plasmid carrying the same sequence twice can recombine across the two copies
and delete everything between them. When what sits between them is the counter-selection marker,
the deletion product survives selection and grows like a real clone, so the repeat has to be
designed out rather than screened out.

`longest_repeat(sequence, min_length)` returns the longest substring occurring at least twice, as
`(motif, first, second)`, or `None`. Pass `circular=True` for a plasmid, where a repeat can
straddle the origin — a linear search cannot see those.

`longest_shared_substring(a, b)` answers the related but different question: the longest sequence
shared between two *regions*, with its position in each. "Is the same sequence present on both
sides of a cassette" is not the same question as "does this plasmid contain a repeat", and only
the first predicts a counter-selection escaper. It binary-searches the length over k-mer sets, so
it works on two 6 kb arms where a quadratic dynamic program does not.

`break_repeats(sequence, editable, usage, max_length=13)` removes repeats with minimal synonymous
edits: one codon at a time, working outwards from the middle of the repeat so it is split rather
than trimmed. `editable` is a list of half-open `(start, end)` spans that may be changed, each
starting on a codon boundary; everything else — a preserved cassette, a backbone — is left alone.
A repeat with no editable position in either copy is returned in `.unfixable` rather than raising,
because the fix is then a design decision.

```python
result = break_repeats(plasmid, editable_spans, usage, max_length=13,
                       forbidden=["GGTCTC"])
result.sequence, result.edits, result.unfixable
```

Why 13: that is roughly where a repeat stops being expected by chance. A ~6 kb plasmid contains on
the order of 7.6 coincidental 10-mers but only ~0.1 coincidental 13-mers, so 13 bp is about the
shortest limit reachable by design. `DEFAULT_REPEAT_LIMIT` holds it.

## Limitations

* Only **direct** repeats on the same strand are found; inverted repeats and hairpins are not
  modelled, by either the repeat search or the complexity rules.
* `homologous_regions` is **ungapped** — right for two spellings of one protein, which stay in
  register, wrong for two homologous genes that do not — and only examines offsets carrying an
  exact seed match, so homology spread evenly enough to leave no 9-mer is missed. Lower `seed` to
  find it, at a cost in speed.
* The vendor rules are read off a vendor's warning messages and public guidance, not from its
  screener. Treat a clean report as "better than the version that was flagged", not as a quote.
* The repeat threshold assumes roughly uniform base composition. A match inside a homopolymer or
  other low-complexity run is far likelier by chance than the 13 bp figure implies, and slipped-
  strand mispairing is a different mechanism from homologous recombination.
* Everything here works on plain strings, not `SeqRecord`s: pass a map read with `read_snapgene`
  in as `str(record.seq)`.
