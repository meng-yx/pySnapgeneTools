# `pysnapgene.design` — codon choice and repeats

Codon-optimise a protein against a usage table, and find and remove the direct repeats that let
a designed plasmid recombine with itself.

[← back to the README](../README.md)

```python
from pysnapgene.design import CodonUsage, codon_optimize, longest_repeat, break_repeats

usage = CodonUsage.from_csv("ecoli_codon_usage.csv")   # codon, amino_acid, fraction
```

## Codon optimisation

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

* Only **direct** repeats on the same strand are found; inverted repeats are not modelled.
* The repeat threshold assumes roughly uniform base composition. A match inside a homopolymer or
  other low-complexity run is far likelier by chance than the 13 bp figure implies, and slipped-
  strand mispairing is a different mechanism from homologous recombination.
* Everything here works on plain strings, not `SeqRecord`s: pass a map read with `read_snapgene`
  in as `str(record.seq)`.
