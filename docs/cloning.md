# `pysnapgene.cloning` — cloning in silico

Type IIS digestion, Golden Gate, and homology assembly (In-Fusion, Gibson, HiFi), so a designed
vector can be checked against the construct it is meant to produce rather than reasoned about.

[← back to the README](../README.md)

## Digestion

```python
from Bio.Restriction import BsaI
from pysnapgene.cloning import digest, excise

fragments = digest(sequence, BsaI, circular=True)
backbone, stuffer = excise(vector, BsaI)        # a vector cut exactly twice
```

`digest` splits on the **top-strand** cut and reports the overhang at each end, so for a 5'
overhang enzyme each fragment's sequence begins with the overhang it shares with the fragment
before it. Overhang length comes from the enzyme, so SapI's three bases work as well as BsaI's
four. `Fragment.start` is the fragment's own offset in the input — note that Biopython returns a
circular digest rotated so the first fragment begins at the *last* cut, which `digest` accounts
for.

`excise` separates a destination vector into `(backbone, stuffer)`, telling them apart by length.
3' overhang enzymes and blunt cutters are refused rather than guessed at.

## Golden Gate

```python
from pysnapgene.cloning import golden_gate, translate_orf

product = golden_gate(vector, insert, BsaI)     # the stuffer replaced by the insert
protein = translate_orf(product, atg_position)  # the fusion, reading frame wrapped
```

`golden_gate` picks the insert fragment whose overhangs match the vector's stuffer and raises if
none does — which is the check worth having, because an insert ordered against the wrong overhangs
looks perfectly fine until it refuses to assemble. The product is returned rotated to the
backbone's start, so a position `p` in the vector moves to `(p - backbone.start) % len(vector)`.
`translate_orf` reads through the origin by default, which an ORF in a plasmid frequently needs.

The reason to assemble rather than measure: between a tag in the vector and the gene in the
insert, the junctions contribute residues belonging to neither, and a different number at each
end. Only the product tells you how long the linker really is.

## Homology assembly — In-Fusion, Gibson, HiFi

All of these join fragments that share terminal sequence identity, differing only in the homology
length they want, so simulating one is the same operation.

```python
from pysnapgene.cloning import infusion, ambiguous_homology

assembly = infusion([linearised_backbone, synthesised_fragment])   # 15 bp minimum
assembly.sequence                 # the circular product
assembly.junctions                # (upstream, downstream, homology, start) per join
assembly.locate(1, 240)           # where base 240 of fragment 1 landed in the product
ambiguous_homology(assembly)      # junctions whose homology occurs twice — unpredictable
```

`homology_assemble` is the same function with the homology length as an argument (Gibson and HiFi
usually want 20–40 bp; `INFUSION_MIN_HOMOLOGY` is 15).

Fragments are taken in the order and orientation given — nothing is reordered or flipped — and a
junction that *would* work reverse-complemented says so, because that is the usual mistake:

```
ValueError: no homology into fragment 1, but its reverse complement shares 20 bp
            -- the fragment is the wrong way round
```

`max_homology` bounds the search; it does not truncate a join. Setting it below the real overlap
makes the junction fail rather than shorten, since a shorter suffix of one fragment is not a
prefix of the next.

`offsets` and `locate` are what make the product annotatable: every feature of every input
fragment can be moved onto the product, which is how an accurate map of an assembly gets built
rather than guessed at.

Worth checking before ordering anything: the product is only predictable while each homology
sequence is unique, and `ambiguous_homology` says whether one occurs twice.

## Limitations

* 3' overhang enzymes and blunt cutters are refused.
* `excise` tells the backbone from the stuffer by length, so a vector whose stuffer is longer than
  its backbone needs `digest` and a choice made by hand.
* Homology assembly does not model efficiency: a four-piece assembly simulates as readily as a
  two-piece one, and only the second is pleasant at the bench.
* Everything here works on plain strings, not `SeqRecord`s: pass a map read with `read_snapgene`
  in as `str(record.seq)`.
