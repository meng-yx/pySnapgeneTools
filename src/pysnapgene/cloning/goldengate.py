"""Type IIS digestion and Golden Gate assembly, in silico.

The reason to simulate rather than trust the design: a Golden Gate junction contributes
residues of its own to the fusion protein, and it contributes a different number at each
end depending on where the cut falls in the reading frame. Counting linker length in the
empty vector therefore gives the wrong answer. Assembling the product and reading the
protein off it gives the right one.

Coordinates follow Biopython: ``digest`` splits on the **top-strand** cut, so for a 5'
overhang enzyme each fragment's top strand begins with the four-base overhang it shares
with the fragment before it.
"""

from __future__ import annotations

from typing import NamedTuple

from Bio.Restriction.Restriction import AbstractCut
from Bio.Seq import Seq

__all__ = [
    "Fragment",
    "digest",
    "excise",
    "golden_gate",
    "translate_orf",
]


class Fragment(NamedTuple):
    """One top-strand fragment of a digest, with the overhangs at either end.

    ``left`` is the overhang this fragment's own top strand starts with; ``right`` is the
    one the next fragment starts with, which this fragment anneals to. Either is empty at
    the end of a linear molecule.
    """

    seq: str
    left: str
    right: str
    start: int

    def __len__(self) -> int:
        return len(self.seq)


def _overhang_length(enzyme: AbstractCut) -> int:
    length = abs(int(enzyme.ovhg))
    if length == 0:
        raise ValueError(f"{enzyme} leaves blunt ends, so fragments cannot be ordered")
    if not enzyme.is_5overhang():
        raise NotImplementedError(f"{enzyme} leaves 3' overhangs, which are not modelled")
    return length


def digest(sequence: str, enzyme: AbstractCut, *, circular: bool = True) -> list[Fragment]:
    """Cut ``sequence`` with ``enzyme`` and return the fragments in order.

    A circular molecule gives one fragment per cut site, the first beginning at the first
    cut; a linear one gives one more, and its outermost ends carry no overhang.
    """
    overhang = _overhang_length(enzyme)
    seq = Seq(sequence.upper())
    cuts = sorted(enzyme.search(seq, linear=not circular))
    if not cuts:
        raise ValueError(f"{enzyme} does not cut this sequence")
    pieces = [str(piece) for piece in enzyme.catalyse(seq, linear=not circular)]

    # Biopython reports cut positions 1-based, and on a circular molecule it returns the
    # fragments rotated: the first one begins at the LAST cut and wraps the origin.
    starts: list[int]
    if circular:
        starts = [cuts[-1] - 1] + [cut - 1 for cut in cuts[:-1]]
    else:
        starts = [0] + [cut - 1 for cut in cuts]
    if len(starts) != len(pieces):  # pragma: no cover - guards a Biopython change
        raise AssertionError(f"{len(pieces)} fragments for {len(starts)} boundaries")

    fragments = []
    for index, (piece, start) in enumerate(zip(pieces, starts)):
        is_first, is_last = index == 0, index == len(pieces) - 1
        if circular:
            left = piece[:overhang]
            right = pieces[(index + 1) % len(pieces)][:overhang]
        else:
            left = "" if is_first else piece[:overhang]
            right = "" if is_last else pieces[index + 1][:overhang]
        fragments.append(Fragment(piece, left, right, start))
    return fragments


def excise(vector: str, enzyme: AbstractCut, *, circular: bool = True) -> tuple[Fragment, Fragment]:
    """Split a destination vector cut exactly twice into ``(backbone, stuffer)``.

    The stuffer is the piece released between the two sites -- the ccdB cassette, in a
    counter-selection vector -- and the backbone is what the insert is ligated into. They
    are told apart by length, the stuffer being the shorter; pass the fragments to
    :func:`digest` yourself if a vector ever has a stuffer longer than its backbone.
    """
    fragments = digest(vector, enzyme, circular=circular)
    if len(fragments) != 2:
        raise ValueError(
            f"a destination vector should be cut twice, giving 2 fragments, not {len(fragments)}"
        )
    first, second = fragments
    return (first, second) if len(first) > len(second) else (second, first)


def golden_gate(
    vector: str,
    insert: str,
    enzyme: AbstractCut,
    *,
    vector_circular: bool = True,
    insert_circular: bool = False,
) -> str:
    """Assemble ``insert`` into ``vector``, returning the circular product.

    The insert is whichever fragment of the digested ``insert`` sequence carries overhangs
    matching the vector's stuffer. The returned sequence starts where the backbone starts,
    i.e. at the overhang 3' of the insert, so it is rotated relative to the input vector.

    Raises ``ValueError`` if no fragment of the insert matches both overhangs -- which is
    the check worth having, because an insert ordered against the wrong overhangs looks
    perfectly fine until it refuses to assemble.
    """
    backbone, stuffer = excise(vector, enzyme, circular=vector_circular)
    candidates = digest(insert, enzyme, circular=insert_circular)
    matches = [f for f in candidates if f.left == stuffer.left and f.right == stuffer.right]
    if not matches:
        offered = ", ".join(f"{f.left}/{f.right}" for f in candidates if f.left or f.right)
        raise ValueError(
            f"no insert fragment carries the vector's overhangs "
            f"{stuffer.left}/{stuffer.right}; the digest offers {offered or 'none'}"
        )
    if len(matches) > 1:
        raise ValueError(f"{len(matches)} insert fragments carry {stuffer.left}/{stuffer.right}")
    return backbone.seq + matches[0].seq


def translate_orf(sequence: str, start: int, *, circular: bool = True, limit: int = 100_000) -> str:
    """Protein encoded from ``start``, up to but not including the first stop codon.

    On a circular sequence the reading frame wraps the origin, which it has to: an ORF in a
    plasmid frequently spans position 0. Raises if no stop codon is reached.
    """
    sequence = sequence.upper()
    if circular:
        turns = -(-limit // max(len(sequence), 1)) + 1
        runway = (sequence * turns)[start : start + limit]
    else:
        runway = sequence[start : start + limit]
    protein = str(Seq(runway[: len(runway) // 3 * 3]).translate())
    stop = protein.find("*")
    if stop < 0:
        raise ValueError(f"no stop codon within {limit} bases of position {start}")
    return protein[:stop]
