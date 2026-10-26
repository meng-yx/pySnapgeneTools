"""Homology-based assembly in silico: In-Fusion, Gibson, HiFi, SLiCE.

All of these join fragments that share terminal sequence identity, differing only in the
homology length they want and in the enzymology. Simulating one is therefore the same
operation: find the overlap at each junction, and splice it out of one side.

The reason to simulate rather than reason it through: the product is what you will sequence,
and it is not quite any of the inputs. A junction that regenerates a restriction site, a
homology arm that lands inside a homopolymer, a fragment that needs reverse-complementing --
each is obvious in the assembled sequence and easy to miss in the design.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

from Bio.Seq import Seq

__all__ = [
    "Junction",
    "Assembly",
    "homology_assemble",
    "infusion",
    "ambiguous_homology",
    "INFUSION_MIN_HOMOLOGY",
]

#: In-Fusion HD's documented minimum. Gibson and HiFi usually want 20-40 bp.
INFUSION_MIN_HOMOLOGY = 15


class Junction(NamedTuple):
    """One homology join, and where it ended up in the product."""

    upstream: int
    """Index of the fragment on the 5' side."""
    downstream: int
    """Index of the fragment on the 3' side."""
    homology: str
    """The shared sequence, present once in the product."""
    start: int
    """Where that sequence begins in the product."""

    def __len__(self) -> int:
        return len(self.homology)


class Assembly(NamedTuple):
    """The product, and how the fragments map onto it."""

    sequence: str
    junctions: list[Junction]
    offsets: list[int]
    """Where each fragment's first base lands in the product. On a circular product the last
    fragment may run past the end and wrap, so a span is read modulo ``len(sequence)``."""
    circular: bool

    def __len__(self) -> int:
        return len(self.sequence)

    def locate(self, fragment: int, position: int) -> int:
        """Where ``position`` within fragment ``fragment`` lands in the product."""
        return (self.offsets[fragment] + position) % len(self.sequence)


def _overlap(left: str, right: str, minimum: int, maximum: int) -> str | None:
    """Longest suffix of ``left`` that is also a prefix of ``right``, within bounds."""
    widest = min(maximum, len(left), len(right))
    for size in range(widest, minimum - 1, -1):
        if left[-size:] == right[:size]:
            return right[:size]
    return None


def _describe_failure(left: str, right: str, index: int, minimum: int, maximum: int) -> str:
    """Say why a junction failed, and whether flipping the fragment would fix it."""
    flipped = str(Seq(right).reverse_complement())
    hint = _overlap(left, flipped, minimum, maximum)
    if hint:
        return (
            f"no homology into fragment {index}, but its reverse complement shares "
            f"{len(hint)} bp -- the fragment is the wrong way round"
        )
    best = _overlap(left, right, 1, maximum)
    if best:
        return (
            f"fragment {index} shares only {len(best)} bp ({best}), fewer than the "
            f"{minimum} bp required"
        )
    return f"fragment {index} shares no terminal sequence with the fragment before it"


def homology_assemble(
    fragments: Sequence[str | Seq],
    *,
    circular: bool = True,
    min_homology: int = INFUSION_MIN_HOMOLOGY,
    max_homology: int = 120,
) -> Assembly:
    """Join ``fragments`` on their terminal homology, in the order given.

    Each fragment must share at least ``min_homology`` bases between its 5' end and the 3' end
    of the fragment before it; with ``circular`` the last must likewise close onto the first.
    The shared sequence appears **once** in the product, so the product is shorter than the sum
    of the fragments by the total homology.

    Fragments are taken in the order given and in the orientation given -- nothing is reordered
    or flipped. A junction that would work reverse-complemented is reported as such, since that
    is the usual mistake.

    Raises ``ValueError`` if a junction has no usable homology.
    """
    pieces = [str(f).upper() for f in fragments]
    if not pieces:
        raise ValueError("no fragments to assemble")
    if len(pieces) == 1 and not circular:
        return Assembly(pieces[0], [], [0], False)

    product = pieces[0]
    offsets = [0]
    junctions: list[Junction] = []
    for index in range(1, len(pieces)):
        homology = _overlap(product, pieces[index], min_homology, max_homology)
        if homology is None:
            raise ValueError(
                _describe_failure(product, pieces[index], index, min_homology, max_homology)
            )
        offsets.append(len(product) - len(homology))
        junctions.append(Junction(index - 1, index, homology, len(product) - len(homology)))
        product += pieces[index][len(homology) :]

    if circular:
        # the last fragment's 3' end must duplicate the first fragment's 5' start
        closing = _overlap(product, pieces[0], min_homology, max_homology)
        if closing is None:
            raise ValueError(_describe_failure(product, pieces[0], 0, min_homology, max_homology))
        product = product[: len(product) - len(closing)]
        # the closing homology physically sits at the start of the product: it is fragment 0's
        # own prefix, and the duplicate was trimmed off the far end
        junctions.append(Junction(len(pieces) - 1, 0, closing, 0))
    return Assembly(product, junctions, offsets, circular)


def infusion(
    fragments: Sequence[str | Seq],
    *,
    circular: bool = True,
    min_homology: int = INFUSION_MIN_HOMOLOGY,
    max_homology: int = 120,
) -> Assembly:
    """:func:`homology_assemble` with In-Fusion HD's 15 bp minimum. See it for the details."""
    return homology_assemble(
        fragments, circular=circular, min_homology=min_homology, max_homology=max_homology
    )


def ambiguous_homology(assembly: Assembly) -> list[Junction]:
    """Junctions whose homology occurs more than once in the product.

    Worth checking before ordering: the assembly is only predictable while each homology
    sequence is unique, and a repeat gives the enzyme a second place to join.
    """
    sequence = assembly.sequence
    doubled = sequence + sequence if assembly.circular else sequence
    repeated = []
    for junction in assembly.junctions:
        count = 0
        start = doubled.find(junction.homology)
        while 0 <= start < len(sequence):
            count += 1
            start = doubled.find(junction.homology, start + 1)
        if count > 1:
            repeated.append(junction)
    return repeated
