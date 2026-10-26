"""The SnapGene packet layer.

A ``.dna`` file is a flat sequence of TLV packets. Each one is

* 1 byte  -- packet type
* 4 bytes -- big-endian unsigned length of the payload
* N bytes -- the payload

Nothing else wraps them: the file is simply packets end to end. This module is the
only place that knows about that layout.
"""

from __future__ import annotations

import struct
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import BinaryIO

__all__ = [
    "Packet",
    "COOKIE",
    "DNA",
    "PRIMERS",
    "NOTES",
    "FEATURES",
    "PACKET_NAMES",
    "iter_packets",
    "pack_packet",
    "pack_packets",
]

#: Packet types this library understands. Any other type is preserved verbatim.
COOKIE = 0x09
DNA = 0x00
PRIMERS = 0x05
NOTES = 0x06
FEATURES = 0x0A

#: Human-readable names, for diagnostics only.
PACKET_NAMES = {
    0x00: "DNA",
    0x01: "compressed DNA",
    0x02: "unknown",
    0x03: "restriction enzymes",
    0x05: "primers",
    0x06: "notes",
    0x07: "history",
    0x08: "additional sequence properties",
    0x09: "cookie",
    0x0A: "features",
    0x0B: "history node",
    0x0D: "enzyme set",
    0x0E: "custom enzyme sets",
    0x10: "alignable sequence",
    0x11: "alignable sequences",
    0x1C: "enzyme visibilities",
}


@dataclass(frozen=True)
class Packet:
    """One TLV packet."""

    type: int
    data: bytes

    @property
    def name(self) -> str:
        return PACKET_NAMES.get(self.type, f"0x{self.type:02x}")

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self.data)


def iter_packets(source: bytes | BinaryIO) -> Iterator[Packet]:
    """Yield every packet in a SnapGene file.

    ``source`` may be the raw bytes or a binary file object. A truncated trailing
    packet raises :class:`ValueError` rather than being silently dropped, because a
    short read usually means the file is damaged.
    """
    data = source if isinstance(source, (bytes, bytearray)) else source.read()
    offset = 0
    end = len(data)
    while offset < end:
        if end - offset < 5:
            raise ValueError(
                f"truncated packet header at byte {offset}: {end - offset} byte(s) left, need 5"
            )
        packet_type = data[offset]
        (length,) = struct.unpack_from(">I", data, offset + 1)
        start = offset + 5
        stop = start + length
        if stop > end:
            raise ValueError(
                f"packet 0x{packet_type:02x} at byte {offset} claims {length} bytes "
                f"but only {end - start} remain"
            )
        yield Packet(packet_type, bytes(data[start:stop]))
        offset = stop


def pack_packet(packet_type: int, payload: bytes) -> bytes:
    """Serialise one packet."""
    if not 0 <= packet_type <= 0xFF:
        raise ValueError(f"packet type out of range: {packet_type}")
    return struct.pack(">BI", packet_type, len(payload)) + payload


def pack_packets(packets: Iterable[Packet]) -> bytes:
    """Serialise a sequence of packets back into file bytes."""
    return b"".join(pack_packet(p.type, p.data) for p in packets)
