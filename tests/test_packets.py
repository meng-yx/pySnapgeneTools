"""The TLV packet layer."""

from __future__ import annotations

import struct

import pytest

from pysnapgene import Packet, iter_packets, pack_packet
from pysnapgene.writer.packets import pack_packets


def test_pack_then_iterate():
    raw = pack_packet(0x00, b"hello") + pack_packet(0x0A, b"<Features/>")
    packets = list(iter_packets(raw))
    assert [(p.type, p.data) for p in packets] == [
        (0x00, b"hello"),
        (0x0A, b"<Features/>"),
    ]


def test_packet_header_is_type_then_big_endian_length():
    raw = pack_packet(0x05, b"abcd")
    assert raw[0] == 0x05
    assert struct.unpack_from(">I", raw, 1)[0] == 4
    assert raw[5:] == b"abcd"


def test_empty_payload_is_legal():
    assert list(iter_packets(pack_packet(0x06, b""))) == [Packet(0x06, b"")]


def test_round_trip_many_packets():
    packets = [Packet(i, bytes([i]) * i) for i in range(1, 20)]
    assert list(iter_packets(pack_packets(packets))) == packets


def test_truncated_header_raises():
    with pytest.raises(ValueError, match="truncated packet header"):
        list(iter_packets(b"\x00\x00\x00"))


def test_truncated_payload_raises():
    raw = pack_packet(0x00, b"12345")[:-2]
    with pytest.raises(ValueError, match="claims 5 bytes"):
        list(iter_packets(raw))


def test_packet_type_out_of_range():
    with pytest.raises(ValueError):
        pack_packet(256, b"")


def test_packet_name_is_readable():
    assert Packet(0x0A, b"").name == "features"
    assert Packet(0xF0, b"").name == "0xf0"


def test_iterate_from_file_object(tmp_path):
    path = tmp_path / "x.bin"
    path.write_bytes(pack_packet(0x09, b"SnapGene"))
    with path.open("rb") as handle:
        assert [p.type for p in iter_packets(handle)] == [0x09]
