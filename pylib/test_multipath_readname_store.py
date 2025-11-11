#!/usr/bin/env python3

import os
import pytest

import LRAA_Globals as LG
from ReadNameStore import ReadNameStore
from MpReadIdStore import MpReadIdStore
from Splice_graph import Splice_graph
from GenomeFeature import Exon
from MultiPath import MultiPath


def build_minimal_sg():
    sg = Splice_graph()
    Exon.reset_counter()
    e1 = Exon("chr1", 100, 150, "+", 1)
    e2 = Exon("chr1", 200, 250, "+", 1)
    sg._node_id_to_node[e1.get_id()] = e1
    sg._node_id_to_node[e2.get_id()] = e2
    return sg, e1, e2


def test_multipath_external_store_read_names_and_count(tmp_path):
    # Initialize external stores (default backend is in-memory)
    base = tmp_path / "readstore"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    LG.MP_READ_ID_STORE = MpReadIdStore(str(base / "mpreads"))

    sg, e1, e2 = build_minimal_sg()

    # Create MultiPath without in-memory names to force external streaming
    mp = MultiPath(sg, [[e1.get_id(), e2.get_id()]], read_names=set(), read_count=0)
    mp_id = mp.get_id()

    # Append associations externally
    names = {"a", "b", "c"}
    for n in names:
        rid = LG.READ_NAME_STORE.get_or_add(n)
        LG.MP_READ_ID_STORE.append(mp_id, rid)

    # Validate count comes from external store
    assert mp.get_read_names_count() == len(names)

    # Validate names stream from store and match set
    streamed = mp.get_read_names()
    assert streamed == names

    # Cleanup globals
    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None
    if LG.MP_READ_ID_STORE is not None:
        LG.MP_READ_ID_STORE.close()
        LG.MP_READ_ID_STORE = None


def test_memory_backend_roundtrip(monkeypatch):
    monkeypatch.setenv("LRAA_READSTORE_BACKEND", "memory")
    try:
        name_store = ReadNameStore("unused")
        rid_a = name_store.get_or_add("readA")
        rid_b = name_store.get_or_add("readB")
        assert rid_a != rid_b
        assert name_store.get_or_add("readA") == rid_a
        assert name_store.get_id("readB") == rid_b
        assert name_store.get_name(rid_a) == "readA"

        mp_store = MpReadIdStore("unused")
        mp_store.append("mp1", rid_a)
        mp_store.append("mp1", rid_b)
        mp_store.append("mp1", rid_a)

        assert mp_store.count("mp1") == 3
        streamed = list(mp_store.iter_read_ids("mp1"))
        assert streamed.count(rid_a) == 2
        assert streamed.count(rid_b) == 1

        name_store.close()
        mp_store.close()
    finally:
        monkeypatch.delenv("LRAA_READSTORE_BACKEND", raising=False)
