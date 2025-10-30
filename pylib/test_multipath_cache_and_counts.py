#!/usr/bin/env python3

import os
import pytest

import LRAA_Globals as LG
from ReadNameStore import ReadNameStore
from MpReadIdStore import MpReadIdStore
from Splice_graph import Splice_graph
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathGraphNode import MultiPathGraphNode
from Scored_path import Scored_path


def _build_sg_with_exons(n=4):
    sg = Splice_graph()
    Exon.reset_counter()
    exons = []
    start = 100
    for i in range(n):
        e = Exon("chr1", start, start + 49, "+", 1)
        sg._node_id_to_node[e.get_id()] = e
        exons.append(e)
        start += 100
    return sg, exons


def test_multipath_names_no_cache_and_count(tmp_path):
    # Initialize external stores (SQLite fallback)
    base = tmp_path / "cachetest"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    LG.MP_READ_ID_STORE = MpReadIdStore(str(base / "mpreads"))

    sg, exons = _build_sg_with_exons(2)
    mp = MultiPath(sg, [[exons[0].get_id(), exons[1].get_id()]], read_names=set(), read_count=0)
    mp_id = mp.get_id()

    initial = {"r1", "r2"}
    for n in initial:
        rid = LG.READ_NAME_STORE.get_or_add(n)
        LG.MP_READ_ID_STORE.append(mp_id, rid)

    # First call resolves names via store
    first_stream = mp.get_read_names()
    assert first_stream == initial

    # Append a new external association after cache is set
    rid3 = LG.READ_NAME_STORE.get_or_add("r3")
    LG.MP_READ_ID_STORE.append(mp_id, rid3)

    # Count reflects new total via store
    assert mp.get_read_names_count() == 3
    # Names now reflect latest store contents (no caching retained)
    assert mp.get_read_names() == initial.union({"r3"})

    # Cleanup globals
    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None
    if LG.MP_READ_ID_STORE is not None:
        LG.MP_READ_ID_STORE.close()
        LG.MP_READ_ID_STORE = None


def test_multipath_get_read_names_count_without_store():
    sg, exons = _build_sg_with_exons(2)
    # No stores configured and no in-memory names
    mp = MultiPath(sg, [[exons[0].get_id(), exons[1].get_id()]], read_names=set(), read_count=0)
    assert mp.get_read_names_count() == 0


def test_scored_path_fallback_scoring_by_counts():
    sg, exons = _build_sg_with_exons(4)

    class DummyMPG:
        def __init__(self, sg):
            self._sg = sg
        def get_splice_graph(self):
            return self._sg

    # two disjoint paths with no read names (forces fallback)
    mpA = MultiPath(sg, [[exons[0].get_id(), exons[1].get_id()]], read_names=set(), read_count=0)
    mpB = MultiPath(sg, [[exons[2].get_id(), exons[3].get_id()]], read_names=set(), read_count=0)

    lendA, rendA = mpA.get_coords()
    lendB, rendB = mpB.get_coords()

    nodeA = MultiPathGraphNode(mpA, count=3, lend=lendA, rend=rendA, mpg=DummyMPG(sg))
    nodeB = MultiPathGraphNode(mpB, count=5, lend=lendB, rend=rendB, mpg=DummyMPG(sg))

    # Register containment so both nodes are represented without requiring spacer merging
    nodeA.add_containment(nodeB)

    spath = Scored_path([nodeA])
    # With no names, initial score falls back to summed counts of represented nodes
    assert spath.get_initial_score() == 8

    # Rescoring with exclusions should NOT apply fallback; score remains 0 (IDs absent, exclusion non-empty)
    spath.rescore(exclude_read_ids={999})
    assert spath.get_score() == 0
