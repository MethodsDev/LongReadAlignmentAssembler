#!/usr/bin/env python3

import os
import pytest

import LRAA_Globals as LG
from ReadNameStore import ReadNameStore

from Splice_graph import Splice_graph
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathGraphNode import MultiPathGraphNode
from Scored_path import Scored_path


def build_minimal_sg_with_exons():
    sg = Splice_graph()
    # Reset exon counter to deterministic IDs if needed
    Exon.reset_counter()
    # Create three simple exons in order
    e1 = Exon("chr1", 100, 150, "+", 1)
    e2 = Exon("chr1", 200, 250, "+", 1)
    e3 = Exon("chr1", 300, 350, "+", 1)
    # Register them in the splice graph so MultiPath can resolve IDs
    sg._node_id_to_node[e1.get_id()] = e1
    sg._node_id_to_node[e2.get_id()] = e2
    sg._node_id_to_node[e3.get_id()] = e3
    return sg, e1, e2, e3


def test_scored_path_reads_without_containments(tmp_path):
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    # initialize name store so string names are mapped to IDs
    base = tmp_path / "spnames1"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    
    class DummyMPG:
        def __init__(self, sg):
            self._sg = sg
        def get_splice_graph(self):
            return self._sg

    # path over two exons; embed read names directly to avoid external store dependency
    mp = MultiPath(sg, [[e1.get_id(), e2.get_id()]], read_names={"r1", "r2"}, read_count=2)
    lend, rend = mp.get_coords()
    mpgn = MultiPathGraphNode(mp, count=2, lend=lend, rend=rend, mpg=DummyMPG(sg))

    sp = Scored_path([mpgn])
    names = sp.get_all_represented_reads()

    assert names == {"r1", "r2"}

    # cleanup
    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None


def test_scored_path_reads_with_containments(tmp_path):
    sg, e1, e2, e3 = build_minimal_sg_with_exons()
    base = tmp_path / "spnames2"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))

    class DummyMPG:
        def __init__(self, sg):
            self._sg = sg
        def get_splice_graph(self):
            return self._sg

    # parent path e1-e3, contained path e1-e2
    mp_parent = MultiPath(sg, [[e1.get_id(), e3.get_id()]], read_names={"p1"}, read_count=1)
    mp_child = MultiPath(sg, [[e1.get_id(), e2.get_id()]], read_names={"c1", "c2"}, read_count=2)

    lend_p, rend_p = mp_parent.get_coords()
    lend_c, rend_c = mp_child.get_coords()

    node_parent = MultiPathGraphNode(mp_parent, count=1, lend=lend_p, rend=rend_p, mpg=DummyMPG(sg))
    node_child = MultiPathGraphNode(mp_child, count=2, lend=lend_c, rend=rend_c, mpg=DummyMPG(sg))

    # manually register containment
    node_parent.add_containment(node_child)

    sp = Scored_path([node_parent])
    names = sp.get_all_represented_reads()

    # union of parent and contained names
    assert names == {"p1", "c1", "c2"}

    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None


class _DummyMPG:
    def __init__(self, sg):
        self._sg = sg

    def get_splice_graph(self):
        return self._sg


def _scored_path(sg, exon_ids, read_names, count):
    mp = MultiPath(sg, [exon_ids], read_names=read_names, read_count=count)
    lend, rend = mp.get_coords()
    mpgn = MultiPathGraphNode(mp, count=count, lend=lend, rend=rend, mpg=_DummyMPG(sg))
    return mp, Scored_path([mpgn])


@pytest.fixture
def name_store(tmp_path):
    base = tmp_path / "synscore"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    preserved = set(LG.SYNTHETIC_READ_IDS)
    LG.SYNTHETIC_READ_IDS.clear()
    yield
    LG.SYNTHETIC_READ_IDS.clear()
    LG.SYNTHETIC_READ_IDS.update(preserved)
    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None


def test_path_supported_only_by_a_synthetic_read_scores_zero(name_store):
    """A structure no read supports must not be selectable as a candidate."""
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    mp, sp = _scored_path(sg, [e1.get_id(), e2.get_id()], {"reftranscript:t1"}, 1)
    assert sp.compute_path_score() == 1  # before registration, it counts

    LG.SYNTHETIC_READ_IDS.update(mp.get_read_ids())
    assert sp.compute_path_score() == 0
    # must not fall through to the node-count fallback, which still holds the read
    assert sp._all_represented_read_ids


def test_synthetic_read_does_not_inflate_a_real_supported_path(name_store):
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    mp, sp = _scored_path(
        sg, [e1.get_id(), e2.get_id()], {"reftranscript:t1", "r1", "r2"}, 3
    )
    synthetic = {
        rid
        for rid in mp.get_read_ids()
        if str(LG.READ_NAME_STORE.get_name(int(rid))).startswith("reftranscript:")
    }
    assert len(synthetic) == 1
    LG.SYNTHETIC_READ_IDS.update(synthetic)
    # scores the two real reads only
    assert sp.compute_path_score() == 2


def test_fallback_still_applies_when_no_read_ids_exist(name_store):
    """The node-count fallback covers missing id data, not zero support."""
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    _, sp = _scored_path(sg, [e1.get_id(), e2.get_id()], set(), 5)
    assert not sp._all_represented_read_ids
    assert sp.compute_path_score() == 5
