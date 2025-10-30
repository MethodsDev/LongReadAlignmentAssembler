#!/usr/bin/env python3

import pytest

from Splice_graph import Splice_graph
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathGraphNode import MultiPathGraphNode
from Scored_path import Scored_path
from LRAA import LRAA
import LRAA_Globals


def build_minimal_sg_with_exons():
    sg = Splice_graph()
    # deterministic exon IDs
    Exon.reset_counter()
    e1 = Exon("chr1", 100, 150, "+", 1)
    e2 = Exon("chr1", 200, 250, "+", 1)
    e3 = Exon("chr1", 300, 350, "+", 1)
    e4 = Exon("chr1", 400, 450, "+", 1)
    for e in (e1, e2, e3, e4):
        sg._node_id_to_node[e.get_id()] = e
    return sg, e1, e2, e3, e4


class DummyMPG:
    def __init__(self, sg):
        self._sg = sg

    def get_splice_graph(self):
        return self._sg


def test_validate_pairwise_prunes_contained_when_collapse_enabled():
    sg, e1, e2, e3, e4 = build_minimal_sg_with_exons()

    # longer path contains the shorter path
    mp_long = MultiPath(sg, [[e1.get_id(), e2.get_id(), e3.get_id(), e4.get_id()]], read_names={"rL"}, read_count=1)
    mp_short = MultiPath(sg, [[e1.get_id(), e2.get_id()]], read_names={"rS1", "rS2"}, read_count=2)

    lendL, rendL = mp_long.get_coords()
    lendS, rendS = mp_short.get_coords()

    node_long = MultiPathGraphNode(mp_long, count=1, lend=lendL, rend=rendL, mpg=DummyMPG(sg))
    node_short = MultiPathGraphNode(mp_short, count=2, lend=lendS, rend=rendS, mpg=DummyMPG(sg))

    sp_long = Scored_path([node_long])
    sp_short = Scored_path([node_short])

    # ensure collapse mode enabled
    LRAA_Globals.config["restrict_asm_to_collapse"] = True

    lraa = LRAA(splice_graph=sg, num_parallel_processes=1)

    paths = [sp_long, sp_short]
    lraa._validate_pairwise_incompatibilities(paths)

    # shorter contained path should be pruned
    assert len(paths) == 1
    assert paths[0].get_multiPath_obj().get_simple_path() == sp_long.get_multiPath_obj().get_simple_path()
