#!/usr/bin/env python3
"""get_exon_predecessor must return the predecessor node, not the Exon class.

The exon-island walk in _prune_low_support_introns only ever tests the result
against None, so a bug returning the class was invisible: islands were still
found correctly. It made the function unusable for walking backwards, which is
a trap for anyone who reaches for the obvious counterpart to
get_exon_successor. Pinned as a contract, since the two must be symmetric.
"""
import os, sys

sys.path.insert(0, os.path.dirname(__file__))

import networkx as nx

import Splice_graph
from GenomeFeature import Exon


def _chain():
    """Three abutting exon segments wired a -> b -> c, as the graph builds them."""
    sg = Splice_graph.Splice_graph()
    sg._splice_graph = nx.DiGraph()
    a = Exon("chr1", 100, 199, "+", 10)
    b = Exon("chr1", 200, 299, "+", 10)
    c = Exon("chr1", 300, 399, "+", 10)
    sg._splice_graph.add_edge(a, b)
    sg._splice_graph.add_edge(b, c)
    return sg, a, b, c


def test_predecessor_returns_the_node():
    sg, a, b, c = _chain()
    assert sg.get_exon_predecessor(b) is a
    assert sg.get_exon_predecessor(c) is b


def test_predecessor_is_none_at_the_head():
    sg, a, _, _ = _chain()
    assert sg.get_exon_predecessor(a) is None


def test_walking_backwards_reaches_the_head():
    """The property the old code silently made impossible."""
    sg, a, _, c = _chain()
    node, seen = c, [c]
    while (node := sg.get_exon_predecessor(node)) is not None:
        seen.append(node)
    assert seen[-1] is a and len(seen) == 3


def test_the_two_directions_are_inverse():
    sg, a, b, c = _chain()
    for n in (a, b):
        assert sg.get_exon_predecessor(sg.get_exon_successor(n)) is n
