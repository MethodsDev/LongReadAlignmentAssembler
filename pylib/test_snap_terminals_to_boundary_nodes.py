#!/usr/bin/env python3
"""
A remapped model whose path ends on a TSS/POLYA node must END at that node.

Regression for a model whose path was trimmed onto a POLYA node 17 bp inside its
terminal exon while its coordinates kept the original end on an A-rich read-end
stack: it reported PolyA "True" and InternalPriming "True" at once, because
has_PolyA() reads the path and the internal-priming test reads the coordinates.
"""

import pytest

import LRAA_Globals
from GenomeFeature import PolyAsite, TSS
from LRAA import LRAA
from Splice_graph import Splice_graph
from Transcript import Transcript


def _lraa_with_nodes(strand, nodes, path):
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = strand
    splice_graph._node_id_to_node = {n.get_id(): n for n in nodes}
    lraa = LRAA(splice_graph)
    # The snap is what is under test, not read-to-graph mapping.
    lraa._map_read_to_graph = lambda *args, **kwargs: list(path)
    return lraa


def _transcript(exons, strand):
    t = Transcript("chr1", exons, strand)
    t.set_transcript_id("draft")
    return t


@pytest.fixture(autouse=True)
def _quiet(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "show_progress_assign_transcripts", False)


def test_plus_strand_polyA_inside_last_exon_snaps_3prime_end():
    polyA = PolyAsite("chr1", 243, 243, "+", 5)
    lraa = _lraa_with_nodes("+", [polyA], ["E:1", "I:1", "E:2", polyA.get_id()])
    t = _transcript([[100, 150], [200, 260]], "+")

    mapped = lraa.assign_transcripts_paths_in_graph([t], snap_terminals_for=[t])

    assert mapped == [t]
    assert t.has_PolyA()
    assert t.get_exon_segments() == [[100, 150], [200, 243]]
    assert t.get_coords() == (100, 243)
    assert t.get_cdna_len() == 51 + 44


def test_minus_strand_polyA_left_and_TSS_right_both_snap():
    polyA = PolyAsite("chr1", 110, 110, "-", 5)
    tss = TSS("chr1", 255, 255, "-", 5)
    lraa = _lraa_with_nodes(
        "-", [polyA, tss], [polyA.get_id(), "E:1", "I:1", "E:2", tss.get_id()]
    )
    t = _transcript([[100, 150], [200, 260]], "-")

    lraa.assign_transcripts_paths_in_graph([t], snap_terminals_for=[t])

    assert t.has_PolyA() and t.has_TSS()
    assert t.get_exon_segments() == [[110, 150], [200, 255]]
    assert t.get_coords() == (110, 255)


def test_models_not_named_for_snapping_keep_coordinates():
    # reference / imported models are remapped without being moved
    polyA = PolyAsite("chr1", 243, 243, "+", 5)
    lraa = _lraa_with_nodes("+", [polyA], ["E:1", "I:1", "E:2", polyA.get_id()])
    t = _transcript([[100, 150], [200, 260]], "+")

    lraa.assign_transcripts_paths_in_graph([t])

    assert t.has_PolyA()
    assert t.get_exon_segments() == [[100, 150], [200, 260]]


@pytest.mark.parametrize("polyA_pos", [260, 275, 190])
def test_node_at_or_outside_the_terminal_exon_is_not_extended_to(polyA_pos):
    polyA = PolyAsite("chr1", polyA_pos, polyA_pos, "+", 5)
    lraa = _lraa_with_nodes("+", [polyA], ["E:1", "I:1", "E:2", polyA.get_id()])
    t = _transcript([[100, 150], [200, 260]], "+")

    lraa.assign_transcripts_paths_in_graph([t], snap_terminals_for=[t])

    assert t.get_exon_segments() == [[100, 150], [200, 260]]


def test_path_without_boundary_nodes_is_untouched():
    lraa = _lraa_with_nodes("+", [], ["E:1", "I:1", "E:2"])
    t = _transcript([[100, 150], [200, 260]], "+")

    lraa.assign_transcripts_paths_in_graph([t], snap_terminals_for=[t])

    assert not t.has_PolyA()
    assert t.get_exon_segments() == [[100, 150], [200, 260]]
