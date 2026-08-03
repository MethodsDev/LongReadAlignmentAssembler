#!/usr/bin/env python3

import networkx as nx
import pytest

import LRAA_Globals
from GenomeFeature import Exon, Intron, PolyAsite
from Splice_graph import Splice_graph


def _build_terminal_spur(length, side="right", with_alternative_intron=True):
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "+"

    if side == "right":
        internal_exon = Exon("chr1", 100, 200, "+", 10)
        spur = Exon("chr1", 201, 200 + length, "+", 1)
        intron = Intron("chr1", 201, 300, "+", 10)
        edges = [(internal_exon, spur)]
        if with_alternative_intron:
            edges.append((internal_exon, intron))
    else:
        spur = Exon("chr1", 101 - length, 100, "+", 1)
        internal_exon = Exon("chr1", 101, 200, "+", 10)
        intron = Intron("chr1", 1, 100, "+", 10)
        edges = [(spur, internal_exon)]
        if with_alternative_intron:
            edges.append((intron, internal_exon))

    splice_graph._splice_graph = nx.DiGraph(edges)
    return splice_graph, spur


@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize(
    ("max_spur_length", "spur_length", "should_prune"),
    [
        (14, 14, True),
        (14, 15, False),
        (13, 13, True),
        (13, 14, False),
    ],
)
def test_terminal_spur_pruning_uses_configured_maximum(
    monkeypatch, side, max_spur_length, spur_length, should_prune
):
    monkeypatch.setitem(LRAA_Globals.config, "max_exon_spur_length", max_spur_length)
    splice_graph, spur = _build_terminal_spur(spur_length, side)

    splice_graph._prune_exon_spurs_at_introns()

    assert (spur not in splice_graph._splice_graph) is should_prune


def test_default_ont_terminal_spur_maximum_is_14():
    assert LRAA_Globals.config["max_exon_spur_length"] == 14


def test_short_terminal_exon_without_alternative_intron_is_retained(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "max_exon_spur_length", 14)
    splice_graph, spur = _build_terminal_spur(14, with_alternative_intron=False)

    splice_graph._prune_exon_spurs_at_introns()

    assert spur in splice_graph._splice_graph


def test_input_transcript_exon_overlap_exempts_terminal_spur(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "max_exon_spur_length", 14)
    splice_graph, spur = _build_terminal_spur(14)
    lend, rend = spur.get_coords()
    splice_graph._input_transcript_exon_coords_itree[lend : rend + 1] = "reference"

    splice_graph._prune_exon_spurs_at_introns()

    assert spur in splice_graph._splice_graph


def test_polya_connection_prevents_right_terminal_spur_classification(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "max_exon_spur_length", 14)
    splice_graph, spur = _build_terminal_spur(14)
    polya = PolyAsite("chr1", spur.get_coords()[1], spur.get_coords()[1], "+", 5)
    splice_graph._splice_graph.add_edge(spur, polya)

    splice_graph._prune_exon_spurs_at_introns()

    assert spur in splice_graph._splice_graph
