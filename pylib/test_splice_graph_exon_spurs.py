#!/usr/bin/env python3

import networkx as nx
import pytest
import pysam

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


def test_public_build_rebuilds_components_after_terminal_spur_pruning(
    tmp_path, monkeypatch
):
    monkeypatch.setitem(LRAA_Globals.config, "infer_TSS", False)
    monkeypatch.setitem(LRAA_Globals.config, "infer_PolyA", False)
    monkeypatch.setitem(LRAA_Globals.config, "max_exon_spur_length", 14)

    bam_path = tmp_path / "terminal_spur.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 500}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_writer:
        for index in range(5):
            alignment = pysam.AlignedSegment()
            alignment.query_name = f"spliced{index}"
            alignment.query_sequence = "A" * 201
            alignment.flag = 0
            alignment.reference_id = 0
            alignment.reference_start = 99
            alignment.mapping_quality = 60
            alignment.cigartuples = [(0, 101), (3, 100), (0, 100)]
            alignment.query_qualities = pysam.qualitystring_to_array("I" * 201)
            bam_writer.write(alignment)

        for index in range(2):
            alignment = pysam.AlignedSegment()
            alignment.query_name = f"unspliced{index}"
            alignment.query_sequence = "A" * 115
            alignment.flag = 0
            alignment.reference_id = 0
            alignment.reference_start = 99
            alignment.mapping_quality = 60
            alignment.cigartuples = [(0, 115)]
            alignment.query_qualities = pysam.qualitystring_to_array("I" * 115)
            bam_writer.write(alignment)

    pysam.index(str(bam_path))
    contig_sequence = list("A" * 500)
    contig_sequence[200:202] = "GT"
    contig_sequence[298:300] = "AG"

    splice_graph = Splice_graph()
    splice_graph.build_splice_graph_for_contig(
        "chr1",
        "+",
        "".join(contig_sequence),
        str(bam_path),
        None,
        None,
        input_transcripts=None,
        quant_mode=False,
    )

    live_nodes = set(splice_graph._splice_graph)
    component_nodes = {
        node for component in splice_graph._components for node in component
    }
    live_node_ids = {node.get_id() for node in live_nodes}

    assert all(node.get_coords() != (201, 214) for node in live_nodes)
    assert component_nodes == live_nodes
    assert set(splice_graph._node_id_to_node) == live_node_ids
    assert set(splice_graph._node_id_to_component) == live_node_ids
    assert len(splice_graph._components) == 1
    assert set(splice_graph._node_id_to_component.values()) == {0}
