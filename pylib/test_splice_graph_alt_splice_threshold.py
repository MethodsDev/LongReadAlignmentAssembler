#!/usr/bin/env python3

import networkx as nx

import LRAA_Globals
from GenomeFeature import Exon, Intron
from Splice_graph import Splice_graph


def _initialized_params(min_alt_splice_freq):
    return {
        "read_aln_gap_merge_int": Splice_graph._read_aln_gap_merge_int,
        "inter_exon_segment_merge_dist": Splice_graph._inter_exon_segment_merge_dist,
        "max_genomic_contig_length": Splice_graph._max_genomic_contig_length,
        "min_alt_splice_freq": min_alt_splice_freq,
        "min_alt_unspliced_freq": Splice_graph._min_alt_unspliced_freq,
        "max_intron_length_for_exon_segment_filtering": (
            Splice_graph._max_intron_length_for_exon_segment_filtering
        ),
        "min_intron_support": Splice_graph._min_intron_support,
        "remove_unspliced_introns": Splice_graph._remove_unspliced_introns,
    }


def _graph_with_introns(*introns):
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "+"
    splice_graph._intron_objs = {
        "{}:{}".format(*intron.get_coords()): intron for intron in introns
    }
    splice_graph._splice_graph = nx.DiGraph()
    splice_graph._splice_graph.add_nodes_from(introns)
    return splice_graph


def test_initialized_alt_splice_threshold_controls_all_intron_pruning(monkeypatch):
    original_params = _initialized_params(Splice_graph._min_alt_splice_freq)
    initialized_threshold = 0.30
    monkeypatch.setitem(LRAA_Globals.config, "min_alt_splice_freq", 0.01)
    monkeypatch.setitem(
        LRAA_Globals.config, "aggregate_adjacent_splice_boundaries", False
    )

    try:
        Splice_graph.init_sg_params(_initialized_params(initialized_threshold))
        assert Splice_graph._min_alt_splice_freq == initialized_threshold
        assert (
            Splice_graph._min_alt_splice_freq
            != LRAA_Globals.config["min_alt_splice_freq"]
        )

        dominant_boundary = Intron("chr1", 101, 200, "+", 8)
        alternative_boundary = Intron("chr1", 101, 300, "+", 2)
        boundary_graph = _graph_with_introns(
            dominant_boundary, alternative_boundary
        )
        boundary_graph._prune_spurious_introns_shared_boundary("left")

        assert dominant_boundary in boundary_graph._splice_graph
        assert alternative_boundary not in boundary_graph._splice_graph

        dominant_island = Intron("chr1", 401, 500, "+", 8)
        alternative_island = Intron("chr1", 601, 700, "+", 2)
        island_graph = _graph_with_introns(dominant_island, alternative_island)
        exon = Exon("chr1", 501, 600, "+", 10)
        island_graph._splice_graph.add_edges_from(
            [(dominant_island, exon), (alternative_island, exon)]
        )
        island_graph._prune_low_support_introns()

        assert dominant_island in island_graph._splice_graph
        assert alternative_island not in island_graph._splice_graph
    finally:
        Splice_graph.init_sg_params(original_params)
