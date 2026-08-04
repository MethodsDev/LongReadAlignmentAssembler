#!/usr/bin/env python3

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import networkx as nx

import LRAA_Globals
from GenomeFeature import Exon, Intron, PolyAsite, TSS
from IsoformReadRescue import (
    _build_transcript_models,
    _merge_contiguous_genomic_segments,
    _normalize_read_identifier,
    _parse_rescue_alignments,
)
from LRAA import LRAA
from Splice_graph import Splice_graph
from Transcript import Transcript


def test_merge_contiguous_genomic_segments_collapses_split_exons():
    segments = [
        (100, 110),
        (111, 125),
        (200, 210),
        (211, 211),
        (300, 350),
    ]

    assert _merge_contiguous_genomic_segments(segments) == [
        (100, 125),
        (200, 211),
        (300, 350),
    ]


def test_merge_contiguous_genomic_segments_preserves_true_introns():
    segments = [
        (100, 150),
        (250, 300),
    ]

    assert _merge_contiguous_genomic_segments(segments) == segments


def _build_minus_strand_boundary_fixture():
    Exon.reset_counter()
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "-"

    polya = PolyAsite("chr1", 100, 100, "-", 10)
    exon_low = Exon("chr1", 100, 149, "-", 10)
    intron = Intron("chr1", 150, 199, "-", 10)
    exon_high = Exon("chr1", 200, 249, "-", 10)
    tss = TSS("chr1", 249, 249, "-", 10)

    splice_graph._splice_graph = nx.DiGraph(
        [
            (polya, exon_low),
            (exon_low, intron),
            (intron, exon_high),
            (exon_high, tss),
        ]
    )
    splice_graph._intron_objs = {"150:199": intron}
    splice_graph._finalize_splice_graph()

    transcript = Transcript("chr1", [[100, 149], [200, 249]], "-")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")
    transcript.set_simple_path(
        [
            polya.get_id(),
            exon_low.get_id(),
            intron.get_id(),
            exon_high.get_id(),
            tss.get_id(),
        ]
    )
    transcript_models = _build_transcript_models(splice_graph, [transcript], "A" * 300)

    return splice_graph, transcript_models, polya, exon_low, intron, exon_high, tss


def test_minus_strand_rescue_soft_clips_follow_transcript_ends(tmp_path, monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "max_soft_clip_at_TSS", 0)
    monkeypatch.setitem(LRAA_Globals.config, "max_soft_clip_at_PolyA", 0)
    (
        splice_graph,
        transcript_models,
        polya,
        exon_low,
        intron,
        exon_high,
        tss,
    ) = _build_minus_strand_boundary_fixture()

    rescue_sam = tmp_path / "minus_rescue.sam"
    sequence = "C" * 100
    rescue_sam.write_text(
        "@HD\tVN:1.6\tSO:unknown\n"
        "@SQ\tSN:g1^t1\tLN:100\n"
        f"clip_tx_3prime\t0\tg1^t1\t1\t60\t95M5S\t*\t0\t0\t{sequence}\t*\tAS:i:95\tNM:i:0\n"
        f"clip_tx_5prime\t0\tg1^t1\t6\t60\t5S95M\t*\t0\t0\t{sequence}\t*\tAS:i:95\tNM:i:0\n"
    )

    _, read_details = _parse_rescue_alignments(
        str(rescue_sam),
        splice_graph,
        transcript_models,
        read_path_mapper=LRAA(splice_graph)._map_read_to_graph,
    )
    paths_by_read = {
        read_name: read_details["read_name_to_multipaths"][
            _normalize_read_identifier(read_name)
        ][0].get_simple_path()
        for read_name in ("clip_tx_3prime", "clip_tx_5prime")
    }

    assert paths_by_read["clip_tx_3prime"] == [
        exon_low.get_id(),
        intron.get_id(),
        exon_high.get_id(),
        tss.get_id(),
    ]
    assert paths_by_read["clip_tx_5prime"] == [
        polya.get_id(),
        exon_low.get_id(),
        intron.get_id(),
        exon_high.get_id(),
    ]
