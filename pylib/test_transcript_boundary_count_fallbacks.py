#!/usr/bin/env python3

from GenomeFeature import Exon, PolyAsite, TSS
from MultiPath import MultiPath
from Splice_graph import Splice_graph


def _build_transcript():
    sg = Splice_graph()
    sg._contig_acc = "chr1"
    sg._contig_strand = "+"

    exon = Exon("chr1", 110, 200, "+", 7)
    tss = TSS("chr1", 100, 100, "+", 5)
    polya = PolyAsite("chr1", 210, 210, "+", 9)

    for node in (tss, exon, polya):
        sg._node_id_to_node[node.get_id()] = node

    mp = MultiPath(sg, [[tss.get_id(), exon.get_id(), polya.get_id()]])
    transcript = mp.toTranscript()
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")

    return transcript


def test_tss_read_count_falls_back_to_boundary_node_support():
    transcript = _build_transcript()
    transcript._TSS_read_count = None

    assert transcript.has_TSS() is True
    assert transcript.get_TSS_read_count() == 5
    assert ' TSS_read_count "5";' in transcript.to_GTF_format(include_TPM=False)


def test_polya_read_count_falls_back_to_boundary_node_support():
    transcript = _build_transcript()
    transcript._PolyA_read_count = None

    assert transcript.has_PolyA() is True
    assert transcript.get_PolyA_read_count() == 9
    assert ' PolyA_read_count "9";' in transcript.to_GTF_format(include_TPM=False)


def test_refresh_boundary_annotations_from_simple_path_clears_stale_imported_flags():
    transcript = _build_transcript()
    transcript._imported_has_TSS = False
    transcript._imported_has_POLYA = False
    transcript.set_TSS_read_count(99)
    transcript.set_PolyA_read_count(88)

    transcript.refresh_boundary_annotations_from_simple_path()

    assert transcript.has_TSS() is True
    assert transcript.has_PolyA() is True
    assert transcript.get_TSS_read_count() == 99
    assert transcript.get_PolyA_read_count() == 88


def test_source_annotated_boundary_flags_survive_refresh():
    transcript = _build_transcript()
    transcript._source_has_annotated_TSS = True
    transcript._source_has_annotated_POLYA = True
    transcript._imported_has_TSS = True
    transcript._imported_has_POLYA = True

    transcript.refresh_boundary_annotations_from_simple_path()

    # Imported flags are remap-derived and can be cleared/reset.
    assert transcript.has_TSS() is True
    assert transcript.has_PolyA() is True
    # Source provenance should remain intact for merge tracking.
    assert transcript.has_source_annotated_TSS() is True
    assert transcript.has_source_annotated_PolyA() is True
