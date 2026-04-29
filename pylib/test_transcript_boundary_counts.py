#!/usr/bin/env python3

from GenomeFeature import Exon, PolyAsite, TSS
from MultiPath import MultiPath
from Splice_graph import Splice_graph
from Transcript import GTF_contig_to_transcripts


def test_multipath_transcript_exports_boundary_read_counts():
    sg = Splice_graph()
    sg._contig_acc = 'chr1'
    sg._contig_strand = '+'

    exon = Exon('chr1', 110, 200, '+', 7)
    tss = TSS('chr1', 100, 100, '+', 5)
    polya = PolyAsite('chr1', 210, 210, '+', 9)

    for node in (tss, exon, polya):
        sg._node_id_to_node[node.get_id()] = node

    mp = MultiPath(sg, [[tss.get_id(), exon.get_id(), polya.get_id()]])
    transcript = mp.toTranscript()
    transcript.set_gene_id('g1')
    transcript.set_transcript_id('t1')

    gtf = transcript.to_GTF_format(include_TPM=False)

    assert ' TSS "True";' in gtf
    assert ' PolyA "True";' in gtf
    assert ' TSS_read_count "5";' in gtf
    assert ' PolyA_read_count "9";' in gtf


def test_gtf_parser_round_trips_boundary_read_counts(tmp_path):
    gtf_path = tmp_path / 'boundary_counts.gtf'
    gtf_path.write_text(
        'chr1\tLRAA\ttranscript\t100\t210\t.\t+\t.\t'
        'gene_id "g1"; transcript_id "t1"; '
        'TSS "True"; PolyA "True"; TSS_read_count "5"; PolyA_read_count "9";\n'
        'chr1\tLRAA\texon\t110\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )

    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(gtf_path))
    transcript = contig_to_transcripts['chr1'][0]

    assert transcript.has_TSS() is True
    assert transcript.has_PolyA() is True
    assert transcript.has_source_annotated_TSS() is True
    assert transcript.has_source_annotated_PolyA() is True
    assert transcript.get_TSS_read_count() == 5
    assert transcript.get_PolyA_read_count() == 9
