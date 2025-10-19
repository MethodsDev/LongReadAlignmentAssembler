#!/usr/bin/env python3

import os
import pytest
from SQANTI_like_annotator import SQANTI_like_annotator
from MockTranscript import MockTranscript


def write_minimal_gtf(tmp_path):
    gtf_path = os.path.join(tmp_path, "test.gtf")
    with open(gtf_path, "w") as fh:
        # Gene
        fh.write("chr1\tsrc\tgene\t50\t400\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\";\n")
        # Transcript t1 (multi-exon): exon1 100-120, exon2 160-200 (intron 121-159)
        fh.write("chr1\tsrc\texon\t100\t120\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n")
        fh.write("chr1\tsrc\texon\t160\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n")
        # Transcript t2 (multi-exon): exon1 150-180, exon2 220-260 (intron 181-219)
        fh.write("chr1\tsrc\texon\t150\t180\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t2\";\n")
        fh.write("chr1\tsrc\texon\t220\t260\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t2\";\n")
    return gtf_path


@pytest.fixture
def annotator(tmp_path):
    gtf = write_minimal_gtf(tmp_path)
    return SQANTI_like_annotator(gtf)


def test_se_exonic_when_overlaps_multiple_exons_same_strand(annotator):
    # Read overlaps exons from both t1 (exon2) and t2 (exon1) but no introns
    read = MockTranscript("read1", "chr1", "+", [(165, 170)])
    result = annotator.classify_alignment_or_isoform("chr1", "+", "read1", read)
    assert result["sqanti_cat"] == "se_exonic"


def test_se_intergenic_when_no_overlap(annotator):
    # Read outside any annotated exon/intron
    read = MockTranscript("read2", "chr1", "+", [(300, 320)])
    result = annotator.classify_alignment_or_isoform("chr1", "+", "read2", read)
    assert result["sqanti_cat"] == "se_intergenic"