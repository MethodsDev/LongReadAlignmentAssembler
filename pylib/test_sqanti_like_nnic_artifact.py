#!/usr/bin/env python3

"""NNIC_artifact category of the SQANTI-like classifier.

With a genome, an NNIC whose novel junctions (introns absent from the reference)
include one without a canonical motif (GT-AG, GC-AG, AT-AC on the feature's strand)
is NNIC_artifact; column noncanonical_novel_junctions lists those junctions as
lend-rend:motif in the transcribed orientation. Reference introns are not tested.
"""

import os
import pytest
import pysam
from SQANTI_like_annotator import SQANTI_like_annotator
from MockTranscript import MockTranscript

# Reference isoforms on both strands, introns 201-299 and 401-499. The genome is all
# A's, so no reference intron has a canonical motif: they must never be tested.
REF = {
    "P": ("+", [(100, 200), (300, 400), (500, 600)]),
    "M": ("-", [(100, 200), (300, 400), (500, 600)]),
}

# Novel second intron 401-459 (exon 2 kept, new acceptor): motif bases 401-402, 458-459
NOVEL = [(100, 200), (300, 400), (460, 600)]


def _write_genome(path, motifs):
    seq = ["A"] * 1000
    for pos, base in motifs.items():
        seq[pos - 1] = base
    with open(path, "w") as fh:
        fh.write(">chr1\n" + "".join(seq) + "\n")
    pysam.faidx(path)
    return path


def _annotator(tmp_path, motifs=None):
    gtf = os.path.join(tmp_path, "ref.gtf")
    with open(gtf, "w") as fh:
        for tid, (strand, exons) in REF.items():
            for lend, rend in exons:
                fh.write(
                    'chr1\tsrc\texon\t{}\t{}\t.\t{}\t.\tgene_id "g{}"; transcript_id "{}";\n'.format(
                        lend, rend, strand, tid, tid
                    )
                )
    genome = (
        None
        if motifs is None
        else _write_genome(os.path.join(tmp_path, "genome.fa"), motifs)
    )
    return SQANTI_like_annotator(gtf, 50, genome)


def classify(annotator, exons, strand="+"):
    q = MockTranscript("q", "chr1", strand, exons)
    return annotator.classify_alignment_or_isoform("chr1", strand, "q", q)


def test_canonical_novel_junction_stays_NNIC(tmp_path):
    a = _annotator(tmp_path, {401: "G", 402: "T", 458: "A", 459: "G"})
    r = classify(a, NOVEL)
    assert r["sqanti_cat"] == "NNIC"
    assert r["noncanonical_novel_junctions"] == ""


@pytest.mark.parametrize("motif", ["GCAG", "ATAC"])
def test_minor_canonical_motifs(tmp_path, motif):
    a = _annotator(tmp_path, dict(zip((401, 402, 458, 459), motif)))
    assert classify(a, NOVEL)["sqanti_cat"] == "NNIC"


def test_noncanonical_novel_junction_is_artifact(tmp_path):
    a = _annotator(tmp_path, {401: "G", 402: "A", 458: "A", 459: "G"})
    r = classify(a, NOVEL)
    assert r["sqanti_cat"] == "NNIC_artifact"
    # the reference intron 201-299 (AA-AA here) is not listed
    assert r["noncanonical_novel_junctions"] == "401-459:GA-AG"


def test_minus_strand_motif_is_read_on_the_transcribed_strand(tmp_path):
    # top strand CT..AC is GT-AG on the minus strand
    a = _annotator(tmp_path, {401: "C", 402: "T", 458: "A", 459: "C"})
    assert classify(a, NOVEL, "-")["sqanti_cat"] == "NNIC"
    # top strand CT..TC: minus-strand GA-AG
    a = _annotator(tmp_path, {401: "C", 402: "T", 458: "T", 459: "C"})
    r = classify(a, NOVEL, "-")
    assert r["sqanti_cat"] == "NNIC_artifact"
    assert r["noncanonical_novel_junctions"] == "401-459:GA-AG"


def test_without_genome_NNIC_is_unchanged(tmp_path):
    r = classify(_annotator(tmp_path), NOVEL)
    assert r["sqanti_cat"] == "NNIC"
    assert r["noncanonical_novel_junctions"] == ""


def test_reference_categories_are_not_tested(tmp_path):
    # FSM and ISM use only reference introns, which are AA-AA in this genome
    a = _annotator(tmp_path, {})
    assert classify(a, REF["P"][1])["sqanti_cat"] == "FSM"
    assert classify(a, [(150, 200), (300, 400)])["sqanti_cat"] == "ISM"
