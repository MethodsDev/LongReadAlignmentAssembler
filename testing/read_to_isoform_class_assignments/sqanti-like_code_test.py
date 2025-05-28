#!/usr/bin/env python3

import sys, os, re
import pytest
from unittest.mock import MagicMock, patch

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../pylib"])
)

from MockTranscript import MockTranscript
from Transcript import Transcript
from SQANTI_like_annotator import SQANTI_like_annotator

import pytest
from SQANTI_like_annotator import SQANTI_like_annotator


# --- Pytest fixtures and helpers --- #
@pytest.fixture
def annotator(monkeypatch):
    # Patch the GTF parser to avoid file I/O
    monkeypatch.setattr(
        "SQANTI_like_annotator.GTF_contig_to_transcripts",
        type("Dummy", (), {"parse_GTF_to_Transcripts": staticmethod(lambda x: {})}),
    )
    return SQANTI_like_annotator("dummy.gtf")


def add_ref_transcript(annotator, t):
    annotator.transcript_id_to_obj[t.get_transcript_id()] = t
    chrom_strand = f"{t._chrom}:{t._strand}"
    # Add exons to interval tree
    for exon in t.get_exon_segments():
        annotator.stranded_chrom_exon_itrees[chrom_strand][
            exon[0] : exon[1] + 1
        ] = t.get_transcript_id()
    # Add introns to interval tree and isoform structures
    if t.has_introns():
        intron_str = t.get_introns_string()
        annotator.splice_patterns_to_isoforms[intron_str].add(t.get_transcript_id())
        for intron in t.get_introns():
            annotator.intron_to_isoforms[f"{t._chrom}:{t._strand}:{intron}"].add(
                t.get_transcript_id()
            )
            annotator.stranded_splice_sites.add(f"{t._chrom}:{t._strand}:{intron[0]}")
            annotator.stranded_splice_sites.add(f"{t._chrom}:{t._strand}:{intron[1]}")
            annotator.stranded_chrom_intron_itrees[chrom_strand][
                intron[0] : intron[1] + 1
            ] = t.get_transcript_id()


# --- Now, the actual tests --- #

###########################
## multi-exon feature tests


def test_FSM(annotator):
    # Multi-exon FSM
    ref = MockTranscript("T3", "chr1", "+", [(100, 150), (200, 250), (300, 350)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F3", "chr1", "+", [(100, 150), (200, 250), (300, 350)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F3", feature)
    assert res["sqanti_cat"] == "FSM"
    assert "T3" in res["matching_isoforms"]


def test_ISM(annotator):
    # Multi-exon ISM: feature has a subset of reference introns
    ref = MockTranscript("T4", "chr1", "+", [(100, 150), (200, 250), (300, 350)])
    add_ref_transcript(annotator, ref)
    # ISM: only first two exons, thus only first intron
    feature = MockTranscript("F4", "chr1", "+", [(100, 150), (200, 250)])
    # Patch restrict_splice_matched_isoforms to always return the reference isoform
    annotator.restrict_splice_matched_isoforms = lambda isoforms, obj: isoforms
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F4", feature)
    assert res["sqanti_cat"] == "ISM"


def test_NIC(annotator):
    # NIC: known splice sites, new combination not in reference
    ref1 = MockTranscript("T5a", "chr1", "+", [(100, 150), (200, 450), (600, 700)])
    add_ref_transcript(annotator, ref1)
    # Feature: combines introns from both refs (which isn't a reference pattern)
    feature = MockTranscript(
        "F5", "chr1", "+", [(100, 150), (200, 250), (400, 450), (600, 700)]
    )
    # Add all possible known splice sites
    for pos in [150 + 1, 200 - 1, 250 + 1, 300 - 1, 350 + 1, 400 - 1, 450 + 1, 600 - 1]:
        annotator.stranded_splice_sites.add(f"chr1:+:{pos}")
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F5", feature)
    assert res["sqanti_cat"] == "NIC"


def test_NNIC(annotator):
    # NNIC: novel splice sites
    ref = MockTranscript(
        "T6", "chr1", "+", [(100, 150), (200, 250), (400, 450), (600, 700)]
    )
    add_ref_transcript(annotator, ref)
    # Feature with at least one novel intron
    feature = MockTranscript("F6", "chr1", "+", [(200, 300), (600, 650)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F6", feature)
    assert res["sqanti_cat"] == "NNIC"


def test_genic(annotator):
    # NNIC: novel splice sites
    ref = MockTranscript(
        "T6", "chr1", "+", [(100, 150), (200, 250), (400, 450), (600, 700)]
    )
    add_ref_transcript(annotator, ref)
    # Feature with at least one novel intron
    feature = MockTranscript("F6", "chr1", "+", [(200, 300), (550, 650)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F6", feature)
    assert res["sqanti_cat"] == "genic"


def test_intronic(annotator):
    # NNIC: novel splice sites
    ref = MockTranscript("T6", "chr1", "+", [(100, 150), (600, 700)])
    add_ref_transcript(annotator, ref)
    # Feature with at least one novel intron
    feature = MockTranscript("F6", "chr1", "+", [(200, 300), (400, 500)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F6", feature)
    assert res["sqanti_cat"] == "intronic"


def test_antisense(annotator):
    # NNIC: novel splice sites
    ref = MockTranscript("T6", "chr1", "+", [(100, 150), (600, 700)])
    add_ref_transcript(annotator, ref)
    # Feature with at least one novel intron
    feature = MockTranscript("F6", "chr1", "-", [(100, 200), (400, 650)])
    res = annotator.classify_alignment_or_isoform("chr1", "-", "F6", feature)
    assert res["sqanti_cat"] == "antisense"


def test_intergenic(annotator):
    # NNIC: novel splice sites
    ref = MockTranscript("T6", "chr1", "+", [(100, 150), (300, 400)])
    add_ref_transcript(annotator, ref)

    ref2 = MockTranscript("T6b", "chr1", "+", [(600, 700), (800, 900)])
    add_ref_transcript(annotator, ref2)

    # Feature with at least one novel intron
    feature = MockTranscript("F6", "chr1", "+", [(450, 475), (500, 550)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F6", feature)
    assert res["sqanti_cat"] == "intergenic"


###########################
# single-exon feature tests


def test_se_FM(annotator):
    # Single-exon FSM
    ref = MockTranscript("T1", "chr1", "+", [(100, 200)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F1", "chr1", "+", [(100, 200)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F1", feature)
    assert res["sqanti_cat"] == "se_FM"
    assert "T1" in res["matching_isoforms"]


def test_se_IM(annotator):
    # Single-exon, partial overlap (should be se_IM)
    ref = MockTranscript("T2", "chr1", "+", [(100, 200)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F2", "chr1", "+", [(150, 200)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F2", feature)
    assert res["sqanti_cat"] == "se_IM"
    assert "T2" in res["matching_isoforms"]


def test_se_intergenic(annotator):
    # Single-exon, overlaps an exon but not enough for se_FM/se_IM
    ref1 = MockTranscript("T7", "chr1", "+", [(100, 200)])
    add_ref_transcript(annotator, ref1)
    ref2 = MockTranscript("T7", "chr1", "+", [(500, 600)])
    add_ref_transcript(annotator, ref1)
    feature = MockTranscript("F7", "chr1", "+", [(300, 400)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F7", feature)
    assert res["sqanti_cat"] == "se_intergenic"


def test_se_intronic(annotator):
    # Single-exon, overlaps an intron (but not an exon)
    ref = MockTranscript("T8", "chr1", "+", [(100, 150), (200, 250)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F8", "chr1", "+", [(151, 199)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F8", feature)
    assert res["sqanti_cat"] == "se_intronic"


def test_se_antisense(annotator):
    # Single-exon, overlaps exon but on opposite strand
    ref = MockTranscript("T9", "chr1", "-", [(100, 200)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F9", "chr1", "+", [(100, 200)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F9", feature)
    assert res["sqanti_cat"] == "se_antisense"


def test_se_genic(annotator):
    # Single-exon, does not overlap any exon or intron
    ref = MockTranscript("T10", "chr1", "+", [(100, 200), (350, 500)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F10", "chr1", "+", [(300, 500)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F10", feature)
    assert res["sqanti_cat"] == "se_genic"


def test_se_exonic(annotator):
    # Single-exon, does not overlap any exon or intron
    ref = MockTranscript("T10", "chr1", "+", [(100, 200), (250, 500)])
    add_ref_transcript(annotator, ref)
    feature = MockTranscript("F10", "chr1", "+", [(300, 350)])
    res = annotator.classify_alignment_or_isoform("chr1", "+", "F10", feature)
    assert res["sqanti_cat"] == "se_exonic"
