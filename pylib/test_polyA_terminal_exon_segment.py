#!/usr/bin/env python3

"""An untrimmed polyA tail can be ALIGNED rather than soft-clipped.

In splice mode the aligner will place such a tail on a genomic A-run kilobases
downstream, joined by a spurious intron. By then the tail is an exon, so the
soft-clip polyA handling in Pretty_alignment never examines it, and the
alignment contributes both a 3' end in the wrong place and a junction no
transcript has.
"""

import sys
from pathlib import Path

import pysam

REPO_ROOT = Path(__file__).resolve().parents[1]
for _import_dir in (str(REPO_ROOT / "pylib"), str(REPO_ROOT / "util")):
    if _import_dir not in sys.path:
        sys.path.insert(0, _import_dir)

import LRAA_Globals
import Util_funcs


M, I, D, N, S = 0, 1, 2, 3, 4


def _aln(cigar, sequence, flag=0, start=1000):
    aln = pysam.AlignedSegment()
    aln.query_name = "read"
    aln.flag = flag
    aln.reference_id = 0
    aln.reference_start = start
    aln.mapping_quality = 60
    aln.cigar = cigar
    aln.query_sequence = sequence
    aln.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    aln.set_tag("NM", 0)
    consumed = sum(length for code, length in cigar if code in (M, I, S, 7, 8))
    assert consumed == len(sequence), "cigar and sequence disagree on query length"
    return aln


def test_aligned_polyA_tail_behind_a_junction_is_detected():
    aln = _aln([(M, 900), (N, 7412), (M, 15)], "C" * 900 + "A" * 15)
    assert Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_tail_separated_from_the_junction_by_an_insertion_is_detected():
    """The shape real alignments take: `...24M 7413N 1I 15M 37S`.

    A short indel does not end an exon, so the terminal exon segment reaches
    across it. Requiring N immediately inside the last aligned block missed this
    entirely, which is the common case rather than a corner one.
    """
    cigar = [(M, 900), (N, 7413), (I, 1), (M, 15), (S, 37)]
    aln = _aln(cigar, "C" * 900 + "A" + "A" * 15 + "G" * 37)
    assert Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_polyT_at_the_left_of_a_reverse_alignment_is_detected():
    """SAM stores the sequence genome-forward, so a polyA tail on a reverse
    alignment reads as polyT at the left end -- the base convention already used
    for soft-clipped tails in Pretty_alignment._set_read_soft_clipping_info."""
    aln = _aln([(M, 15), (N, 7412), (M, 900)], "T" * 15 + "C" * 900, flag=16)
    assert Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_terminal_segment_with_no_junction_behind_it_is_left_alone():
    """Where the transcript simply ends there is no spurious junction, and a
    soft-clipped tail there is already handled downstream."""
    aln = _aln([(M, 400), (D, 2), (M, 500), (S, 20)], "C" * 900 + "A" * 20)
    assert not Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_a_long_A_rich_terminal_exon_is_not_a_tail():
    aln = _aln([(M, 400), (N, 1000), (M, 40)], "C" * 400 + "A" * 40)
    assert not Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_a_short_terminal_segment_of_mixed_bases_is_not_a_tail():
    aln = _aln([(M, 400), (N, 1000), (M, 15)], "C" * 400 + "ACGTACGTACGTACG")
    assert not Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_a_tail_shorter_than_min_PolyA_ident_length_is_not_a_tail():
    aln = _aln([(M, 400), (N, 1000), (M, 4)], "C" * 400 + "AAAA")
    assert not Util_funcs.has_polyA_terminal_exon_segment(aln)


def test_quantification_discards_the_alignment():
    aln = _aln([(M, 900), (N, 7412), (M, 15)], "C" * 900 + "A" * 15)
    assert Util_funcs.quant_discard_reason(aln) == "polyA_terminal_segment"


def test_the_exclusion_can_be_switched_off():
    aln = _aln([(M, 900), (N, 7412), (M, 15)], "C" * 900 + "A" * 15)
    prior = LRAA_Globals.config["no_exclude_polyA_terminal_segment"]
    LRAA_Globals.config["no_exclude_polyA_terminal_segment"] = True
    try:
        assert Util_funcs.quant_discard_reason(aln) is None
    finally:
        LRAA_Globals.config["no_exclude_polyA_terminal_segment"] = prior
