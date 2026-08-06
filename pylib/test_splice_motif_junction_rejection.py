#!/usr/bin/env python3

"""A read junction enters the graph only when its splice motif implies the same strand
the graph is being built for. Junctions declined for the two possible reasons are
recorded, because a read whose junctions are all declined has no path through the graph
and its locus then cannot be reconstructed from reads alone.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from Splice_graph import Splice_graph


def _contig_with_junctions(motifs, exon_len=10, intron_len=20):
    """Build a sequence whose successive gaps carry the requested motifs.

    Returns (sequence, alignment_segments) with 1-based inclusive segment coords.
    """
    seq = []
    segments = []
    pos = 1
    for i, motif in enumerate(motifs + [None]):
        segments.append((pos, pos + exon_len - 1))
        seq.append("A" * exon_len)
        pos += exon_len
        if motif is None:
            break
        left, right = motif[:2], motif[2:]
        seq.append(left + "T" * (intron_len - 4) + right)
        pos += intron_len
    return "".join(seq), segments


def _classify(strand, motifs):
    sequence, segments = _contig_with_junctions(motifs)
    sg = Splice_graph()
    sg._contig_acc = "chrT"
    sg._contig_strand = strand
    sg._contig_seq_str = sequence
    admitted = sg._get_introns_matching_splicing_consensus(segments)
    return (
        admitted,
        sg._junctions_rejected_wrong_strand,
        sg._junctions_rejected_no_motif,
    )


def test_plus_strand_admits_only_top_strand_motifs():
    # GTAG implies '+', CTAC implies '-', AAAA implies neither
    admitted, wrong_strand, no_motif = _classify("+", ["GTAG", "CTAC", "AAAA"])

    assert [(lend, rend, "+") for lend, rend, _ in admitted] == admitted
    assert len(admitted) == 1
    assert len(wrong_strand) == 1
    assert len(no_motif) == 1


def test_minus_strand_admits_only_bottom_strand_motifs():
    admitted, wrong_strand, no_motif = _classify("-", ["GTAG", "CTAC", "AAAA"])

    assert len(admitted) == 1
    assert admitted[0][2] == "-"
    # the GTAG junction is canonical, but in the opposite orientation
    assert len(wrong_strand) == 1
    assert len(no_motif) == 1


def test_junctions_without_a_motif_are_recorded_separately():
    """The no-motif case previously produced no record at all."""
    admitted, wrong_strand, no_motif = _classify("+", ["AAAA", "AAAA"])

    assert admitted == []
    assert wrong_strand == set()
    assert len(no_motif) == 2


def test_all_canonical_junctions_leave_no_rejections():
    admitted, wrong_strand, no_motif = _classify("+", ["GTAG", "GCAG", "ATAC"])

    assert len(admitted) == 3
    assert wrong_strand == set()
    assert no_motif == set()
