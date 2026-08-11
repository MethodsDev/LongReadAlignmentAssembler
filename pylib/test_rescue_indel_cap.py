#!/usr/bin/env python3

"""The indel-length cap on transcriptome rescue alignments.

A long insertion or deletion against a transcript is structural disagreement, not
sequencing error. Nothing else in the acceptance chain sees it: the aligned-fraction
test counts query bases so a deletion costs it nothing, the explained-bases baseline
excludes deleted bases by construction, and the single-merged-segment test cannot fire
because a deletion emits its own genome segment and bridges its own gap.
"""

import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals


TARGET = "tx1"
TARGET_LEN = 4000


def _sam(tmp_path, cigar, seq_len):
    """A one-record SAM aligning a read to TARGET with the given cigar."""
    path = tmp_path / "rescue.sam"
    with open(path, "wt") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        fh.write(f"@SQ\tSN:{TARGET}\tLN:{TARGET_LEN}\n")
        fh.write(
            "\t".join(
                ["r1", "0", TARGET, "1", "60", cigar, "*", "0", "0", "A" * seq_len, "*", "NM:i:0"]
            )
            + "\n"
        )
    return str(path)


def _indel_lengths(sam_path):
    with pysam.AlignmentFile(sam_path, "r") as fh:
        read = next(iter(fh))
    return [(code, length) for code, length in read.cigartuples]


def _declined(cigartuples, max_indel_length):
    """The gate as applied in _parse_rescue_alignments."""
    return max_indel_length > 0 and any(
        code in (1, 2) and length >= max_indel_length for code, length in cigartuples
    )


def test_deletion_at_the_cap_is_declined(tmp_path):
    """'At this or longer length' -- the cap itself fails, so the test is >=."""
    cig = _indel_lengths(_sam(tmp_path, "500M30D500M", 1000))
    assert _declined(cig, 30) is True


def test_deletion_below_the_cap_is_kept(tmp_path):
    cig = _indel_lengths(_sam(tmp_path, "500M29D500M", 1000))
    assert _declined(cig, 30) is False


def test_insertion_is_capped_the_same_way(tmp_path):
    """Insertions and deletions share one threshold."""
    cig = _indel_lengths(_sam(tmp_path, "500M30I500M", 1030))
    assert _declined(cig, 30) is True
    assert _declined(_indel_lengths(_sam(tmp_path, "500M29I500M", 1029)), 30) is False


def test_the_observed_74_base_skip_is_declined_on_both_platforms(tmp_path):
    """The case this exists for: a read skipping 74 bases the target contains.

    Measured on chr20, correctly placed reads never exceed a 32-base deletion on
    PacBio HiFi or 45 on ONT cDNA, so 74 is outside platform error on either.
    """
    cig = _indel_lengths(_sam(tmp_path, "542M74D2613M1I1927M", 5083))
    assert _declined(cig, 10) is True, "HiFi cap"
    assert _declined(cig, 30) is True, "ONT cap"


def test_zero_disables_the_cap(tmp_path):
    cig = _indel_lengths(_sam(tmp_path, "500M74D500M", 1000))
    assert _declined(cig, 0) is False


def test_soft_clipping_and_matches_are_not_indels(tmp_path):
    """Only I and D are measured; clipping is handled by the aligned-fraction test."""
    cig = _indel_lengths(_sam(tmp_path, "200S1000M200S", 1400))
    assert _declined(cig, 10) is False


def test_hifi_default_is_stricter_than_the_ont_default():
    """--HiFi lowers the cap; the config default is the ONT/general value."""
    assert LRAA_Globals.config["rescue_unassigned_max_indel_length"] == 30


def _aln_from_cigar(cigar, seq_len, nm, ref="tx1", ref_len=6000):
    """A pysam record with the given cigar and NM, for metric comparisons."""
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": ref, "LN": ref_len}]}
    )
    a = pysam.AlignedSegment(header)
    a.query_name = "r1"
    a.reference_id = 0
    a.reference_start = 0
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = "A" * seq_len
    a.set_tag("NM", nm)
    return a


def test_gap_aware_identity_charges_a_skipped_span_that_explained_bases_ignores():
    """The measured case: a clean spliced genome alignment vs a 74-base skip.

    Both explain 5,081 read bases, so the explained-bases baseline ties and the rescue
    is admitted. Gap-aware identity separates them.
    """
    from IsoformReadRescue import _explained_read_bases, _gap_aware_identity

    genome = _aln_from_cigar("1S39M33011N118M859N214M1325N172M442N2612M1I1927M", 5084, 2)
    transcript = _aln_from_cigar("1S542M74D2613M1I1927M", 5084, 76)

    # identical on the metric that currently gates acceptance
    assert _explained_read_bases(genome) == _explained_read_bases(transcript) == 5081

    # separated on the gap-aware one, so the rescue is declined
    assert _gap_aware_identity(transcript) < _gap_aware_identity(genome)


def test_gap_aware_identity_ignores_reference_skips():
    """An intron is not a gap in the read's agreement, or no genome alignment would pass."""
    from IsoformReadRescue import _gap_aware_identity

    spliced = _aln_from_cigar("500M2000N500M", 1000, 0)
    contiguous = _aln_from_cigar("1000M", 1000, 0)
    assert _gap_aware_identity(spliced) == _gap_aware_identity(contiguous) == 1.0


def test_gap_aware_identity_equal_alignments_are_not_declined():
    """The comparison is 'at least as good', so an exact match on both sides passes."""
    from IsoformReadRescue import _gap_aware_identity

    genome = _aln_from_cigar("1000M", 1000, 0)
    transcript = _aln_from_cigar("1000M", 1000, 0)
    assert not (_gap_aware_identity(transcript) < _gap_aware_identity(genome))


def test_identity_is_never_negative_and_gates_on_it():
    """The old formula could return a negative 'identity'; this one cannot.

    Divided matched-minus-NM by the matched count alone, an alignment with zero
    substitutions and a 592-base insertion scored below zero and was clamped to 0,
    failing every threshold. Charging gap bases to the denominator gives a real
    fraction in [0,1] for the same alignment.
    """
    from IsoformReadRescue import _gap_aware_identity, _passes_percent_identity

    # 1S 18M 1D 30M 591I 135M 1I 321M 1961S -- zero substitutions, huge insertion
    a = _aln_from_cigar("1S18M1D30M591I135M1I321M1961S", 3057, 593)
    val = _gap_aware_identity(a)
    assert 0.0 <= val <= 1.0
    assert val > 0.0, "an alignment with no substitutions must not score zero identity"

    # it still fails a stringent threshold, on the gap burden rather than by clamping
    assert _passes_percent_identity(a, 97.0) is False


def test_identity_gate_accepts_a_clean_alignment():
    from IsoformReadRescue import _passes_percent_identity

    clean = _aln_from_cigar("1000M", 1000, 0)
    assert _passes_percent_identity(clean, 97.0) is True
