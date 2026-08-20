#!/usr/bin/env python3

"""The single retention policy, and the thresholds it must be told rather than read.

Util_funcs.quant_discard_reason is the one place that decides whether
quantification keeps an alignment.  Bam_alignment_extractor consumes it; cut
selection consumes it for emission.  The hazard these tests guard is a consumer
reading a config default where the run in front of it uses something else --
--HiFi raises min_per_id to 97, and LRAA --quant_only swaps
min_mapping_quality_for_final_quant into min_mapping_quality before filtering.
Both of those defaults coincide with the plain ones, so a wrong read is silent.
"""

import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
import Util_funcs


def _aln(nm=0, mapq=60, length=1000, flag=0, tags=True):
    aln = pysam.AlignedSegment()
    aln.query_name = "r"
    aln.flag = flag
    aln.reference_id = 0
    aln.reference_start = 100
    aln.mapping_quality = mapq
    aln.cigar = [(0, length)]
    aln.query_sequence = "A" * length
    aln.query_qualities = pysam.qualitystring_to_array("I" * length)
    if tags:
        aln.set_tag("NM", nm)
    return aln


def test_thresholds_are_taken_from_the_caller_not_the_config():
    """The whole point of the arguments: a caller may know better than the config.

    Fails if the predicate reads LRAA_Globals for a value it was handed, which is
    the bug shape both --HiFi and the final-quant MAPQ swap produce.
    """

    read = _aln(nm=100)  # 90% identity over 1000 aligned bases

    assert Util_funcs.quant_discard_reason(read, min_per_id=80.0) is None
    assert Util_funcs.quant_discard_reason(read, min_per_id=97.0) == "low_perID"

    strong = _aln(mapq=5)
    assert Util_funcs.quant_discard_reason(strong, min_mapping_quality=0) is None
    assert (
        Util_funcs.quant_discard_reason(strong, min_mapping_quality=30)
        == "min_mapping_quality"
    )


def test_the_final_quant_mapping_quality_is_a_distinct_setting():
    """Both default to 0, so a consumer reading the wrong one fails silently.

    Pinned with the two made to differ, which is the only configuration where the
    mistake is observable.
    """

    saved = (
        LRAA_Globals.config["min_mapping_quality"],
        LRAA_Globals.config["min_mapping_quality_for_final_quant"],
    )
    try:
        LRAA_Globals.config["min_mapping_quality"] = 0
        LRAA_Globals.config["min_mapping_quality_for_final_quant"] = 30
        read = _aln(mapq=5)

        # discovery keeps it; final quant does not.  A consumer that must match the
        # quant step has to ask for the final-quant value explicitly.
        assert (
            Util_funcs.quant_discard_reason(
                read,
                min_mapping_quality=LRAA_Globals.config["min_mapping_quality"],
            )
            is None
        )
        assert (
            Util_funcs.quant_discard_reason(
                read,
                min_mapping_quality=LRAA_Globals.config[
                    "min_mapping_quality_for_final_quant"
                ],
            )
            == "min_mapping_quality"
        )
    finally:
        LRAA_Globals.config["min_mapping_quality"] = saved[0]
        LRAA_Globals.config["min_mapping_quality_for_final_quant"] = saved[1]


def test_a_record_with_no_reference_is_rejected():
    """The strand split rejects reference_id < 0; this policy must not be looser."""

    read = _aln()
    read.reference_id = -1
    assert Util_funcs.quant_discard_reason(read) in ("unmapped", "no_chromosome")


@pytest.mark.parametrize(
    "flag,reason",
    [(256, "secondary"), (2048, "supplementary"), (1024, "duplicate"), (512, "qcfail")],
)
def test_structural_rejections_report_their_reason(flag, reason):
    """Reasons are the extractor's discard-counter keys, so they must stay stable."""

    assert Util_funcs.quant_discard_reason(_aln(flag=flag)) == reason


def test_an_alignment_without_nm_is_not_judged_on_identity():
    """No tag, nothing to measure: the extractor keeps it, so this must too."""

    assert Util_funcs.quant_discard_reason(_aln(tags=False), min_per_id=99.0) is None


def _write_bam(path, records):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for name, flag, mapq, nm in records:
            aln = pysam.AlignedSegment(fh.header)
            aln.query_name = name
            aln.flag = flag
            aln.reference_id = 0
            aln.reference_start = 1000
            aln.mapping_quality = mapq
            aln.cigar = [(0, 1000)]
            aln.query_sequence = "A" * 1000
            aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
            aln.set_tag("NM", nm)
            fh.write(aln)
    pysam.index(str(path))
    return path


def test_records_parsed_from_a_bam_are_judged_the_same_as_hand_built_ones(tmp_path):
    """Every branch, against records that went through pysam's BAM parser.

    The cases above construct AlignedSegment objects directly, which is the cheap way to
    reach a branch but not the way records actually arrive. A parsed record differs in the
    places this predicate reads: NM comes back with whatever type serialization gave it,
    get_cigar_stats() runs over parsed CIGAR ops, and reference_id and the pairing flags are
    resolved against a real header. A policy that judged hand-built segments correctly and
    parsed ones differently would be wrong exactly where it is used.

    Worth having beyond the round-trip: the real corpora carry zero secondary records, so
    these flags are the only way that branch is reached at all.
    """
    # name, flag, mapq, NM
    records = [
        ("keeper", 0, 60, 0),
        ("secondary", 256, 60, 0),
        ("supplementary", 2048, 60, 0),
        ("duplicate", 1024, 60, 0),
        ("qcfail", 512, 60, 0),
        # paired (1) and mapped, but not flagged proper-pair (2)
        ("improper", 1, 60, 0),
        ("lowmapq", 0, 3, 0),
        # 200 mismatches over 1000 aligned bases -> 80% identity
        ("lowid", 0, 60, 200),
    ]
    bam = _write_bam(tmp_path / "cases.bam", records)

    expected = {
        "keeper": None,
        "secondary": "secondary",
        "supplementary": "supplementary",
        "duplicate": "duplicate",
        "qcfail": "qcfail",
        "improper": "improper_pair",
        "lowmapq": "min_mapping_quality",
        "lowid": "low_perID",
    }

    seen = {}
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        for read in fh.fetch(until_eof=True):
            seen[read.query_name] = Util_funcs.quant_discard_reason(
                read,
                max_intron_length=200000,
                min_mapping_quality=10,
                min_per_id=90.0,
            )

    assert seen == expected


def test_rdna_masked_reads_are_rejected_when_their_blocks_overlap_the_mask(tmp_path):
    """The mask is a real IntervalTree, exercised through a parsed record.

    Positioned deliberately after the structural/mapq checks and before the
    intron/per-id ones: a read failing only the mask (mapq 60, no mismatches, no
    long intron) must be rejected as "rdna_masked" specifically, not fall through
    to None or to some other reason.
    """
    import intervaltree

    records = [
        ("in_mask", 0, 60, 0),  # spans [1000, 2000)
        ("outside_mask", 0, 60, 0),
    ]
    bam = tmp_path / "cases.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header) as fh:
        for i, (name, flag, mapq, nm) in enumerate(records):
            aln = pysam.AlignedSegment(fh.header)
            aln.query_name = name
            aln.flag = flag
            aln.reference_id = 0
            aln.reference_start = 1000 + i * 5000
            aln.mapping_quality = mapq
            aln.cigar = [(0, 1000)]
            aln.query_sequence = "A" * 1000
            aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
            aln.set_tag("NM", nm)
            fh.write(aln)
    pysam.index(str(bam))

    mask = {"chr1": intervaltree.IntervalTree()}
    mask["chr1"][1000:2000] = True  # exactly the "in_mask" read's span

    seen = {}
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        for read in fh.fetch(until_eof=True):
            seen[read.query_name] = Util_funcs.quant_discard_reason(
                read, rdna_mask=mask
            )
    assert seen == {"in_mask": "rdna_masked", "outside_mask": None}


def test_rdna_mask_is_read_from_config_when_not_passed_explicitly(tmp_path):
    """Same "None reads the config default" convention as every other threshold.

    normalize_bam_by_strand.py runs as a separate process and threads its mask
    through explicitly (see _record_is_evidence); every in-process caller --
    Bam_alignment_extractor, Pretty_alignment_manager, StreamingQuant -- relies
    on this fallback instead, so it is load-bearing for the masking feature as a
    whole, not just a convenience.
    """
    import intervaltree

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}],
    }
    bam = tmp_path / "one.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as fh:
        aln = pysam.AlignedSegment(fh.header)
        aln.query_name = "r"
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = 1000
        aln.mapping_quality = 60
        aln.cigar = [(0, 1000)]
        aln.query_sequence = "A" * 1000
        aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
        aln.set_tag("NM", 0)
        fh.write(aln)
    pysam.index(str(bam))

    mask = {"chr1": intervaltree.IntervalTree()}
    mask["chr1"][1000:2000] = True

    saved = LRAA_Globals.config.get("rdna_mask_intervals")
    try:
        with pysam.AlignmentFile(str(bam), "rb") as fh:
            read = next(fh.fetch(until_eof=True))

        LRAA_Globals.config["rdna_mask_intervals"] = None
        assert Util_funcs.quant_discard_reason(read) is None

        LRAA_Globals.config["rdna_mask_intervals"] = mask
        assert Util_funcs.quant_discard_reason(read) == "rdna_masked"

        # An explicit argument still wins over config, matching every other
        # threshold in this predicate: config says masked, the explicit
        # argument here says no mask at all, and the explicit one must win.
        assert Util_funcs.quant_discard_reason(read, rdna_mask={}) is None
    finally:
        LRAA_Globals.config["rdna_mask_intervals"] = saved
