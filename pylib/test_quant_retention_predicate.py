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
