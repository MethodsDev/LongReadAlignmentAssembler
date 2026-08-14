#!/usr/bin/env python3

import logging
import sys
from pathlib import Path

import pysam

REPO_ROOT = Path(__file__).resolve().parents[1]
for _import_dir in (str(REPO_ROOT / "pylib"), str(REPO_ROOT / "util")):
    if _import_dir not in sys.path:
        sys.path.insert(0, _import_dir)

import LRAA_Globals
import Util_funcs
import separate_bam_by_strand as sbs
from Bam_alignment_extractor import Bam_alignment_extractor


def _alignment(read_name, flag, ref_id, start):
    aln = pysam.AlignedSegment()
    aln.query_name = read_name
    aln.flag = flag
    aln.reference_id = ref_id
    aln.reference_start = start
    aln.mapping_quality = 60
    aln.cigar = [(0, 50)]
    aln.query_sequence = "A" * 50
    aln.query_qualities = pysam.qualitystring_to_array("I" * 50)
    aln.set_tag("NM", 0)
    aln.set_tag("AS", 50)
    return aln


def _write_test_bam(path):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as outf:
        outf.write(_alignment("read1", 0, 0, 100))
        outf.write(_alignment("read1", 256, 0, 150))
        outf.write(_alignment("read1", 2048, 0, 200))
    pysam.index(str(path))


def test_secondary_and_supplementary_alignments_are_always_discarded(tmp_path):
    bam_path = tmp_path / "input.bam"
    _write_test_bam(bam_path)

    read_flags = [
        read.flag
        for read in Bam_alignment_extractor(str(bam_path)).get_read_alignments("chr1")
    ]

    assert read_flags == [0]


def test_an_alignment_ending_at_the_region_start_is_retained(tmp_path):
    """`--region` boundaries are 1-based inclusive; pysam.fetch is 0-based half-open.

    Passing region_lend straight through started the window one base late, so an
    alignment whose only overlap with the region was its final base was silently
    dropped. Read A below ends exactly at 1-based 200 and read B starts at 201;
    a region of 200-300 must return both, and returned only B before the fix.

    This is reachable without chunking: any region-restricted run lost alignments
    ending at its start, and those runs get compared against whole-contig ones,
    which take the unrestricted fetch branch and lose nothing.
    """
    bam_path = tmp_path / "boundary.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        outf.write(_alignment("A", 0, 0, 150))  # 0-based 150-199 -> 1-based 151-200
        outf.write(_alignment("B", 0, 0, 200))  # 0-based 200-249 -> 1-based 201-250
    pysam.index(str(bam_path))

    names = [
        read.query_name
        for read in Bam_alignment_extractor(str(bam_path)).get_read_alignments(
            "chr1", region_lend=200, region_rend=300
        )
    ]

    assert names == ["A", "B"]


def test_discard_counters_name_the_reason_for_each_dropped_record(tmp_path, caplog):
    """The per-reason tally is the only place a dropped record is accounted for.

    Nothing returns the counters, so a run's sole evidence that 131k of 457k
    records were secondary rather than, say, below the mapping-quality floor is
    this one log line.
    """
    bam_path = tmp_path / "input.bam"
    _write_test_bam(bam_path)

    with caplog.at_level(logging.INFO, logger="Bam_alignment_extractor"):
        Bam_alignment_extractor(str(bam_path)).get_read_alignments("chr1")

    messages = " ".join(record.getMessage() for record in caplog.records)
    assert "'secondary': 1" in messages
    assert "'supplementary': 1" in messages


# Long-intron filtering during read INGESTION, ie. what read assignment may use.
#
# The same rule also runs at the strand split that opens depth normalization
# (util/separate_bam_by_strand.py), but that split only produces the splice-graph
# input, and it is skipped entirely under --no_norm or --normalize_max_cov_level 0.
# The check asserted here is on the path every quantified read takes, so it holds
# regardless of normalization settings.
#
# Threshold convention, shared with the splitter: the comparison is 'longer than',
# so an intron of exactly max_intron_length is KEPT and only max_intron_length + 1
# and above is discarded, matching pylib/Pretty_alignment.py:330,341 and
# pylib/Splice_graph.py:1515.

CONTIG = "chr1"
CONTIG_LENGTH = 3000000

MATCH = 0
INSERTION = 1
DELETION = 2
INTRON = 3
SOFTCLIP = 4

MAX_INTRON = LRAA_Globals.config["max_intron_length"]


def _spliced_alignment(header, read_name, cigartuples, reference_start=1000):
    aln = pysam.AlignedSegment(header)
    aln.query_name = read_name
    aln.flag = 0
    aln.reference_id = 0
    aln.reference_start = reference_start
    aln.mapping_quality = 60
    query_length = sum(
        length for op, length in cigartuples if op in (MATCH, INSERTION, SOFTCLIP)
    )
    aln.query_sequence = "A" * query_length
    aln.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    aln.cigartuples = list(cigartuples)
    aln.set_tag("NM", 0)
    return aln


def _run_extractor(tmp_path, cigars_by_read_name, caplog, max_intron_length=None):
    """extracts a synthetic bam; returns (retained read names, extractor log text)"""

    header = pysam.AlignmentHeader.from_references([CONTIG], [CONTIG_LENGTH])
    bam_path = tmp_path / "input.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as writer:
        for read_name, cigartuples in cigars_by_read_name.items():
            writer.write(_spliced_alignment(header, read_name, cigartuples))
    pysam.index(str(bam_path))

    original_max_intron_length = LRAA_Globals.config["max_intron_length"]
    if max_intron_length is not None:
        LRAA_Globals.config["max_intron_length"] = max_intron_length
    try:
        caplog.clear()
        with caplog.at_level(logging.INFO, logger="Bam_alignment_extractor"):
            reads = Bam_alignment_extractor(str(bam_path)).get_read_alignments(CONTIG)
    finally:
        LRAA_Globals.config["max_intron_length"] = original_max_intron_length

    log_text = " ".join(record.getMessage() for record in caplog.records)

    return {read.query_name for read in reads}, log_text


def test_read_with_long_internal_intron_is_discarded(tmp_path, caplog):
    """the long intron is neither the first nor the last, so terminal-intron
    trimming in Pretty_alignment could never have removed it"""

    retained, _log_text = _run_extractor(
        tmp_path,
        {
            "long_internal_intron": [
                (MATCH, 50),
                (INTRON, 1000),
                (MATCH, 50),
                (INTRON, MAX_INTRON + 1),
                (MATCH, 50),
                (INTRON, 1000),
                (MATCH, 50),
            ]
        },
        caplog,
    )

    assert retained == set()


def test_read_of_the_same_reference_span_without_a_long_intron_is_kept(
    tmp_path, caplog
):
    """span alone does not disqualify a read; only an individual overlong intron does"""

    long_intron = MAX_INTRON + 1
    first_half = long_intron // 2
    second_half = long_intron - first_half - 2
    control_cigar = [
        (MATCH, 50),
        (INTRON, first_half),
        (MATCH, 2),
        (INTRON, second_half),
        (MATCH, 50),
    ]
    discarded_cigar = [(MATCH, 50), (INTRON, long_intron), (MATCH, 50)]

    def _reference_span(cigartuples):
        return sum(
            length
            for op, length in cigartuples
            if op in (MATCH, DELETION, INTRON)
        )

    assert _reference_span(control_cigar) == _reference_span(discarded_cigar)

    retained, _log_text = _run_extractor(
        tmp_path,
        {"same_span_short_introns": control_cigar, "long_intron": discarded_cigar},
        caplog,
    )

    assert retained == {"same_span_short_introns"}


def test_read_with_long_terminal_intron_is_discarded(tmp_path, caplog):
    """a terminal overlong intron is otherwise merely trimmed off the alignment
    (Pretty_alignment.py:330,341); here the whole record goes"""

    retained, _log_text = _run_extractor(
        tmp_path,
        {
            "long_terminal_intron": [
                (MATCH, 50),
                (INTRON, MAX_INTRON + 1),
                (MATCH, 50),
                (INTRON, 1000),
                (MATCH, 50),
            ]
        },
        caplog,
    )

    assert retained == set()


def test_intron_of_exactly_max_intron_length_is_kept(tmp_path, caplog):
    retained, _log_text = _run_extractor(
        tmp_path,
        {
            "exactly_max": [(MATCH, 50), (INTRON, MAX_INTRON), (MATCH, 50)],
            "one_base_longer": [(MATCH, 50), (INTRON, MAX_INTRON + 1), (MATCH, 50)],
        },
        caplog,
    )

    assert retained == {"exactly_max"}


def test_long_deletion_is_not_treated_as_an_intron(tmp_path, caplog):
    """regression guarding the CIGAR-versus-get_blocks() decision.

    read.get_blocks() reports an identical gap for a deletion as for an intron of
    the same length, so a filter written against get_blocks() would discard this
    read.  The filter scans CIGAR N operations, which are the introns by
    specification, so a long genomic deletion is retained.
    """

    retained, log_text = _run_extractor(
        tmp_path,
        {
            "long_deletion": [
                (MATCH, 50),
                (DELETION, MAX_INTRON + 1),
                (MATCH, 50),
            ]
        },
        caplog,
    )

    assert retained == {"long_deletion"}
    assert "long_intron" not in log_text


def test_max_intron_length_of_zero_disables_intron_filtering(tmp_path, caplog):
    retained, _log_text = _run_extractor(
        tmp_path,
        {"very_long_intron": [(MATCH, 50), (INTRON, MAX_INTRON + 1), (MATCH, 50)]},
        caplog,
        max_intron_length=0,
    )

    assert retained == {"very_long_intron"}


def test_user_max_intron_length_override_reaches_the_extractor(tmp_path, caplog):
    """--max_intron_length lands in LRAA_Globals.config, which is where the
    extractor reads it, so no separate plumbing is needed"""

    retained, _log_text = _run_extractor(
        tmp_path,
        {
            "intron_501": [(MATCH, 50), (INTRON, 501), (MATCH, 50)],
            "intron_500": [(MATCH, 50), (INTRON, 500), (MATCH, 50)],
        },
        caplog,
        max_intron_length=500,
    )

    assert retained == {"intron_500"}


def test_discard_counter_names_the_long_intron_reason(tmp_path, caplog):
    """the counters are logged and are a run's only record of what left the
    pipeline, so the drop must be attributed rather than silent"""

    _retained, log_text = _run_extractor(
        tmp_path,
        {"long_intron": [(MATCH, 50), (INTRON, MAX_INTRON + 1), (MATCH, 50)]},
        caplog,
    )

    assert "'long_intron': 1" in log_text


def test_intron_length_rule_has_a_single_implementation(tmp_path, caplog):
    """The splice-graph input and read assignment must not be able to disagree
    about which reads carry a modellable intron, so both call sites resolve the
    same function object at call time; replacing it moves both.
    """

    header = pysam.AlignmentHeader.from_references([CONTIG], [CONTIG_LENGTH])
    unspliced_read = _spliced_alignment(header, "unspliced", [(MATCH, 100)])

    # the splitter re-exports the CIGAR scan rather than reimplementing it
    assert sbs.get_longest_intron_length is Util_funcs.get_longest_intron_length

    # and the rule is defined exactly once in the tree.  The needle is assembled so
    # that this file is not itself a match.
    needle = "def " + "get_longest_intron_length"
    definitions = sorted(
        path.relative_to(REPO_ROOT).as_posix()
        for path in REPO_ROOT.rglob("*.py")
        if needle in path.read_text(errors="replace")
    )
    assert definitions == ["pylib/Util_funcs.py"]

    # replacing the shared predicate must move both call sites, which it can only
    # do if neither holds its own copy of the rule
    original = Util_funcs.has_disqualifying_long_intron
    try:
        Util_funcs.has_disqualifying_long_intron = (
            lambda read, max_intron_length: True
        )
        assert sbs.get_discard_reason(unspliced_read, MAX_INTRON) == "long_intron"

        retained, log_text = _run_extractor(
            tmp_path, {"unspliced": [(MATCH, 100)]}, caplog
        )
        assert retained == set()
        assert "'long_intron': 1" in log_text
    finally:
        Util_funcs.has_disqualifying_long_intron = original

    assert sbs.get_discard_reason(unspliced_read, MAX_INTRON) is None
