#!/usr/bin/env python3

"""Record filtering and accounting of util/separate_bam_by_strand.py.

The splitter is the pipeline's entry point: whatever it drops is gone from every
downstream stage, including depth normalization, so the drops are asserted here
record by record and the counters are asserted to account for every input record.

Threshold convention asserted below: the comparison is 'longer than', so an
intron of exactly max_intron_length is KEPT and only max_intron_length + 1 and
above is discarded. This matches the existing enforcement points in
pylib/Pretty_alignment.py and pylib/Splice_graph.py, which both discard on
intron_len > LRAA_Globals.config["max_intron_length"].
"""

import subprocess
import sys
from pathlib import Path

import pysam

UTIL_DIR = Path(__file__).resolve().parent
REPO_ROOT = UTIL_DIR.parent
for _import_dir in (str(UTIL_DIR), str(REPO_ROOT / "pylib")):
    if _import_dir not in sys.path:
        sys.path.insert(0, _import_dir)

import LRAA_Globals
import separate_bam_by_strand as sbs

SPLITTER = UTIL_DIR / "separate_bam_by_strand.py"

CONTIG = "chr1"
CONTIG_LENGTH = 3000000

MATCH = 0
INSERTION = 1
DELETION = 2
INTRON = 3
SOFTCLIP = 4

MAX_INTRON = 1000


def _header():
    return pysam.AlignmentHeader.from_references([CONTIG], [CONTIG_LENGTH])


def _alignment(header, read_name, cigartuples, flag=0, reference_start=1000):
    aln = pysam.AlignedSegment(header)
    aln.query_name = read_name
    aln.flag = flag
    aln.reference_id = -1 if flag & 0x4 else 0
    aln.reference_start = reference_start
    aln.mapping_quality = 60
    query_length = sum(
        length
        for op, length in cigartuples
        if op in (MATCH, INSERTION, SOFTCLIP)
    )
    aln.query_sequence = "A" * query_length
    aln.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    aln.cigartuples = list(cigartuples)
    return aln


def _write_bam(path, alignments):
    header = _header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as writer:
        for make_alignment in alignments:
            writer.write(make_alignment(header))
    return path


def _run_splitter(tmp_path, alignments, max_intron_length=MAX_INTRON):
    """runs the split over a synthetic bam; returns (retained read names, counters)"""

    input_bam = _write_bam(tmp_path / "input.bam", alignments)
    top_bam = tmp_path / "out.+.bam"
    bottom_bam = tmp_path / "out.-.bam"

    counters = sbs.split_bam_by_strand(
        str(input_bam), str(top_bam), str(bottom_bam), max_intron_length
    )

    retained = set()
    for output_bam in (top_bam, bottom_bam):
        with pysam.AlignmentFile(str(output_bam), "rb") as reader:
            retained.update(read.query_name for read in reader)

    return retained, counters


def _spliced(intron_length, num_leading_introns=0):
    """cigar of an alignment whose longest intron is intron_length.

    num_leading_introns short introns are placed ahead of it, so the long intron
    is internal rather than terminal.
    """

    cigartuples = []
    for _ in range(num_leading_introns):
        cigartuples += [(MATCH, 100), (INTRON, 50)]
    cigartuples += [(MATCH, 100), (INTRON, intron_length), (MATCH, 100)]
    return cigartuples


def _cigar_of_span(target_span, intron_length):
    """cigar spanning target_span reference bases using introns of intron_length"""

    cigartuples = []
    span = 0
    while target_span - span > intron_length + 200:
        cigartuples += [(MATCH, 100), (INTRON, intron_length)]
        span += 100 + intron_length
    cigartuples += [(MATCH, target_span - span)]

    return cigartuples


def _reference_span(cigartuples):
    return sum(
        length
        for op, length in cigartuples
        if op in (MATCH, DELETION, INTRON)
    )


def _primary(header, read_name, cigartuples=((MATCH, 100),)):
    return _alignment(header, read_name, cigartuples)


def test_secondary_record_is_discarded(tmp_path):
    retained, counters = _run_splitter(
        tmp_path,
        [
            lambda h: _primary(h, "keeper"),
            lambda h: _alignment(h, "secondary", [(MATCH, 100)], flag=256),
        ],
    )

    assert retained == {"keeper"}
    assert counters["num_records_discarded_secondary"] == 1


def test_supplementary_record_is_discarded(tmp_path):
    retained, counters = _run_splitter(
        tmp_path,
        [
            lambda h: _primary(h, "keeper"),
            lambda h: _alignment(h, "supplementary", [(MATCH, 100)], flag=2048),
        ],
    )

    assert retained == {"keeper"}
    assert counters["num_records_discarded_supplementary"] == 1


def test_duplicate_and_qcfail_records_are_discarded(tmp_path):
    retained, counters = _run_splitter(
        tmp_path,
        [
            lambda h: _primary(h, "keeper"),
            lambda h: _alignment(h, "duplicate", [(MATCH, 100)], flag=1024),
            lambda h: _alignment(h, "qcfail", [(MATCH, 100)], flag=512),
        ],
    )

    assert retained == {"keeper"}
    assert counters["num_records_discarded_duplicate"] == 1
    assert counters["num_records_discarded_qcfail"] == 1


def test_long_internal_intron_is_discarded_while_equal_span_is_kept(tmp_path):
    # both reads span the same reference distance; only one of them gets there by
    # way of a single long intron, and that intron is internal to the alignment
    long_internal = _spliced(MAX_INTRON + 1, num_leading_introns=2)
    equal_span = _cigar_of_span(_reference_span(long_internal), MAX_INTRON // 2)

    assert _reference_span(equal_span) == _reference_span(long_internal)
    assert max(length for op, length in equal_span if op == INTRON) <= MAX_INTRON

    retained, counters = _run_splitter(
        tmp_path,
        [
            lambda h: _alignment(h, "long_internal_intron", long_internal),
            lambda h: _alignment(h, "same_span_no_long_intron", equal_span),
        ],
    )

    assert retained == {"same_span_no_long_intron"}
    assert counters["num_records_discarded_long_intron"] == 1


def test_long_terminal_intron_is_also_discarded(tmp_path):
    retained, counters = _run_splitter(
        tmp_path,
        [
            # long intron is the leftmost/rightmost one, ie. the case the
            # terminal-intron pruning elsewhere in LRAA already trims
            lambda h: _alignment(h, "long_terminal_intron", _spliced(MAX_INTRON + 1)),
            lambda h: _primary(h, "keeper", _spliced(MAX_INTRON - 1)),
        ],
    )

    assert retained == {"keeper"}
    assert counters["num_records_discarded_long_intron"] == 1


def test_intron_exactly_at_threshold_is_kept(tmp_path):
    # 'longer than' comparison: max_intron_length itself passes, one more does not
    retained, counters = _run_splitter(
        tmp_path,
        [
            lambda h: _alignment(h, "at_threshold", _spliced(MAX_INTRON)),
            lambda h: _alignment(h, "one_over_threshold", _spliced(MAX_INTRON + 1)),
        ],
    )

    assert retained == {"at_threshold"}
    assert counters["num_records_discarded_long_intron"] == 1


def test_zero_max_intron_length_disables_intron_filtering(tmp_path):
    long_intron_cigar = _spliced(500000, num_leading_introns=1)

    retained, counters = _run_splitter(
        tmp_path,
        [lambda h: _alignment(h, "very_long_intron", long_intron_cigar)],
        max_intron_length=0,
    )

    assert retained == {"very_long_intron"}
    assert counters["num_records_discarded_long_intron"] == 0

    # ... and the same read is discarded once the filter is in force
    retained, counters = _run_splitter(
        tmp_path,
        [lambda h: _alignment(h, "very_long_intron", long_intron_cigar)],
        max_intron_length=200000,
    )

    assert retained == set()
    assert counters["num_records_discarded_long_intron"] == 1


def test_negative_max_intron_length_disables_intron_filtering(tmp_path):
    retained, _counters = _run_splitter(
        tmp_path,
        [lambda h: _alignment(h, "very_long_intron", _spliced(500000))],
        max_intron_length=-1,
    )

    assert retained == {"very_long_intron"}


def test_deletion_is_not_treated_as_an_intron(tmp_path):
    # read.get_blocks() reports the same gap for a deletion as for an intron of
    # the same length; the filter reads CIGAR N operations, so a long deletion is
    # not mistaken for a long intron
    retained, counters = _run_splitter(
        tmp_path,
        [lambda h: _alignment(h, "long_deletion", [(MATCH, 100), (DELETION, MAX_INTRON + 1), (MATCH, 100)])],
    )

    assert retained == {"long_deletion"}
    assert counters["num_records_discarded_long_intron"] == 0


def test_counters_account_for_every_record_read(tmp_path):
    alignments = [
        lambda h: _primary(h, "fwd"),
        lambda h: _alignment(h, "rev", [(MATCH, 100)], flag=16),
        lambda h: _alignment(h, "secondary", [(MATCH, 100)], flag=256),
        lambda h: _alignment(h, "supplementary", [(MATCH, 100)], flag=2048),
        lambda h: _alignment(h, "duplicate", [(MATCH, 100)], flag=1024),
        lambda h: _alignment(h, "qcfail", [(MATCH, 100)], flag=512),
        lambda h: _alignment(h, "unmapped", [(MATCH, 100)], flag=4),
        lambda h: _alignment(h, "long_intron", _spliced(MAX_INTRON + 1, 1)),
    ]

    retained, counters = _run_splitter(tmp_path, alignments)

    assert retained == {"fwd", "rev"}
    assert counters["num_records"] == len(alignments)

    num_discarded = sum(
        counters[sbs.discard_counter_name(reason)] for reason in sbs.DISCARD_REASONS
    )
    num_written = counters["num_forward"] + counters["num_reverse"]

    assert num_discarded == 6
    assert num_written == 2
    assert (
        num_discarded + num_written + counters["num_neither"]
        == counters["num_records"]
    )

    # every count is reported against the number of records read
    report = sbs.build_report(counters)
    for reason in sbs.DISCARD_REASONS:
        assert "Num records discarded as {}:".format(reason) in report
    assert "of {}".format(counters["num_records"]) in report


def test_longest_intron_length_finds_internal_introns():
    header = _header()
    read = _alignment(header, "internal", _spliced(7777, num_leading_introns=3))

    assert sbs.get_longest_intron_length(read) == 7777
    assert sbs.get_longest_intron_length(_alignment(header, "unspliced", [(MATCH, 100)])) == 0


def test_help_reports_the_new_options_with_defaults():
    help_text = subprocess.run(
        [sys.executable, str(SPLITTER), "--help"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    assert "--max_intron_length" in help_text
    assert str(LRAA_Globals.config["max_intron_length"]) in help_text
    assert "disable intron length filtering" in help_text
