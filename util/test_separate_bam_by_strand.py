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

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

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


# ---------------------------------------------------------------------------
# the per-contig fan-out over a whole bam
#
# The serial sweep is one thread whatever the machine, so the whole-bam callers
# -- stage 1 of the chunked pipeline and normalize_bam_by_strand.py -- pay one
# single-threaded pass over every record.  The fan-out reads each reference
# through the index instead and concatenates the parts.  It splits the bam every
# cluster's quantification builds its splice graph from, so the record set it
# produces has to be the sweep's record set and the counters it reports have to
# be the sweep's counters: nothing downstream would notice a quietly different
# splice graph, it would just be different.

FAN_OUT_REFERENCES = (
    ("chr1", 200000),
    ("chr2", 150000),
    ("chr3", 100000),
    ("chrM", 16569),
    # in the header, holding nothing: a whole-genome header carries hundreds of
    # these and the fan-out must not spend a pass on them
    ("chrEmpty", 50000),
)

POPULATED_REFERENCES = tuple(
    name for name, _length in FAN_OUT_REFERENCES if name != "chrEmpty"
)

# deliberately unequal, and NOT descending in header order: the fan-out hands the
# biggest scope to a worker first, so with these counts the order the passes are
# scheduled in (chr3, chr1, chr2, chrM) is not the order the parts must be
# concatenated in (chr1, chr2, chr3, chrM).  A concatenation that followed the
# work rather than the header would come out unsorted here.
RECORDS_PER_REFERENCE = {"chr1": 210, "chr2": 140, "chr3": 280, "chrM": 70}

PLACED_UNMAPPED_PER_REFERENCE = 3
NUM_UNPLACED = 5


def _fan_out_header():
    return pysam.AlignmentHeader.from_references(
        [name for name, _length in FAN_OUT_REFERENCES],
        [length for _name, length in FAN_OUT_REFERENCES],
    )


def _fan_out_record(
    header,
    read_name,
    reference_id,
    reference_start,
    flag=0,
    cigartuples=((MATCH, 100),),
):
    """one record, placed on the reference given rather than always on the first.

    An unmapped record keeps its reference_id when it has one: a mate-placed
    unmapped read carries RNAME and POS with FLAG 4, and a coordinate fetch does
    return it.  reference_id -1 is the coordinate-less kind, which no fetch by
    reference can reach.
    """

    aln = pysam.AlignedSegment(header)
    aln.query_name = read_name
    aln.flag = flag
    aln.reference_id = reference_id
    aln.reference_start = reference_start
    aln.mapping_quality = 0 if flag & 0x4 else 60
    aln.query_sequence = "A" * 100
    aln.query_qualities = pysam.qualitystring_to_array("I" * 100)
    aln.cigartuples = None if flag & 0x4 else list(cigartuples)
    return aln


def _write_fan_out_bam(path, with_unmapped=True):
    """a coordinate-sorted, indexed, multi-reference bam holding every discard reason.

    Every reason appears on every populated reference, so a merge that adopted
    one scope's counters instead of summing them is wrong in every counter rather
    than only in the total.
    """

    header = _fan_out_header()
    long_intron = [(MATCH, 50), (INTRON, MAX_INTRON + 1), (MATCH, 50)]
    placed = []

    for reference_id, (name, _length) in enumerate(FAN_OUT_REFERENCES):
        if name not in POPULATED_REFERENCES:
            continue

        for index in range(RECORDS_PER_REFERENCE[name]):
            reference_start = 1000 + index * 17
            flag, cigartuples = 0, ((MATCH, 100),)
            if index % 7 == 0:
                flag = 16
            elif index % 7 == 1:
                flag = 256
            elif index % 7 == 2:
                flag = 2048
            elif index % 7 == 3:
                flag = 1024
            elif index % 7 == 4:
                flag = 512
            elif index % 7 == 5:
                cigartuples = long_intron
            placed.append(
                (
                    reference_id,
                    reference_start,
                    _fan_out_record(
                        header,
                        "{}_{}".format(name, index),
                        reference_id,
                        reference_start,
                        flag=flag,
                        cigartuples=cigartuples,
                    ),
                )
            )

        if with_unmapped:
            for which in range(PLACED_UNMAPPED_PER_REFERENCE):
                reference_start = 900 + which
                placed.append(
                    (
                        reference_id,
                        reference_start,
                        _fan_out_record(
                            header,
                            "{}_placed_unmapped_{}".format(name, which),
                            reference_id,
                            reference_start,
                            flag=4,
                        ),
                    )
                )

    placed.sort(key=lambda entry: (entry[0], entry[1]))

    with pysam.AlignmentFile(str(path), "wb", header=header) as writer:
        for _reference_id, _reference_start, record in placed:
            writer.write(record)
        if with_unmapped:
            # coordinate-less records live at the end of a coordinate-sorted bam
            for which in range(NUM_UNPLACED):
                writer.write(
                    _fan_out_record(header, "unplaced_{}".format(which), -1, -1, flag=4)
                )

    pysam.index(str(path))
    return path


def _record_identities(bam_path):
    """name, flag, reference and position of every record written.

    Not the bam's bytes: the fan-out concatenates per-reference parts, so two
    equal record sets may legitimately differ in where compression blocks fall.
    """

    with pysam.AlignmentFile(str(bam_path), "rb") as reader:
        return sorted(
            (read.query_name, read.flag, read.reference_id, read.reference_start)
            for read in reader.fetch(until_eof=True)
        )


def _split_both_ways(tmp_path, bam, num_workers=4, max_intron_length=MAX_INTRON):
    """the sweep and the fan-out over one bam: {arm: (counters, output paths)}"""

    results = dict()

    for arm in ("sweep", "fan_out"):
        outputs = tuple(
            str(tmp_path / "{}.{}.bam".format(arm, strand)) for strand in ("+", "-")
        )
        if arm == "sweep":
            counters = sbs.split_bam_by_strand(
                str(bam), outputs[0], outputs[1], max_intron_length
            )
        else:
            counters = sbs.split_bam_by_strand_parallel(
                str(bam),
                outputs[0],
                outputs[1],
                max_intron_length,
                num_workers=num_workers,
            )
        results[arm] = (counters, outputs)

    return results


def test_fan_out_writes_the_same_records_the_sweep_writes(tmp_path):
    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    results = _split_both_ways(tmp_path, bam)

    sweep_counters, sweep_outputs = results["sweep"]
    fan_out_counters, fan_out_outputs = results["fan_out"]

    # the fixture has to be able to fail this: several populated references, both
    # orientations, and every discard reason that a valid bam can carry
    assert sweep_counters["num_forward"] > 0
    assert sweep_counters["num_reverse"] > 0
    for reason in sbs.DISCARD_REASONS:
        if reason == "no_chromosome":
            # a record with no reference is unmapped in any valid bam, and
            # "unmapped" is tested first
            continue
        assert sweep_counters[sbs.discard_counter_name(reason)] > 0, reason

    for strand in (0, 1):
        assert _record_identities(fan_out_outputs[strand]) == _record_identities(
            sweep_outputs[strand]
        )

    assert fan_out_counters == sweep_counters
    assert sbs.build_report(fan_out_counters) == sbs.build_report(sweep_counters)


def test_fan_out_keeps_unmapped_records_in_the_denominator(tmp_path):
    """Both kinds of unmapped record, counted and then discarded as the sweep does.

    A coordinate fetch cannot return a coordinate-less record, so a fan-out that
    visited only the references would drop them from num_records -- and
    build_report divides every percentage by that count, so the report would read
    as a complete accounting while describing a smaller file.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    num_in_file = int(
        subprocess.check_output(["samtools", "view", "-c", str(bam)], text=True)
    )
    num_unmapped = (
        len(POPULATED_REFERENCES) * PLACED_UNMAPPED_PER_REFERENCE + NUM_UNPLACED
    )

    results = _split_both_ways(tmp_path, bam)
    for arm in ("sweep", "fan_out"):
        counters = results[arm][0]
        assert counters["num_records"] == num_in_file, arm
        assert counters[sbs.discard_counter_name("unmapped")] == num_unmapped, arm

    scopes, total = sbs.index_record_counts(str(bam))
    assert total == num_in_file
    assert sum(count for _scope, count in scopes) == num_in_file
    assert (sbs.UNPLACED_SCOPE, NUM_UNPLACED) in scopes
    # header order, which is the order the parts are concatenated in, and no pass
    # spent on a reference holding nothing
    assert [scope for scope, _count in scopes] == list(POPULATED_REFERENCES) + [
        sbs.UNPLACED_SCOPE
    ]


def test_a_bam_with_no_unmapped_records_needs_no_unplaced_pass(tmp_path):
    bam = _write_fan_out_bam(tmp_path / "mapped_only.bam", with_unmapped=False)
    scopes, _total = sbs.index_record_counts(str(bam))

    assert sbs.UNPLACED_SCOPE not in [scope for scope, _count in scopes]

    results = _split_both_ways(tmp_path, bam)
    assert results["fan_out"][0] == results["sweep"][0]
    assert results["fan_out"][0][sbs.discard_counter_name("unmapped")] == 0


def test_merge_counters_sums_every_counter_rather_than_adopting_one(tmp_path):
    """The negative control, stated on the merge itself.

    A per-scope dict carries that scope's records, so keeping one of them -- the
    last to finish, say, which is right there when the loop ends -- leaves every
    counter at one reference's value while the report presents it as the file's.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    per_scope = [
        sbs.split_bam_by_strand(
            str(bam),
            str(tmp_path / "scope_{}.+.bam".format(index)),
            str(tmp_path / "scope_{}.-.bam".format(index)),
            MAX_INTRON,
            contig=scope,
        )
        for index, (scope, _count) in enumerate(sbs.index_record_counts(str(bam))[0])
    ]
    assert len(per_scope) > 1

    merged = sbs.merge_counters(per_scope)

    for counter_name in sbs.new_counters():
        assert merged[counter_name] == sum(
            counters[counter_name] for counters in per_scope
        ), counter_name

    # the negative control: adopting one scope's dict differs in the denominator
    # and in the per-reason breakdown alike
    adopted = per_scope[-1]
    assert adopted["num_records"] != merged["num_records"]
    for counter_name in (
        "num_forward",
        "num_reverse",
        sbs.discard_counter_name("secondary"),
        sbs.discard_counter_name("long_intron"),
    ):
        assert adopted[counter_name] != merged[counter_name], counter_name


def test_a_wrong_counter_merge_is_refused_rather_than_reported(tmp_path, monkeypatch):
    """The same negative control end to end: the fan-out must not ship the bug.

    With the merge replaced by exactly that mistake, the run has to fail. The
    guard is the index's own record total, checked against the merged accounting
    before any output is concatenated.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    monkeypatch.setattr(sbs, "merge_counters", lambda per_scope: dict(per_scope[-1]))

    with pytest.raises(RuntimeError) as raised:
        sbs.split_bam_by_strand_parallel(
            str(bam),
            str(tmp_path / "out.+.bam"),
            str(tmp_path / "out.-.bam"),
            MAX_INTRON,
            num_workers=4,
        )

    assert str(sbs.index_record_counts(str(bam))[1]) in str(raised.value)


def test_merge_counters_refuses_a_foreign_accounting():
    counters = sbs.new_counters()
    del counters["num_forward"]

    with pytest.raises(ValueError):
        sbs.merge_counters([sbs.new_counters(), counters])


def test_fan_out_output_can_be_indexed_and_reports_what_the_sweep_reports(tmp_path):
    """End to end through main(), which indexes both outputs.

    `samtools index` refuses a bam whose records are not in coordinate order, so
    parts concatenated in completion order rather than header order fail here.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    reports = dict()

    for arm, workers in (("sweep", "1"), ("fan_out", "4")):
        prefix = str(tmp_path / arm)
        completed = subprocess.run(
            [
                sys.executable,
                str(SPLITTER),
                "--bam",
                str(bam),
                "--output_prefix",
                prefix,
                "--max_intron_length",
                str(MAX_INTRON),
                "--num_workers",
                workers,
            ],
            capture_output=True,
            text=True,
        )

        assert completed.returncode == 0, completed.stderr
        for strand in ("+", "-"):
            assert Path("{}.{}.bam.bai".format(prefix, strand)).exists()
            with pysam.AlignmentFile("{}.{}.bam".format(prefix, strand), "rb") as reader:
                positions = [
                    (read.reference_id, read.reference_start)
                    for read in reader.fetch(until_eof=True)
                ]
            assert positions == sorted(positions), (arm, strand)

        reports[arm] = completed.stderr[
            completed.stderr.index("Num input bam records:") :
        ]

    assert reports["fan_out"] == reports["sweep"]
    assert "Num records discarded as unmapped:" in reports["fan_out"]


def test_an_unindexed_input_is_indexed_rather_than_swept(tmp_path):
    """The fan-out addresses a reference by name, so it needs the index.

    Neither whole-bam caller guarantees one -- a WDL task localizes no .bai it
    did not declare -- and the alternative to building it is the single-threaded
    pass this exists to remove.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    unindexed = tmp_path / "unindexed.bam"
    unindexed.write_bytes(bam.read_bytes())
    assert not Path(str(unindexed) + ".bai").exists()

    completed = subprocess.run(
        [
            sys.executable,
            str(SPLITTER),
            "--bam",
            str(unindexed),
            "--output_prefix",
            str(tmp_path / "built"),
            "--max_intron_length",
            str(MAX_INTRON),
        ],
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert Path(str(unindexed) + ".bai").exists()
    assert "has no index" in completed.stderr

    sweep = _split_both_ways(tmp_path, bam)["sweep"][1]
    for strand, swept in zip(("+", "-"), sweep):
        assert _record_identities(
            str(tmp_path / "built.{}.bam".format(strand))
        ) == _record_identities(swept)


@pytest.mark.skipif(
    os.geteuid() == 0, reason="root writes into a read-only directory anyway"
)
def test_an_unindexable_input_is_refused_rather_than_swept_silently(tmp_path):
    """Refusal, not a quiet hour on one core.

    A caller that means to take the sweep says so with --num_workers 1, and the
    refusal names that flag.
    """

    read_only = tmp_path / "read_only"
    read_only.mkdir()
    bam = _write_fan_out_bam(read_only / "multi.bam")
    Path(str(bam) + ".bai").unlink()
    read_only.chmod(0o500)

    try:
        completed = subprocess.run(
            [
                sys.executable,
                str(SPLITTER),
                "--bam",
                str(bam),
                "--output_prefix",
                str(tmp_path / "refused"),
                "--max_intron_length",
                str(MAX_INTRON),
            ],
            capture_output=True,
            text=True,
        )
    finally:
        read_only.chmod(0o700)

    assert completed.returncode != 0
    assert "--num_workers 1" in completed.stderr
    assert not (tmp_path / "refused.+.bam").exists()


def test_one_worker_is_the_unchanged_whole_file_sweep(tmp_path):
    """The documented way out, and it may not quietly become the fan-out."""

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    completed = subprocess.run(
        [
            sys.executable,
            str(SPLITTER),
            "--bam",
            str(bam),
            "--output_prefix",
            str(tmp_path / "swept"),
            "--max_intron_length",
            str(MAX_INTRON),
            "--num_workers",
            "1",
        ],
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr
    assert "single whole-file pass" in completed.stderr
    assert "scope(s) on" not in completed.stderr


def test_a_bam_holding_no_records_still_produces_both_outputs(tmp_path):
    bam = tmp_path / "empty.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=_fan_out_header()):
        pass
    pysam.index(str(bam))

    outputs = [str(tmp_path / "empty.{}.bam".format(strand)) for strand in ("+", "-")]
    counters = sbs.split_bam_by_strand_parallel(
        str(bam), outputs[0], outputs[1], MAX_INTRON, num_workers=4
    )

    assert counters == sbs.new_counters()
    for output in outputs:
        with pysam.AlignmentFile(output, "rb") as reader:
            assert list(reader.references) == [
                name for name, _length in FAN_OUT_REFERENCES
            ]
            assert not list(reader.fetch(until_eof=True))


def test_workers_resolve_against_the_work_and_the_cpus_available():
    assert sbs.resolve_num_workers(4, num_tasks=10) == 4
    # never more processes than there are scopes to read
    assert sbs.resolve_num_workers(32, num_tasks=3) == 3
    # 0 means "this machine", bounded the same way
    assert sbs.resolve_num_workers(0, num_tasks=1) == 1
    assert sbs.resolve_num_workers(0, num_tasks=10**6) == sbs.available_cpu_count()
    assert sbs.resolve_num_workers(-1) == sbs.available_cpu_count()


def test_help_reports_the_worker_count_with_its_default():
    help_text = subprocess.run(
        [sys.executable, str(SPLITTER), "--help"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    assert "--num_workers" in help_text
    assert "one worker per CPU" in help_text



# ---------------------------------------------------------------------------
# coordinate order of the concatenated output
#
# Record-set equality cannot see this, and it is the property a concatenation is
# most likely to lose: the parts finish in whatever order the workers get through
# them.  The consumers of these two bams do care.  Cut selection fetches them BY
# REGION (util/misc/select_contig_cut_points.py:569,624,739), where a record in
# the wrong place is not an error, it is simply missing from the answer.  Depth
# normalization, the immediate consumer, anchors each contig's window grid on the
# first aligned base it sees on that contig and REFUSES any later record starting
# before it (normalize_bam_by_strand.py:_reject_read_before_anchor) -- so out of
# order the run either dies there or, given --window_origin, silently measures
# depth over a different grid.


def _positions(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb") as reader:
        return [
            (read.reference_id, read.reference_start)
            for read in reader.fetch(until_eof=True)
        ]


def test_fan_out_output_is_coordinate_sorted_in_the_inputs_reference_order(tmp_path):
    """Walk the records: (reference_id, reference_start) may never decrease.

    reference_id ascending IS the input's header order, since both outputs are
    written from the input's header; so this asserts the contig order as well as
    the positions within each contig.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    results = _split_both_ways(tmp_path, bam)

    for output in results["fan_out"][1]:
        positions = _positions(output)
        assert positions, output
        assert positions == sorted(positions), output
        # and the references appear in header order, each one contiguously
        reference_ids = [reference_id for reference_id, _start in positions]
        assert reference_ids == sorted(reference_ids)

    # the scheduling order really is not the header order, so header-order
    # concatenation is doing work here rather than coinciding with completion
    scopes = sbs.index_record_counts(str(bam))[0]
    scheduled = [scope for scope, _count in sorted(scopes, key=lambda s: -s[1])]
    assert scheduled != [scope for scope, _count in scopes]


def test_fan_out_output_survives_quickcheck_and_an_indexed_late_contig_fetch(tmp_path):
    """The consumers' own access pattern, on the contig a wrong order breaks first.

    chrM is last in the header, so a concatenation in completion order puts its
    records ahead of records that sort before them; an index built over that
    either fails outright or answers a region query with the wrong records.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    prefix = str(tmp_path / "fanned")
    completed = subprocess.run(
        [
            sys.executable,
            str(SPLITTER),
            "--bam",
            str(bam),
            "--output_prefix",
            prefix,
            "--max_intron_length",
            str(MAX_INTRON),
            "--num_workers",
            "4",
        ],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    late_reference = POPULATED_REFERENCES[-1]
    for strand in ("+", "-"):
        output = "{}.{}.bam".format(prefix, strand)

        assert Path(output + ".bai").exists()
        subprocess.run(
            ["samtools", "quickcheck", "-v", output], check=True, capture_output=True
        )

        with pysam.AlignmentFile(output, "rb") as reader:
            fetched = sorted(
                read.query_name for read in reader.fetch(late_reference, 900, 3000)
            )
        # a fresh handle: an until_eof iteration resumes from wherever the region
        # iterator left the file, and would otherwise report nothing
        with pysam.AlignmentFile(output, "rb") as reader:
            expected = sorted(
                read.query_name
                for read in reader.fetch(until_eof=True)
                if read.reference_name == late_reference
                and 900 <= read.reference_start < 3000
            )

        # the fetch has to be able to return something, or it proves nothing
        assert fetched
        assert all(name.startswith(late_reference) for name in fetched)
        assert fetched == expected


def test_concatenating_the_parts_out_of_order_is_what_verify_part_order_catches(
    tmp_path,
):
    """The negative control for the ordering guarantee.

    Same parts, concatenated in reverse header order -- what a concatenation
    driven by completion order can produce.  The record SET is unchanged, which is
    the point: record-set equality cannot see this at all.

    And neither can samtools: asserted below, `samtools index` builds an index
    over the reversed file with exit 0.  That is why the guarantee is asserted
    against the parts, before they are joined, rather than left to the index.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam", with_unmapped=False)
    scopes = sbs.index_record_counts(str(bam))[0]

    tasks = []
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        for index, (scope, count) in enumerate(scopes):
            task = {
                "scope": scope,
                "num_records": count,
                "reference_id": reader.get_tid(scope),
                "+": str(tmp_path / "part_{}.+.bam".format(index)),
                "-": str(tmp_path / "part_{}.-.bam".format(index)),
            }
            sbs.split_bam_by_strand(
                str(bam), task["+"], task["-"], MAX_INTRON, contig=scope
            )
            tasks.append(task)

    in_order = str(tmp_path / "in_order.bam")
    reversed_order = str(tmp_path / "reversed.bam")
    sbs.concatenate_bam_parts(in_order, [task["+"] for task in tasks])
    sbs.concatenate_bam_parts(
        reversed_order, [task["+"] for task in reversed(tasks)]
    )

    # identical record sets: exactly what record-set equality misses
    assert _record_identities(reversed_order) == _record_identities(in_order)

    assert _positions(in_order) == sorted(_positions(in_order))
    assert _positions(reversed_order) != sorted(_positions(reversed_order))

    # the guard fires on the reversed order and stays quiet on the right one
    sbs.verify_part_order(tasks, "+")
    with pytest.raises(RuntimeError) as raised:
        sbs.verify_part_order(list(reversed(tasks)), "+")
    assert "not in header order" in str(raised.value)

    # and the reason the guard cannot be left to samtools
    indexed = subprocess.run(
        ["samtools", "index", reversed_order], capture_output=True, text=True
    )
    assert indexed.returncode == 0, "if samtools started refusing this, simplify"


def test_a_part_holding_the_wrong_reference_is_refused(tmp_path):
    """The other half of the guard: a part that is not the scope it claims.

    Cheap to check and worth checking, because the part paths are derived from a
    scope's position rather than its name.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam", with_unmapped=False)
    scopes = sbs.index_record_counts(str(bam))[0]

    task = {
        "scope": scopes[0][0],
        "num_records": scopes[0][1],
        # deliberately wrong: the part will hold chr1 records
        "reference_id": 2,
        "+": str(tmp_path / "mismatched.+.bam"),
        "-": str(tmp_path / "mismatched.-.bam"),
    }
    sbs.split_bam_by_strand(
        str(bam), task["+"], task["-"], MAX_INTRON, contig=scopes[0][0]
    )

    with pytest.raises(RuntimeError) as raised:
        sbs.verify_part_order([task], "+")

    assert "rather than 2" in str(raised.value)


def test_an_input_declaring_a_non_coordinate_sort_order_is_refused(tmp_path):
    """The fan-out's order is the input's order, so it will not guess.

    A queryname-sorted input would give per-reference parts that are not sorted
    within themselves, which nothing downstream would report.
    """

    bam = _write_fan_out_bam(tmp_path / "multi.bam")
    header = dict(_fan_out_header().to_dict())
    header["HD"] = {"VN": "1.6", "SO": "queryname"}

    mislabelled = tmp_path / "queryname.bam"
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        with pysam.AlignmentFile(str(mislabelled), "wb", header=header) as writer:
            for read in reader.fetch(until_eof=True):
                writer.write(read)
    shutil.copy(str(bam) + ".bai", str(mislabelled) + ".bai")

    with pytest.raises(RuntimeError) as raised:
        sbs.split_bam_by_strand_parallel(
            str(mislabelled),
            str(tmp_path / "out.+.bam"),
            str(tmp_path / "out.-.bam"),
            MAX_INTRON,
            num_workers=4,
        )

    assert "queryname" in str(raised.value)
    assert "--num_workers 1" in str(raised.value)