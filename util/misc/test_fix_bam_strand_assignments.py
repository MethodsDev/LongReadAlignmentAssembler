#!/usr/bin/env python3

"""What fix_bam_strand_assignments.py promises about the bam it writes.

Two claims are tested, and they pull in opposite directions.  The tool must
CHANGE the strand of every record whose introns contradict its flag, and it
must change NOTHING ELSE: not the record count, not the order, not a
coordinate, not a record it had no evidence about.  A fixture where no read
flips would pass the second claim while testing nothing, so the flips are
asserted by name.
"""

import os
import subprocess
import sys

import pysam
import pytest

UTIL_MISC_DIR = os.path.dirname(os.path.realpath(__file__))
SCRIPT = os.path.join(UTIL_MISC_DIR, "fix_bam_strand_assignments.py")

sys.path.insert(0, UTIL_MISC_DIR)

import fix_bam_strand_assignments as fixer


CONTIG = "chrT"
LENGTH = 30000

MATCH = 0
REF_SKIP = 3

# 1-based inclusive intron coordinates -> dinucleotides planted there.
FORWARD = ("GT", "AG")  # canonical on the top strand
REVERSE = ("CT", "AC")  # its reverse complement
NONCANONICAL = ("AA", "CC")


class Corpus:
    """A contig whose splice sites are chosen rather than incidental."""

    def __init__(self, tmp_path):
        self.dir = tmp_path
        self.sequence = list(("ACGT" * (LENGTH // 4 + 1))[:LENGTH])
        self.reads = []  # (name, flag, blocks, tags)

    def read(self, name, flag, blocks, dinucs=(), tags=()):
        for i, dinuc_pair in enumerate(dinucs, start=1):
            intron_lend = blocks[i - 1][1] + 1
            intron_rend = blocks[i][0] - 1
            self.sequence[intron_lend - 1] = dinuc_pair[0][0]
            self.sequence[intron_lend] = dinuc_pair[0][1]
            self.sequence[intron_rend - 2] = dinuc_pair[1][0]
            self.sequence[intron_rend - 1] = dinuc_pair[1][1]
        self.reads.append((name, flag, blocks, tuple(tags)))
        return self

    def unmapped_read(self, name):
        self.reads.append((name, 4, None, ()))
        return self

    def build(self):
        self.fasta = str(self.dir / "genome.fa")
        sequence = "".join(self.sequence)
        with open(self.fasta, "wt") as ofh:
            print(">{}".format(CONTIG), file=ofh)
            for i in range(0, LENGTH, 60):
                print(sequence[i : i + 60], file=ofh)
        pysam.faidx(self.fasta)

        self.bam = str(self.dir / "reads.bam")
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": LENGTH}],
        }
        with pysam.AlignmentFile(self.bam, "wb", header=header) as ofh:
            for read in self.reads:
                ofh.write(self._alignment(*read))
        return self

    def _alignment(self, name, flag, blocks, tags):
        aln = pysam.AlignedSegment()
        aln.query_name = name
        aln.flag = flag
        if blocks is None:
            aln.reference_id = -1
            aln.reference_start = -1
            aln.query_sequence = "A" * 50
            aln.query_qualities = pysam.qualitystring_to_array("I" * 50)
            return aln
        aln.reference_id = 0
        aln.reference_start = blocks[0][0] - 1
        aln.mapping_quality = 60
        cigar = []
        for i, (lend, rend) in enumerate(blocks):
            if i:
                cigar.append((REF_SKIP, lend - blocks[i - 1][1] - 1))
            cigar.append((MATCH, rend - lend + 1))
        aln.cigartuples = cigar
        length = sum(rend - lend + 1 for lend, rend in blocks)
        aln.query_sequence = "A" * length
        aln.query_qualities = pysam.qualitystring_to_array("I" * length)
        aln.set_tag("NM", 0)
        for tag, value, value_type in tags:
            aln.set_tag(tag, value, value_type=value_type)
        return aln


@pytest.fixture
def corpus(tmp_path):
    """One read per branch of the decision, in coordinate order."""

    return (
        Corpus(tmp_path)
        # aligner says forward, the splice sites say reverse
        .read("flip_to_minus", 0, [(2000, 2200), (2600, 3000)], [REVERSE])
        # aligner says reverse, the splice sites say forward; carries minimap2's ts
        .read(
            "flip_to_plus",
            16,
            [(4000, 4200), (4600, 5000)],
            [FORWARD],
            tags=[("ts", "+", "A")],
        )
        # aligner and splice sites agree
        .read("agree_plus", 0, [(6000, 6200), (6600, 7000)], [FORWARD])
        # a secondary alignment is evidence about itself, and is corrected too
        .read("secondary_flip", 16 | 256, [(8000, 8200), (8600, 9000)], [FORWARD])
        # spliced, but the motifs are not canonical in either orientation
        .read("noncanonical", 0, [(10000, 10200), (10600, 11000)], [NONCANONICAL])
        # one intron votes each way: no majority, so nothing moves
        .read(
            "conflicting",
            0,
            [(12000, 12200), (12600, 13000), (13400, 13600)],
            [FORWARD, REVERSE],
        )
        # nothing spliced: no evidence at all
        .read("unspliced", 16, [(14000, 14500)])
        .unmapped_read("never_placed")
        .build()
    )


def run_fixer(corpus, tmp_path, **kwargs):
    out_bam = str(tmp_path / "fixed.bam")
    counters = fixer.fix_bam_strand_assignments(
        corpus.bam,
        out_bam,
        corpus.fasta,
        kwargs.pop("flip_tag", fixer.DEFAULT_FLIP_TAG),
    )
    return out_bam, counters


def records_by_name(bam_filename):
    with pysam.AlignmentFile(bam_filename, "rb", check_sq=False) as reader:
        return {read.query_name: read for read in reader}


def test_a_read_is_oriented_by_its_own_splice_motifs(corpus, tmp_path):

    out_bam, counters = run_fixer(corpus, tmp_path)
    records = records_by_name(out_bam)

    assert records["flip_to_minus"].is_reverse is True  # CT..AC beat a forward flag
    assert records["flip_to_plus"].is_reverse is False  # GT..AG beat a reverse flag
    assert records["agree_plus"].is_reverse is False
    assert records["unspliced"].is_reverse is True  # no evidence, flag stands

    assert counters["num_records_strand_flipped"] == 3
    assert counters["num_inferred_by_splice_dinucs"] == 4  # 3 flipped + agree_plus


def test_a_flipped_record_carries_the_strand_it_came_from(corpus, tmp_path):

    out_bam, _ = run_fixer(corpus, tmp_path)
    records = records_by_name(out_bam)

    assert records["flip_to_minus"].get_tag("XD") == "+"
    assert records["flip_to_plus"].get_tag("XD") == "-"

    for untouched in ("agree_plus", "noncanonical", "conflicting", "unspliced"):
        assert not records[untouched].has_tag("XD"), untouched


def test_the_flip_tag_is_selectable(corpus, tmp_path):

    out_bam, _ = run_fixer(corpus, tmp_path, flip_tag="ZS")
    records = records_by_name(out_bam)

    assert records["flip_to_minus"].get_tag("ZS") == "+"
    assert not records["flip_to_minus"].has_tag("XD")


def test_ts_moves_with_the_flag_so_the_genomic_orientation_it_encodes_holds(
    corpus, tmp_path
):
    """ts:A is read-relative: reverse + ts '+' and forward + ts '-' are the same
    genomic minus strand.  Flipping one without the other would move it."""

    out_bam, counters = run_fixer(corpus, tmp_path)
    records = records_by_name(out_bam)

    flipped = records["flip_to_plus"]
    assert flipped.is_reverse is False
    assert flipped.get_tag("ts") == "-"
    assert counters["num_ts_tags_flipped"] == 1


def test_a_secondary_alignment_is_corrected_rather_than_skipped(corpus, tmp_path):
    """Unlike the quantification intake, nothing is dropped for being secondary:
    the record has introns, so it has an orientation."""

    out_bam, _ = run_fixer(corpus, tmp_path)
    secondary = records_by_name(out_bam)["secondary_flip"]

    assert secondary.is_secondary is True
    assert secondary.is_reverse is False
    assert secondary.get_tag("XD") == "-"


def test_evidence_that_does_not_decide_leaves_the_record_alone(corpus, tmp_path):

    out_bam, counters = run_fixer(corpus, tmp_path)
    records = records_by_name(out_bam)

    assert records["noncanonical"].is_reverse is False
    assert records["conflicting"].is_reverse is False

    assert counters["num_records_unchanged_no_canonical_intron"] == 1
    assert counters["num_records_unchanged_conflicting_introns"] == 1
    assert counters["num_records_unchanged_unspliced"] == 1
    assert counters["num_records_unchanged_unmapped"] == 1
    assert counters["num_records_spliced"] == 6


def test_every_record_is_written_once_in_input_order(corpus, tmp_path):
    """A repair tool that dropped or reordered reads would be a silent library
    edit, so the output is held to the input record for record."""

    out_bam, counters = run_fixer(corpus, tmp_path)

    def identities(bam_filename):
        with pysam.AlignmentFile(bam_filename, "rb", check_sq=False) as reader:
            return [
                (read.query_name, read.reference_id, read.reference_start,
                 read.cigarstring, read.query_sequence)
                for read in reader
            ]

    assert identities(out_bam) == identities(corpus.bam)
    assert counters["num_records"] == len(corpus.reads)


def test_a_paired_record_is_left_alone_rather_than_half_flipped(tmp_path):
    """Flipping 0x10 on one mate would leave the other mate's 0x20 stale, and
    this tool never sees both mates at once."""

    built = (
        Corpus(tmp_path)
        .read("pair", 1 | 2 | 64, [(2000, 2200), (2600, 3000)], [REVERSE])
        .build()
    )

    out_bam, counters = run_fixer(built, tmp_path)
    record = records_by_name(out_bam)["pair"]

    assert record.is_reverse is False
    assert not record.has_tag("XD")
    assert counters["num_records_unchanged_paired"] == 1


def test_the_command_line_writes_an_indexed_bam(corpus, tmp_path):

    out_bam = str(tmp_path / "cli.bam")

    subprocess.check_call(
        [
            sys.executable,
            SCRIPT,
            "--input_bam",
            corpus.bam,
            "--output_bam",
            out_bam,
            "--genome",
            corpus.fasta,
        ]
    )

    assert os.path.exists(out_bam + ".bai")

    with pysam.AlignmentFile(out_bam, "rb") as reader:
        fetched = {read.query_name: read for read in reader.fetch(CONTIG)}

    assert fetched["flip_to_minus"].is_reverse is True
    assert fetched["flip_to_minus"].get_tag("XD") == "+"
