#!/usr/bin/env python3

"""What fix_bam_strand_assignments.py promises about the bam it writes.

Two claims are tested, and they pull in opposite directions.  The tool must
CHANGE the strand of every record whose evidence contradicts its flag, and it
must change NOTHING ELSE: not the record count, not the order, not a
coordinate, not a record it had no evidence about.  A fixture where no read
flips would pass the second claim while testing nothing, so the flips are
asserted by name.

The annotation fallback is held to its rank as well as its arithmetic: a read
whose own splice motifs decide it keeps that verdict even when the annotation
underneath says otherwise.
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
    """A contig whose splice sites and annotation are chosen, not incidental."""

    def __init__(self, tmp_path):
        self.dir = tmp_path
        self.sequence = list(("ACGT" * (LENGTH // 4 + 1))[:LENGTH])
        self.reads = []  # (name, flag, blocks, tags)
        self.transcripts = []  # (transcript_id, strand, exons)

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

    def transcript(self, transcript_id, strand, exons):
        self.transcripts.append((transcript_id, strand, exons))
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

        self.gtf = str(self.dir / "annot.gtf")
        with open(self.gtf, "wt") as ofh:
            for transcript_id, strand, exons in self.transcripts:
                attrs = 'gene_id "{}"; transcript_id "{}";'.format(
                    transcript_id, transcript_id
                )
                for lend, rend in exons:
                    print(
                        "\t".join(
                            (
                                CONTIG,
                                "test",
                                "exon",
                                str(lend),
                                str(rend),
                                ".",
                                strand,
                                ".",
                                attrs,
                            )
                        ),
                        file=ofh,
                    )

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
        # aligner and splice sites agree, and the annotation under it disagrees
        .read("agree_plus", 0, [(6000, 6200), (6600, 7000)], [FORWARD])
        # a secondary alignment is evidence about itself, and is corrected too
        .read("secondary_flip", 16 | 256, [(8000, 8200), (8600, 9000)], [FORWARD])
        # spliced, but the motifs are not canonical in either orientation
        .read("noncanonical", 0, [(10000, 10200), (10600, 11000)], [NONCANONICAL])
        # one intron votes each way: no majority from the read itself
        .read(
            "conflicting",
            0,
            [(12000, 12200), (12600, 13000), (13400, 13600)],
            [FORWARD, REVERSE],
        )
        # nothing spliced and nothing annotated: no evidence at all
        .read("unspliced", 16, [(14000, 14500)])
        # nothing spliced, but it sits in an annotated minus strand exon
        .read("annot_to_minus", 0, [(16100, 16400)])
        # the plus transcript has more exons here, the minus one more bases
        .read("annot_best_by_bases", 0, [(18000, 18300)])
        # annotated on both strands, equally: still nothing to conclude
        .read("annot_tie", 0, [(20000, 20100)])
        .unmapped_read("never_placed")
        # contradicts agree_plus, and must lose to it
        .transcript("shadowM", "-", [(6000, 6200), (6600, 7000)])
        .transcript("geneM", "-", [(16000, 16500)])
        .transcript("manyExonP", "+", [(18000, 18010), (18100, 18110), (18200, 18210)])
        .transcript("oneExonM", "-", [(18000, 18300)])
        .transcript("tieP", "+", [(20000, 20100)])
        .transcript("tieM", "-", [(20000, 20100)])
        .build()
    )


def run_fixer(corpus, tmp_path, gtf=False, flip_tag=fixer.DEFAULT_FLIP_TAG):
    out_bam = str(tmp_path / ("fixed.annot.bam" if gtf else "fixed.bam"))
    contig_to_exon_itree = (
        fixer.build_contig_exon_itrees(corpus.gtf) if gtf else None
    )
    counters = fixer.fix_bam_strand_assignments(
        corpus.bam,
        out_bam,
        corpus.fasta,
        contig_to_exon_itree,
        flip_tag,
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

    assert counters["num_splice_undecided_no_canonical_intron"] == 1
    assert counters["num_splice_undecided_conflicting_introns"] == 1
    # unspliced, and the three reads the annotation would speak for
    assert counters["num_splice_undecided_unspliced"] == 4
    assert counters["num_records_unchanged_unmapped"] == 1
    assert counters["num_records_spliced"] == 6
    assert counters["num_inferred_by_annot_overlap"] == 0  # no gtf was supplied


def test_without_a_gtf_an_unannotated_read_is_simply_uncertain(corpus, tmp_path):

    _, counters = run_fixer(corpus, tmp_path)

    # noncanonical, conflicting, unspliced, annot_to_minus, annot_best_by_bases,
    # annot_tie
    assert counters["num_records_strand_uncertain"] == 6


# ------------------------------------------------------- the annotation fallback


def test_annotation_orients_a_read_its_own_cigar_cannot(corpus, tmp_path):

    out_bam, counters = run_fixer(corpus, tmp_path, gtf=True)
    records = records_by_name(out_bam)

    annotated = records["annot_to_minus"]
    assert annotated.is_reverse is True  # aligned forward, geneM is on minus
    assert annotated.get_tag("XD") == "+"
    assert counters["num_inferred_by_annot_overlap"] == 2
    assert counters["num_records_strand_flipped"] == 5  # 3 by motif, 2 by annotation


def test_annotation_never_overrules_the_reads_own_splice_motifs(corpus, tmp_path):
    """shadowM covers agree_plus exactly and is on the other strand.  The read's
    GT..AG says forward and is direct evidence, so the annotation is not asked."""

    out_bam, counters = run_fixer(corpus, tmp_path, gtf=True)
    records = records_by_name(out_bam)

    assert records["agree_plus"].is_reverse is False
    assert not records["agree_plus"].has_tag("XD")
    assert counters["num_inferred_by_splice_dinucs"] == 4


def test_the_best_transcript_is_the_one_sharing_the_most_bases(corpus, tmp_path):
    """manyExonP puts three exons under this read and oneExonM puts one, but the
    one covers 301 bases against 33.  Counting exon records would pick plus."""

    out_bam, _ = run_fixer(corpus, tmp_path, gtf=True)
    record = records_by_name(out_bam)["annot_best_by_bases"]

    assert record.is_reverse is True
    assert record.get_tag("XD") == "+"


def test_annotation_on_both_strands_alike_decides_nothing(corpus, tmp_path):

    out_bam, counters = run_fixer(corpus, tmp_path, gtf=True)
    record = records_by_name(out_bam)["annot_tie"]

    assert record.is_reverse is False
    assert not record.has_tag("XD")
    # noncanonical, conflicting, unspliced, annot_tie
    assert counters["num_records_strand_uncertain"] == 4


def test_a_read_off_every_annotated_contig_is_not_an_error(corpus, tmp_path):
    """The gtf need not describe the whole bam; an unannotated contig is simply
    a contig the fallback has nothing to say about."""

    itrees = fixer.build_contig_exon_itrees(corpus.gtf)
    assert fixer.infer_orient_via_annotation([(100, 200)], itrees.get("chrOther")) == "?"


@pytest.mark.parametrize("with_gtf", (False, True))
def test_every_record_is_accounted_for_exactly_once(corpus, tmp_path, with_gtf):
    """The outcome counters partition the file, so a record cannot be silently
    counted twice or not at all."""

    _, counters = run_fixer(corpus, tmp_path, gtf=with_gtf)

    accounted = sum(
        counters["num_records_{}".format(outcome)] for outcome in fixer.OUTCOMES
    )
    assert accounted == counters["num_records"] == len(corpus.reads)

    assert (
        counters["num_inferred_by_splice_dinucs"]
        + counters["num_inferred_by_annot_overlap"]
        == counters["num_records_strand_agree"]
        + counters["num_records_strand_flipped"]
    )


def test_every_record_is_written_once_in_input_order(corpus, tmp_path):
    """A repair tool that dropped or reordered reads would be a silent library
    edit, so the output is held to the input record for record."""

    out_bam, counters = run_fixer(corpus, tmp_path, gtf=True)

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


def test_a_gzipped_gtf_is_read_as_annotation(corpus, tmp_path):

    import gzip

    gzipped = str(tmp_path / "annot.gtf.gz")
    with open(corpus.gtf, "rt") as fh, gzip.open(gzipped, "wt") as ofh:
        ofh.write(fh.read())

    assert fixer.build_contig_exon_itrees(gzipped) == fixer.build_contig_exon_itrees(
        corpus.gtf
    )


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
            "--gtf",
            corpus.gtf,
        ]
    )

    assert os.path.exists(out_bam + ".bai")

    with pysam.AlignmentFile(out_bam, "rb") as reader:
        fetched = {read.query_name: read for read in reader.fetch(CONTIG)}

    assert fetched["flip_to_minus"].is_reverse is True
    assert fetched["flip_to_minus"].get_tag("XD") == "+"
    assert fetched["annot_to_minus"].is_reverse is True
