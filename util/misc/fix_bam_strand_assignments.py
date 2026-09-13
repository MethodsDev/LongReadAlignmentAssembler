#!/usr/bin/env python3

"""Rewrite read alignment strands to the orientation the splice motifs imply.

The bam counterpart to fix_gtf_strand_assignments.py.  An aligner sets a read's
reverse bit from which genomic strand the read SEQUENCE matched, which for an
unstranded cDNA library says nothing about which strand was TRANSCRIBED.  The
introns do say it: a canonical donor/acceptor pair reads GT..AG on the
transcribed strand and CT..AC on the other one, so a spliced read carries its
own transcribed orientation regardless of how the aligner happened to place it.

Every input record is written exactly once, in input order.  A record is left
alone unless its own introns vote for the orientation opposite its flag; then
bit 0x10 is set to the transcribed orientation and the ORIGINAL aligned strand
is recorded in a tag (XD by default), so the change is both visible and
reversible.

What a flip does and does not mean:

  - SEQ, QUAL, CIGAR and POS are untouched, so the alignment stays valid
    against the genome and the file stays coordinate sorted.  Only the claim
    the flag makes about the source molecule changes.
  - After this tool, 0x10 means "the transcript is on the minus strand", not
    "SEQ is the reverse complement of the read as sequenced".  Anything that
    reverse complements SEQ by the flag to recover the original read must use
    the XD tag instead.
  - minimap2's ts:A is read-relative, so flipping the flag without flipping ts
    would change the genomic orientation ts encodes.  ts is flipped alongside,
    leaving that orientation where it was.
  - SA:Z strand fields describe the aligner's mapping and are left as found.

Reads with no intron, no canonical intron, or an even split between canonical
orientations carry no splice evidence and are passed through unchanged.
"""

import sys, os
import argparse
import logging
from collections import defaultdict

import pysam

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import Util_funcs
from GenomeFeature import Intron
from Pretty_alignment import Pretty_alignment


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


BAM_CREF_SKIP = 3  # the cigar N operation, the only source of an intron

DEFAULT_FLIP_TAG = "XD"

# Counted for every record, so the tallies below add up to num_records.
PASS_THROUGH_REASONS = (
    "unmapped",
    "paired",
    "unspliced",
    "no_canonical_intron",
    "conflicting_introns",
)


def main():

    parser = argparse.ArgumentParser(
        description="fix bam read alignment strands using intron splice motifs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam", type=str, required=True, help="input bam filename"
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output bam file, strands set to the transcribed orientation",
    )

    parser.add_argument(
        "--genome", type=str, required=True, help="genome fasta file (faidx indexed)"
    )

    parser.add_argument(
        "--flip_tag",
        type=str,
        default=DEFAULT_FLIP_TAG,
        help="tag set on flipped records, holding the original aligned strand",
    )

    args = parser.parse_args()

    if len(args.flip_tag) != 2:
        sys.exit("Error - --flip_tag must be a two character bam tag")

    counters = fix_bam_strand_assignments(
        args.input_bam,
        args.output_bam,
        args.genome,
        args.flip_tag,
    )

    report_counters(counters)

    index_if_coordinate_sorted(args.output_bam)

    logger.info("Done.")

    sys.exit(0)


def new_counters():

    counters = defaultdict(int)

    counters["num_records"] = 0
    counters["num_records_spliced"] = 0
    counters["num_inferred_by_splice_dinucs"] = 0
    counters["num_records_strand_agree"] = 0
    counters["num_records_strand_flipped"] = 0
    counters["num_ts_tags_flipped"] = 0
    for reason in PASS_THROUGH_REASONS:
        counters["num_records_unchanged_{}".format(reason)] = 0

    return counters


def fix_bam_strand_assignments(
    input_bam_filename,
    output_bam_filename,
    genome_fasta,
    flip_tag=DEFAULT_FLIP_TAG,
):
    """Stream every record, flipping the strand bit where the introns disagree.

    Records are handled one at a time and written immediately, so memory does
    not grow with the bam.  Bam_alignment_extractor is deliberately not the
    intake here: it DISCARDS records (unmapped, secondary, supplementary, low
    identity), materializes a list, and lightens away the pysam record a
    rewrite has to write back out.  A bam repair tool that dropped reads would
    quietly shrink the library, so the pysam iterator is the intake and
    Pretty_alignment is built per record for its intron structure alone.
    """

    counters = new_counters()

    contig_seq = None
    contig_seq_acc = None
    contigs_loaded = set()
    warned_paired = False

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")
    bamfile_writer = pysam.AlignmentFile(
        output_bam_filename, "wb", template=bamfile_reader
    )

    for read in bamfile_reader:

        counters["num_records"] += 1

        reason = pass_through_reason(read)
        if reason is not None:
            if reason == "paired" and not warned_paired:
                warned_paired = True
                logger.warning(
                    "paired records found (eg. {}) and left untouched: flipping one "
                    "mate's 0x10 would strand the other mate's 0x20".format(
                        read.query_name
                    )
                )
            counters["num_records_unchanged_{}".format(reason)] += 1
            bamfile_writer.write(read)
            continue

        contig_acc = read.reference_name

        if contig_acc != contig_seq_acc:
            if contig_acc in contigs_loaded:
                logger.warning(
                    "contig {} revisited: input is not grouped by contig, so the "
                    "genome sequence is being re-read".format(contig_acc)
                )
            contig_seq = Util_funcs.retrieve_contig_seq_from_fasta_file(
                contig_acc, genome_fasta
            )
            contig_seq_acc = contig_acc
            contigs_loaded.add(contig_acc)

        pretty_alignment = Pretty_alignment.get_pretty_alignment(read)
        introns = pretty_alignment.get_introns()

        if not introns:
            # an N shorter than read_aln_gap_merge_int merged back into an exon
            counters["num_records_unchanged_unspliced"] += 1
            bamfile_writer.write(read)
            continue

        counters["num_records_spliced"] += 1

        num_top, num_bottom = count_canonical_intron_orients(introns, contig_seq)

        if num_top == num_bottom:
            reason = (
                "conflicting_introns" if num_top > 0 else "no_canonical_intron"
            )
            counters["num_records_unchanged_{}".format(reason)] += 1
            bamfile_writer.write(read)
            continue

        transcribed_orient = "+" if num_top > num_bottom else "-"
        counters["num_inferred_by_splice_dinucs"] += 1

        aligned_orient = "-" if read.is_reverse else "+"

        if transcribed_orient == aligned_orient:
            counters["num_records_strand_agree"] += 1
        else:
            flip_record_strand(read, aligned_orient, flip_tag, counters)
            counters["num_records_strand_flipped"] += 1

        bamfile_writer.write(read)

    bamfile_writer.close()
    bamfile_reader.close()

    return counters


def pass_through_reason(read):
    """Why this record carries no usable splice evidence, or None if it may.

    Unlike the quantification filters, nothing here is about read quality:
    secondary, supplementary, duplicate, qcfail and low identity records all
    have introns of their own and are corrected like any other.  Only records
    whose strand cannot be decided from their own cigar, or cannot be flipped
    without breaking something else, are passed through.
    """

    if read.is_unmapped or read.reference_id < 0:
        return "unmapped"

    if read.is_paired:
        # flipping 0x10 here would leave the mate's 0x20 describing the old strand
        return "paired"

    cigartuples = read.cigartuples
    if not cigartuples:
        return "unspliced"

    for opcode, _ in cigartuples:
        if opcode == BAM_CREF_SKIP:
            return None

    # no N: deletions merge back into their neighbors, so there is no intron to read
    return "unspliced"


def count_canonical_intron_orients(introns, contig_seq):
    """Votes for each transcribed orientation among this read's introns.

    Intron.check_canonical_splicing is the same classifier the splice graph
    admits junctions with, so a read is oriented by the motifs LRAA already
    believes, rather than by a fourth private copy of the dinucleotide sets.
    """

    num_top_strand = 0
    num_bottom_strand = 0

    for intron_lend, intron_rend in introns:
        orient = Intron.check_canonical_splicing(intron_lend, intron_rend, contig_seq)

        if orient == "+":
            num_top_strand += 1
        elif orient == "-":
            num_bottom_strand += 1

    return num_top_strand, num_bottom_strand


def flip_record_strand(read, orig_aligned_orient, flip_tag, counters):
    """Set 0x10 to the transcribed orientation, keeping the record readable back.

    The original strand goes into the flip tag rather than a bare flag, because
    the flip is only reversible if the strand it replaced survives, and a
    consumer asking "was this read moved" and "where was it" should not need
    two tags.
    """

    read.is_reverse = orig_aligned_orient == "+"

    read.set_tag(flip_tag, orig_aligned_orient, value_type="A")

    if read.has_tag("ts"):
        ts = read.get_tag("ts")
        if ts in ("+", "-"):
            # ts is relative to the read, so it moves with the read's strand
            read.set_tag("ts", "-" if ts == "+" else "+", value_type="A")
            counters["num_ts_tags_flipped"] += 1

    return


def report_counters(counters):

    num_records = counters["num_records"]

    def as_pct(count):
        if num_records < 1:
            return 0.0
        return 100.0 * count / num_records

    logger.info("-records examined: {}".format(num_records))
    for counter_name in sorted(counters.keys()):
        if counter_name == "num_records":
            continue
        count = counters[counter_name]
        logger.info(
            "-{}: {} ({:.2f}%)".format(counter_name, count, as_pct(count))
        )

    return


def index_if_coordinate_sorted(bam_filename):
    """Index the output when the input ordering survived, which it always does.

    A strand flip changes no coordinate, so a coordinate sorted input yields a
    coordinate sorted output.  An input in any other order cannot be indexed
    and is left as written.
    """

    with pysam.AlignmentFile(bam_filename, "rb", check_sq=False) as reader:
        sort_order = reader.header.to_dict().get("HD", {}).get("SO")

    if sort_order != "coordinate":
        logger.info(
            "-output not indexed: header sort order is {}".format(sort_order)
        )
        return

    pysam.index(bam_filename)
    logger.info("-output indexed.")

    return


if __name__ == "__main__":
    main()
