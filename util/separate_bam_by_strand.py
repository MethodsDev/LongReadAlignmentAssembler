#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
import logging
import subprocess
import intervaltree as itree
import gzip
from collections import defaultdict

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)
from Pretty_alignment import Pretty_alignment
import Util_funcs
import LRAA_Globals

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# The intron-length rule lives in pylib/Util_funcs.py, which is the single
# implementation shared with pylib/Bam_alignment_extractor.py so that the
# splice-graph input and read assignment always see one record set.  Re-bound here
# because it is part of this module's filtering vocabulary.
get_longest_intron_length = Util_funcs.get_longest_intron_length

# reasons for discarding an input record, in the order they are tested.
# Each discarded record is counted under exactly one reason, so that the
# discard counts plus the written records sum to the records read.
DISCARD_REASONS = (
    "unmapped",
    "no_chromosome",
    "secondary",
    "supplementary",
    "duplicate",
    "qcfail",
    "long_intron",
)


def discard_counter_name(discard_reason):
    return "num_records_discarded_" + discard_reason


def get_discard_reason(read, max_intron_length):
    """returns the reason for discarding the read, or None if the read is retained.

    max_intron_length <= 0 disables intron length filtering.
    """

    if read.is_unmapped:
        return "unmapped"

    if read.reference_id < 0:
        return "no_chromosome"

    if read.is_secondary:
        return "secondary"

    if read.is_supplementary:
        return "supplementary"

    if read.is_duplicate:
        return "duplicate"

    if read.is_qcfail:
        return "qcfail"

    if Util_funcs.has_disqualifying_long_intron(read, max_intron_length):
        return "long_intron"

    return None


def main():

    parser = argparse.ArgumentParser(
        description="separate bam into strand-specific bam files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--bam", type=str, required=True, help="input bam filename")

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix: files named ${output_prefix}.${strand}.bam",
    )

    parser.add_argument(
        "--infer_read_orient",
        action="store_true",
        default=False,
        help="infer read orientation based on splicing evidence or overlapping reference annotations (in that order)",
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=False,
        help="reference annotation used for inferring transcribed orientation of read (required if --infer_read_orient",
    )

    parser.add_argument(
        "--genome",
        type=str,
        required=False,
        help="genome fasta file, required if --infer_read_orient",
    )

    parser.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
        help="alignments containing any intron (CIGAR 'N' operation) longer than this are discarded; set to 0 or a negative value to disable intron length filtering",
    )

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix

    infer_read_orient_flag = args.infer_read_orient
    genome_fasta = args.genome
    gtf_file = args.gtf
    max_intron_length = args.max_intron_length

    #########
    ### begin

    chrom_to_itree = None

    if infer_read_orient_flag:
        if genome_fasta is None:
            sys.exit("Error - with --infer_read_orient need --genome")

        if gtf_file:
            chrom_to_itree = build_chrom_itrees(gtf_file)

    top_strand_bam_filename = output_prefix + ".+.bam"
    bottom_strand_bam_filename = output_prefix + ".-.bam"

    output_bam_files = (top_strand_bam_filename, bottom_strand_bam_filename)

    if max_intron_length > 0:
        logger.info(
            "-discarding alignments having any intron longer than {}".format(
                max_intron_length
            )
        )
    else:
        logger.info("-intron length filtering disabled")

    counters = split_bam_by_strand(
        input_bam_filename,
        top_strand_bam_filename,
        bottom_strand_bam_filename,
        max_intron_length,
        infer_read_orient_flag=infer_read_orient_flag,
        genome_fasta=genome_fasta,
        chrom_to_itree=chrom_to_itree,
    )

    if counters["num_records"] == 0:
        logger.warning("No records detected for {}".format(input_bam_filename))
    else:
        logger.info(build_report(counters, infer_read_orient_flag))

    # index the bams
    for output_bam_file in output_bam_files:
        subprocess.check_call("samtools index {}".format(output_bam_file), shell=True)

    sys.exit(0)


def split_bam_by_strand(
    input_bam_filename,
    top_strand_bam_filename,
    bottom_strand_bam_filename,
    max_intron_length,
    infer_read_orient_flag=False,
    genome_fasta=None,
    chrom_to_itree=None,
):
    """writes the retained records of the input bam to the strand-specific bam files
    and returns the record accounting as a dict of counters.

    Only primary, non-supplementary alignments passing the intron length
    criterion are retained; every other input record is discarded here and so
    excluded from all downstream processing.
    """

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")

    top_strand_bamfile_writer = pysam.AlignmentFile(
        top_strand_bam_filename, "wb", template=bamfile_reader
    )
    bottom_strand_bamfile_writer = pysam.AlignmentFile(
        bottom_strand_bam_filename, "wb", template=bamfile_reader
    )

    chrom_seq = None
    prev_chrom = None

    counters = {
        # general stats
        "num_records": 0,
        "num_forward": 0,
        "num_reverse": 0,
        "num_neither": 0,
        # infer stats
        "num_inferred_by_splice_dinucs": 0,
        "num_inferred_by_annot_overlap": 0,
        "num_records_strand_flipped": 0,
        "num_records_strand_uncertain": 0,
    }

    # discard stats, one counter per discard reason
    for discard_reason in DISCARD_REASONS:
        counters[discard_counter_name(discard_reason)] = 0

    for read in bamfile_reader:

        counters["num_records"] += 1

        discard_reason = get_discard_reason(read, max_intron_length)
        if discard_reason is not None:
            counters[discard_counter_name(discard_reason)] += 1
            continue

        chrom = bamfile_reader.get_reference_name(read.reference_id)

        strand = "+" if read.is_forward else "-"
        init_strand = strand

        if infer_read_orient_flag:

            if prev_chrom is None or prev_chrom != chrom:
                prev_chrom = chrom
                chrom_seq = Util_funcs.retrieve_contig_seq_from_fasta_file(
                    chrom, genome_fasta
                )

            pretty_alignment = Pretty_alignment.get_pretty_alignment(read)

            # first try by dinuc splice sites of spliced introns from read
            strand = infer_spliced_orient(pretty_alignment, chrom_seq)
            if strand != "?":
                counters["num_inferred_by_splice_dinucs"] += 1
            elif chrom_to_itree is not None:
                # try by annotation mapping
                strand = infer_transcribed_orient_via_annotation_mapping(
                    pretty_alignment, chrom_to_itree[chrom]
                )
                if strand != "?":
                    counters["num_inferred_by_annot_overlap"] += 1

            # fix strand setting in the read alignment record
            if strand != "?" and strand != init_strand:
                read.is_reverse = True if strand == "-" else False
                counters["num_records_strand_flipped"] += 1

        if strand == "?":
            counters["num_records_strand_uncertain"] += 1
            # set to aligned orientation
            strand = init_strand

        # write strand-specific records
        if strand == "+":
            top_strand_bamfile_writer.write(read)
            counters["num_forward"] += 1

        elif strand == "-":
            bottom_strand_bamfile_writer.write(read)
            counters["num_reverse"] += 1

        else:
            counters["num_neither"] += 1

    top_strand_bamfile_writer.close()
    bottom_strand_bamfile_writer.close()
    bamfile_reader.close()

    return counters


def build_report(counters, infer_read_orient_flag=False):
    """record accounting summary, every count reported against the number of
    records read as denominator"""

    num_records = counters["num_records"]

    def as_pct_of_total(count):
        return "{} = {:.1f}% of {}".format(
            count, count / num_records * 100, num_records
        )

    report_vals = ["Num input bam records: {}".format(num_records)]

    num_records_discarded = 0
    for discard_reason in DISCARD_REASONS:
        count = counters[discard_counter_name(discard_reason)]
        num_records_discarded += count
        report_vals.append(
            "Num records discarded as {}: {}".format(
                discard_reason, as_pct_of_total(count)
            )
        )

    report_vals += [
        "Num records discarded (all reasons): {}".format(
            as_pct_of_total(num_records_discarded)
        ),
        "Num records retained: {}".format(
            as_pct_of_total(num_records - num_records_discarded)
        ),
        "Num top strand: {}".format(as_pct_of_total(counters["num_forward"])),
        "Num bottom strand: {}".format(as_pct_of_total(counters["num_reverse"])),
        "Num neither strand and ignored: {}".format(
            as_pct_of_total(counters["num_neither"])
        ),
    ]

    if infer_read_orient_flag:
        report_vals += [
            "Num read orientations inferred by dinuc splice sites: {}".format(
                as_pct_of_total(counters["num_inferred_by_splice_dinucs"])
            ),
            "Num read orientations inferred by annotation overlap: {}".format(
                as_pct_of_total(counters["num_inferred_by_annot_overlap"])
            ),
            "Num read orientations flipped: {}".format(
                as_pct_of_total(counters["num_records_strand_flipped"])
            ),
            "Num reads with uncertain orientation: {}".format(
                as_pct_of_total(counters["num_records_strand_uncertain"])
            ),
        ]

    return "\n".join(report_vals)


def infer_spliced_orient(pretty_alignment, contig_seq):

    introns_coordsets = pretty_alignment.get_introns()

    if len(introns_coordsets) < 1:
        return "?"

    return majority_vote_intron_orient(introns_coordsets, contig_seq)


def majority_vote_intron_orient(intron_coordsets, contig_seq):

    splice_dinucs_top_strand = {"GTAG", "GCAG", "ATAC"}
    splice_dinucs_bottom_strand = {
        "CTAC",
        "CTGC",
        "GTAT",
    }  # revcomp of top strand dinucs

    orient_counts = {"+": 0, "-": 0}

    for intron_coordset in intron_coordsets:
        intron_lend, intron_rend = intron_coordset
        dinuc_left = contig_seq[intron_lend - 1] + contig_seq[intron_lend - 1 + 1]
        dinuc_right = contig_seq[intron_rend - 1 - 1] + contig_seq[intron_rend - 1]
        dinuc_combo = dinuc_left + dinuc_right

        if dinuc_combo in splice_dinucs_top_strand:
            orient_counts["+"] += 1
        elif dinuc_combo in splice_dinucs_bottom_strand:
            orient_counts["-"] += 1

    # check tie or not match
    if orient_counts["+"] == orient_counts["-"]:
        return "?"
    elif orient_counts["+"] > orient_counts["-"]:
        return "+"
    else:
        return "-"


def build_chrom_itrees(gtf_file):

    logger.info("-building chrom itrees from: " + gtf_file)

    chrom_to_itree = defaultdict(lambda: itree.IntervalTree())

    if re.search(gtf_file, "\\.gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue

            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 8:
                continue

            chrom = vals[0]
            lend = int(vals[3])
            rend = int(vals[4])
            feature_type = vals[2]

            if feature_type != "exon":
                continue

            strand = vals[6]

            chrom_to_itree[chrom][lend : rend + 1] = strand

    return chrom_to_itree


def infer_transcribed_orient_via_annotation_mapping(pretty_alignment, chrom_itree):

    alignment_segments = pretty_alignment.get_pretty_alignment_segments()

    strand_counter = defaultdict(int)

    for segment in alignment_segments:
        lend, rend = segment

        for overlapping_exon in chrom_itree[lend : rend + 1]:
            strand = overlapping_exon.data
            strand_counter[strand] += 1

    if strand_counter["+"] == strand_counter["-"]:
        return "?"
    elif strand_counter["+"] > strand_counter["-"]:
        return "+"
    else:
        return "-"


if __name__ == "__main__":
    main()
