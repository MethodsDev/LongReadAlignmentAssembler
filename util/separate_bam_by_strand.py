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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


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

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix

    infer_read_orient_flag = args.infer_read_orient
    genome_fasta = args.genome
    gtf_file = args.gtf

    #########
    ### begin

    chrom_to_itree = None

    if infer_read_orient_flag:
        if genome_fasta is None or gtf_file is None:
            sys.exit("Error - with --infer_read_orient need both --gtf and --genome")

        chrom_to_itree = build_chrom_itrees(gtf_file)

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")

    top_strand_bam_filename = output_prefix + ".+.bam"
    bottom_strand_bam_filename = output_prefix + ".-.bam"

    top_strand_bamfile_writer = pysam.AlignmentFile(
        top_strand_bam_filename, "wb", template=bamfile_reader
    )
    bottom_strand_bamfile_writer = pysam.AlignmentFile(
        bottom_strand_bam_filename, "wb", template=bamfile_reader
    )

    output_bam_files = (top_strand_bam_filename, bottom_strand_bam_filename)

    chrom_seq = None
    prev_chrom = None

    # general stats
    num_records = 0
    num_forward = 0
    num_reverse = 0
    num_neither = 0

    # infer stats
    num_inferred_by_splice_dinucs = 0
    num_inferred_by_annot_overlap = 0
    num_records_strand_flipped = 0
    num_records_strand_uncertain = 0

    for read in bamfile_reader:

        chrom = bamfile_reader.get_reference_name(read.reference_id)

        if chrom is None:
            continue
            # raise RuntimeError("Error, read has no chromosome assignment: " + str(read))

        num_records += 1

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
                num_inferred_by_splice_dinucs += 1
            else:
                # try by annotation mapping
                strand = infer_transcribed_orient_via_annotation_mapping(
                    pretty_alignment, chrom_to_itree[chrom]
                )
                if strand != "?":
                    num_inferred_by_annot_overlap += 1

            # fix strand setting in the read alignment record
            if strand != "?" and strand != init_strand:
                read.is_reverse = True if strand == "-" else False
                num_records_strand_flipped += 1

        if strand == "?":
            num_records_strand_uncertain += 1
            # set to aligned orientation
            strand = init_strand

        # write strand-specific records
        if strand == "+":
            top_strand_bamfile_writer.write(read)
            num_forward += 1

        elif strand == "-":
            bottom_strand_bamfile_writer.write(read)
            num_reverse += 1

        else:
            num_neither += 1

    top_strand_bamfile_writer.close()
    bottom_strand_bamfile_writer.close()

    # assert num_records > 0, "No records read from input bam file: {}".format(input_bam_filename)

    if num_records == 0:
        logger.warning("No aligned reads detected for {}".format(input_bam_filename))

    else:

        report_vals = [
            "Num input bam records: {}".format(num_records),
            "Num top strand: {} = {:.1f}%".format(
                num_forward, num_forward / num_records * 100
            ),
            "Num bottom strand: {} = {:.1f}%".format(
                num_reverse, num_reverse / num_records * 100
            ),
            "Num neither strand and ignored: {} = {:.1f}%".format(
                num_neither, num_neither / num_records
            ),
        ]

        if infer_read_orient_flag:
            report_vals += [
                "Num read orientations inferred by dinuc splice sites: {} = {:.1f}%".format(
                    num_inferred_by_splice_dinucs,
                    num_inferred_by_splice_dinucs / num_records * 100,
                ),
                "Num read orientations inferred by annotation overlap: {} = {:.1f}%".format(
                    num_inferred_by_annot_overlap,
                    num_inferred_by_annot_overlap / num_records * 100,
                ),
                "Num read orientations flipped: {} = {:.1f}%".format(
                    num_records_strand_flipped,
                    num_records_strand_flipped / num_records * 100,
                ),
                "Num reads with uncertaion orientation: {} = {:.1f}%".format(
                    num_records_strand_uncertain,
                    num_records_strand_uncertain / num_records * 100,
                ),
            ]

        report = "\n".join(report_vals)

        logger.info(report)

    # index the bams
    for output_bam_file in output_bam_files:
        subprocess.check_call("samtools index {}".format(output_bam_file), shell=True)

    sys.exit(0)


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
