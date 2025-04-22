#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict
import intervaltree as itree
import pysam
import csv
import subprocess

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Transcript import Transcript, GTF_contig_to_transcripts
from Pretty_alignment import Pretty_alignment
from SQANTI_like_annotator import SQANTI_like_annotator

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="Assign reads (bam) or isoform features (gtf) to sqanti categories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--ref_gtf",
        type=str,
        required=True,
        help="reference GTF to be used for feature comparisons and category SQANTI assignments",
    )

    parser.add_argument(
        "--input_bam",
        type=str,
        required=False,
        help="input bam with long read alignments",
    )

    parser.add_argument(
        "--input_gtf",
        type=str,
        required=False,
        help="input gtf file containing isoform structures",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix for bam and tsv files",
    )

    args = parser.parse_args()

    ref_annot_gtf = args.ref_gtf
    input_gtf = args.input_gtf
    input_bam = args.input_bam
    output_prefix = args.output_prefix

    if input_gtf is None and input_bam is None:
        exit("Error, must specify --input_gtf or --input_bam")

    if input_gtf is not None and input_bam is not None:
        exit("Error, must specify --input_gtf or --input_bam, not both together")

    sqanti_classifier = SQANTI_like_annotator(ref_annot_gtf)

    tsv_output_filename = output_prefix + ".iso_cats.tsv"
    tsv_ofh = open(tsv_output_filename, "wt")
    tsv_writer = csv.DictWriter(
        tsv_ofh,
        fieldnames=[
            "feature_name",
            "sqanti_cat",
            "feature_length",
            "num_exon_segments",
            "structure",
            "matching_isoforms",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    tsv_writer.writeheader()

    feature_counter = 0
    feature_category_counter = defaultdict(int)

    if input_bam is not None:

        ## Examine aligned reads

        logger.info("Classifying reads from bam: {}".format(input_bam))
        bamfile_reader = pysam.AlignmentFile(input_bam, "rb")

        bam_output_filename = output_prefix + ".iso_cats.bam"
        bamwriter = pysam.AlignmentFile(
            bam_output_filename, "wb", template=bamfile_reader
        )

        for read in bamfile_reader:

            feature_counter += 1
            if feature_counter % 1000 == 0:
                print("\r[{}]  ".format(feature_counter), file=sys.stderr, end="")

            if (
                read.is_mapped
                and (not read.is_secondary)
                and (not read.is_supplementary)
            ):

                read_class_info = classify_read(read, bamfile_reader, sqanti_classifier)
                read_class_info["feature_length"] = len(read.query_sequence)

                read.set_tag("CL", read_class_info["sqanti_cat"], "Z")
                read.set_tag("CI", read_class_info["matching_isoforms"], "Z")

                tsv_writer.writerow(read_class_info)
                feature_category_counter[read_class_info["sqanti_cat"]] += 1

            bamwriter.write(read)

        bamwriter.close()

    else:
        ## Examine isoforms
        contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
            input_gtf
        )
        for contig, transcript_list in contig_to_transcripts.items():
            for transcript_obj in transcript_list:
                transcript_strand = transcript_obj.get_strand()

                feature_counter += 1
                if feature_counter % 1000 == 0:
                    print("\r[{}]  ".format(feature_counter), file=sys.stderr, end="")

                transcript_id = transcript_obj.get_transcript_id()
                read_class_info = sqanti_classifier.classify_alignment_or_isoform(
                    contig, transcript_strand, transcript_id, transcript_obj
                )
                read_class_info["feature_length"] = transcript_obj.get_cdna_len()

                tsv_writer.writerow(read_class_info)
                feature_category_counter[read_class_info["sqanti_cat"]] += 1

    tsv_ofh.close()

    # write summary counts
    summary_counts_tsv = output_prefix + ".iso_cats.summary_counts.tsv"
    with open(summary_counts_tsv, "wt") as ofh:
        print("\t".join(["Category", "Count"]), file=ofh)
        for feature_category, count in feature_category_counter.items():
            print("\t".join([feature_category, str(count)]), file=ofh)

    # make barplot of cat counts.
    summary_counts_plot_name = output_prefix + ".iso_cats.summary_counts.pdf"
    cmd = " ".join(
        [
            os.path.join(os.path.dirname(__file__), "misc/plot_SQANTI_cats.Rscript"),
            summary_counts_tsv,
            summary_counts_plot_name,
        ]
    )
    subprocess.check_call(cmd, shell=True)

    logger.info("\nDone. See files: {}.*".format(output_prefix))

    sys.exit(0)


def classify_read(read, bamfile_reader, sqanti_classifier):

    chrom = bamfile_reader.get_reference_name(read.reference_id)

    read_name = read.query_name
    read_strand = "+" if read.is_forward else "-"

    stranded_chrom = "{}:{}".format(chrom, read_strand)

    pretty_alignment = Pretty_alignment.get_pretty_alignment(read)

    read_class_info = sqanti_classifier.classify_alignment_or_isoform(
        chrom, read_strand, read_name, pretty_alignment
    )

    return read_class_info


if __name__ == "__main__":
    main()
