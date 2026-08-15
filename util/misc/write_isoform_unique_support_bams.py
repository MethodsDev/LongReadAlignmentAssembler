#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict
import csv
import gzip
import logging
import argparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def open_maybe_gzip(path, mode="rt"):
    """Open a tracking file whether or not it is gzipped.

    LRAA writes quant.tracking gzipped under --gzip_tracking, and the tracking file is the
    one artifact that scales with library size, so at scale it is the compressed form that
    exists. Detected by suffix, matching every other tracking reader in this repo
    (annotate_bam_with_read_tracking_info.py, extract_tracked_reads_from_bam.py,
    reassign_multigene_tracking_reads.py, singlecell_tracking_to_sparse_matrix.py).
    """
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")

def iter_non_comment_lines(fh):
    for line in fh:
        if line.startswith("#"):
            continue
        yield line


def main():
    parser = argparse.ArgumentParser(
        description="extract reads providing unique isoform support from bam",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="input bam file",
    )

    parser.add_argument(
        "--tracking", type=str, required=True, help="LRAA read tracking file"
    )

    args = parser.parse_args()

    bam_file = args.bam
    tracking_file = args.tracking

    logger.info(
        "-extracting unique read to transcript assignments from {}".format(
            tracking_file
        )
    )
    unique_read_names_to_transcripts = dict()
    # get unique read names to transcript ids
    with open_maybe_gzip(tracking_file, "rt") as fh:
        reader = csv.DictReader(iter_non_comment_lines(fh), delimiter="\t")
        for row in reader:
            if float(row["frac_assigned"]) == 1.0:
                read_name = row["read_name"]
                transcript_id = row["transcript_id"]
                unique_read_names_to_transcripts[read_name] = transcript_id

    logger.info(
        "-extracting the read alignments for those uniquely assigned to transcripts from: {}".format(
            bam_file
        )
    )
    # get the read alignments for the unique reads.
    bamreader = pysam.AlignmentFile(bam_file, "rb")

    transcripts_to_unique_read_alignments = defaultdict(list)
    for read in bamreader.fetch():
        read_name = read.query_name
        if read_name in unique_read_names_to_transcripts:
            transcript_id = unique_read_names_to_transcripts[read_name]
            transcripts_to_unique_read_alignments[transcript_id].append(read)

    # write transcript-specific unique read bam files.
    logger.info("-writing transcript-specific unique bam read alignment files.")
    for transcript_id, bam_reads in transcripts_to_unique_read_alignments.items():
        bam_output_filename = "{}.unique_reads.bam".format(transcript_id)
        logger.info("-writing {}".format(bam_output_filename))
        bamwriter = pysam.AlignmentFile(bam_output_filename, "wb", template=bamreader)
        for bam_read in bam_reads:
            bamwriter.write(bam_read)
        bamwriter.close()

    logger.info("Done")

    sys.exit(0)


if __name__ == "__main__":
    main()
