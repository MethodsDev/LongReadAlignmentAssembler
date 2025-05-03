#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
import logging
import subprocess

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="extract reads from bam file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--bam", type=str, required=True, help="input bam filename")

    parser.add_argument(
        "--read_names_file", type=str, required=True, help="file containing read names"
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output prefix: files named ${output_prefix}.${strand}.bam",
    )

    args = parser.parse_args()

    input_bam_filename = args.bam
    read_names_file = args.read_names_file
    output_bam_filename = args.output_bam

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")

    bamfile_writer = pysam.AlignmentFile(
        output_bam_filename, "wb", template=bamfile_reader
    )

    reads_want = set()
    with open(read_names_file, "rt") as fh:
        for line in fh:
            read_name = line.strip()
            reads_want.add(read_name)

    for read in bamfile_reader:
        if read.query_name in reads_want:
            bamfile_writer.write(read)

    bamfile_writer.close()

    subprocess.check_call("samtools index {}".format(output_bam_filename), shell=True)

    sys.exit(0)


if __name__ == "__main__":
    main()
