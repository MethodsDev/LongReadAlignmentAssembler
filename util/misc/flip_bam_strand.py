#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam


def main():

    parser = argparse.ArgumentParser(
        description="flip bam read alignment strand",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam", type=str, required=True, help="input bam filename"
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output bam file with strand flipped",
    )

    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam

    #########
    ### begin

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")
    bamfile_writer = pysam.AlignmentFile(
        output_bam_filename, "wb", template=bamfile_reader
    )

    for read in bamfile_reader:

        if read.is_forward:
            read.is_forward = False
        else:
            read.is_forward = True

        bamfile_writer.write(read)

    print("Done.")

    sys.exit(0)


if __name__ == "__main__":
    main()
