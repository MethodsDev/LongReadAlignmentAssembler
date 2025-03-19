#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="separate bam according to cell clustering info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam", type=str, required=True, help="input bam filename"
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output prefix: files named ${output_prefix}.${cluster_name}.bam",
    )

    parser.add_argument(
        "--cells",
        type=str,
        required=True,
        help="file containing cell barcodes to extract reads for.",
    )

    args = parser.parse_args()

    input_bam_filename = args.input_bam
    output_bam = args.output_bam
    cells_filename = args.cells

    #########
    ### begin

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb", check_sq=False)

    # get cell cluster info
    cell_barcodes_want = set()
    with open(cells_filename, "rt") as fh:
        for line in fh:
            cell_barcode = line.rstrip()
            cell_barcodes_want.add(cell_barcode)

    bam_writer = pysam.AlignmentFile(output_bam, "wb", template=bamfile_reader)

    cell_barcodes_seen = set()
    reads_kept_counter = 0
    for read in bamfile_reader:
        if read.has_tag("CB"):
            cell_barcode = read.get_tag("CB")
        else:
            logger.warn("read lacks cell barcode info: " + str(read))
            continue

        if cell_barcode in cell_barcodes_want:
            bam_writer.write(read)
            cell_barcodes_seen.add(cell_barcode)
            reads_kept_counter += 1
            print(
                "\r[{}] reads written.".format(reads_kept_counter),
                file=sys.stderr,
                end="",
            )

    logger.info("Done.")
    bam_writer.close()

    fraction_cells_found = len(cell_barcodes_seen) / len(cell_barcodes_want)

    print(
        "Found {} fraction of {} cells of interest in bam.".format(
            fraction_cells_found, len(cell_barcodes_want)
        )
    )

    sys.exit(0)


if __name__ == "__main__":
    main()
