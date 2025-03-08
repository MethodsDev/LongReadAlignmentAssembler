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

    parser.add_argument("--bam", type=str, required=True, help="input bam filename")

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix: files named ${output_prefix}.${cluster_name}.bam",
    )

    parser.add_argument(
        "--cell_clusters",
        type=str,
        required=True,
        help="cell cluster file, format: cell_barcode (tab) cluster_name",
    )

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix
    cell_clusters_filename = args.cell_clusters

    #########
    ### begin

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb", check_sq=False)

    # get cell cluster info
    cell_barcode_to_cluster = dict()
    with open(cell_clusters_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) != 2:
                logger.warn("Skipping line from cell clusters file: {}".format(line))
                continue
            cell_barcode, cluster_name = vals
            cell_barcode_to_cluster[cell_barcode] = cluster_name

    cell_cluster_to_ofh = dict()

    def bam_opener(cluster_name):
        bam_output_filename = output_prefix + "." + cluster_name + ".bam"
        bamfile_writer = pysam.AlignmentFile(
            bam_output_filename, "wb", template=bamfile_reader
        )
        cell_cluster_to_ofh[cluster_name] = bamfile_writer

        return bamfile_writer

    for read in bamfile_reader:
        if read.has_tag("CB"):
            cell_barcode = read.get_tag("CB")
        else:
            logger.warn("read lacks cell barcode info: " + str(read))
            continue

        if cell_barcode in cell_barcode_to_cluster:
            cluster_name = cell_barcode_to_cluster[cell_barcode]
            if cluster_name in cell_cluster_to_ofh:
                bam_writer = cell_cluster_to_ofh[cluster_name]
            else:
                bam_writer = bam_opener(cluster_name)

            bam_writer.write(read)

        else:
            logger.warn(
                "cell barcode not recognized in a cluster: {}".format(cell_barcode)
            )

    logger.info("Done.")

    for bam_writer in cell_cluster_to_ofh.values():
        bam_writer.close()

    sys.exit(0)


if __name__ == "__main__":
    main()
