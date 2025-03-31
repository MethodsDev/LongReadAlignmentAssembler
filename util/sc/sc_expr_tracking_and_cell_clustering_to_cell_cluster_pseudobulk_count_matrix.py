#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
import argparse
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="build cell cluster pseudobulk count matrix file fron cell cluster file and full sample LRAA.quznt.expr.tracking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--cell_clusters_info",
        type=str,
        required=True,
        help="file containing cell cluster info with format:  cell_barcode(tab)cell_cluster",
    )

    parser.add_argument(
        "--LRAA_tracking_file",
        type=str,
        required=True,
        help="LRAA tracking file for the single cell data set",
    )

    parser.add_argument(
        "--output_matrix",
        type=str,
        required=True,
        help="output cluster.pseudobulk.count.matrix",
    )

    args = parser.parse_args()

    cell_clusters_info_file = args.cell_clusters_info
    LRAA_tracking_file = args.LRAA_tracking_file
    output_matrix = args.output_matrix

    # parse cell cluster info
    logger.info("-parsing {} for cell cluster info".format(cell_clusters_info_file))
    cell_clusters = set()
    cell_to_cluster = dict()
    with open(cell_clusters_info_file, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) != 2:
                continue
            cell_barcode, cluster_id = vals
            cell_to_cluster[cell_barcode] = cluster_id
            cell_clusters.add(cluster_id)

    # build matrix
    logger.info("-parsing {}, constructing matrix".format(LRAA_tracking_file))
    total_reads = 0
    read_counter = 0

    transcript_ids = set()
    transcript_id_cluster_read_count = defaultdict(lambda: defaultdict(int))

    transcript_id_to_gene_id = dict()

    with open(LRAA_tracking_file, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            read_name = row["read_name"]
            frac_read_assigned = float(row["frac_assigned"])

            total_reads += frac_read_assigned

            cell_barcode, umi, core_read_name = read_name.split("^")

            if cell_barcode in cell_to_cluster:
                cell_cluster = cell_to_cluster[cell_barcode]
                transcript_id_cluster_read_count[transcript_id][
                    cell_cluster
                ] += frac_read_assigned
                read_counter += frac_read_assigned
                transcript_ids.add(transcript_id)
                transcript_id_to_gene_id[transcript_id] = gene_id

    # write matrix
    logger.info("-writing {}".format(output_matrix))
    with open(output_matrix, "wt") as ofh:
        cluster_ids = sorted(list(cell_clusters))
        print("\t".join(["gene_id", "transcript_id"] + cluster_ids), file=ofh)
        transcript_ids = sorted(list(transcript_id_to_gene_id.keys()))
        for transcript_id in transcript_ids:
            gene_id = transcript_id_to_gene_id[transcript_id]
            read_counts = list()
            for cluster_id in cluster_ids:
                read_counts.append(
                    str(transcript_id_cluster_read_count[transcript_id][cluster_id])
                )
            print("\t".join([gene_id, transcript_id] + read_counts), file=ofh)

    logger.info("-done")

    sys.exit(0)


if __name__ == "__main__":
    main()
