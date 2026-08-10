#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
from collections import defaultdict
import logging
import psutil


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

    parser.add_argument(
        "--cell_barcode_tag",
        type=str,
        default="CB",
        help="BAM tag containing the cell barcode",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="number of threads for BAM compression/decompression",
    )

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix
    cell_clusters_filename = args.cell_clusters
    cell_barcode_tag = args.cell_barcode_tag
    num_threads = args.threads

    #########
    ### begin

    logger.info(f"Starting BAM partitioning by cell cluster")
    logger.info(f"Input BAM: {input_bam_filename}")
    logger.info(f"Output prefix: {output_prefix}")
    logger.info(f"Cell clusters file: {cell_clusters_filename}")
    logger.info(f"Cell barcode tag: {cell_barcode_tag}")
    logger.info(f"Threads: {num_threads}")

    # Initialize memory tracking
    process = psutil.Process()
    initial_memory = process.memory_info().rss / 1024 / 1024  # MB
    logger.info(f"Initial memory usage: {initial_memory:.2f} MB")

    logger.info("Opening input BAM file...")
    bamfile_reader = pysam.AlignmentFile(
        input_bam_filename, "rb", check_sq=False, threads=num_threads
    )

    # get cell cluster info
    logger.info("Loading cell cluster assignments...")
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

    logger.info(f"Loaded {len(cell_barcode_to_cluster)} cell barcodes across clusters")

    # Emit each cluster's called-cell barcodes alongside its BAM. Downstream isoform
    # identification needs the roster of real cells for that cluster: barcodes are
    # what distinguish called cells from empty droplets and ambient signal, and the
    # BAM alone cannot convey cells that contributed no read.
    cluster_to_cell_barcodes = defaultdict(list)
    for cell_barcode, cluster_name in cell_barcode_to_cluster.items():
        cluster_to_cell_barcodes[cluster_name].append(cell_barcode)
    for cluster_name, cell_barcodes in cluster_to_cell_barcodes.items():
        cell_list_filename = f"{output_prefix}.{cluster_name}.cells.txt"
        with open(cell_list_filename, "wt") as cell_list_fh:
            for cell_barcode in sorted(cell_barcodes):
                print(cell_barcode, file=cell_list_fh)
        logger.info(
            f"Wrote {len(cell_barcodes)} cell barcodes for cluster {cluster_name} "
            f"to {cell_list_filename}"
        )
    memory_mb = process.memory_info().rss / 1024 / 1024
    logger.info(f"Memory after loading clusters: {memory_mb:.2f} MB")

    cell_cluster_to_ofh = dict()

    # Initialize read counters
    total_reads_processed = 0
    reads_without_barcode = 0
    reads_unrecognized_cluster = 0
    reads_written = 0
    report_interval = 100000  # Report progress every 100k reads

    logger.info("Starting to process reads...")

    def bam_opener(cluster_name):
        bam_output_filename = output_prefix + "." + cluster_name + ".bam"
        logger.info(f"Creating output BAM for cluster: {cluster_name}")
        bamfile_writer = pysam.AlignmentFile(
            bam_output_filename, "wb", template=bamfile_reader, threads=num_threads
        )
        cell_cluster_to_ofh[cluster_name] = bamfile_writer
        logger.info(f"Total clusters opened so far: {len(cell_cluster_to_ofh)}")

        return bamfile_writer

    for read in bamfile_reader:
        total_reads_processed += 1

        # Progress reporting
        if total_reads_processed % report_interval == 0:
            memory_mb = process.memory_info().rss / 1024 / 1024
            logger.info(
                f"Processed {total_reads_processed:,} reads | "
                f"Written: {reads_written:,} | "
                f"No {cell_barcode_tag}: {reads_without_barcode:,} | "
                f"Unrecognized: {reads_unrecognized_cluster:,} | "
                f"Clusters: {len(cell_cluster_to_ofh)} | "
                f"Memory: {memory_mb:.2f} MB"
            )

        if read.has_tag(cell_barcode_tag):
            cell_barcode = read.get_tag(cell_barcode_tag)
        else:
            reads_without_barcode += 1
            if reads_without_barcode <= 10:  # Only warn for first 10
                logger.warning(
                    f"read lacks cell barcode tag {cell_barcode_tag}: {read}"
                )
            continue

        if cell_barcode in cell_barcode_to_cluster:
            cluster_name = cell_barcode_to_cluster[cell_barcode]
            if cluster_name in cell_cluster_to_ofh:
                bam_writer = cell_cluster_to_ofh[cluster_name]
            else:
                bam_writer = bam_opener(cluster_name)

            bam_writer.write(read)
            reads_written += 1

        else:
            reads_unrecognized_cluster += 1
            if reads_unrecognized_cluster <= 10:  # Only warn for first 10
                logger.warn(
                    "cell barcode not recognized in a cluster: {}".format(cell_barcode)
                )

    # Final summary
    final_memory = process.memory_info().rss / 1024 / 1024
    logger.info("=" * 80)
    logger.info("Processing complete - Summary:")
    logger.info(f"  Total reads processed: {total_reads_processed:,}")
    logger.info(f"  Reads written to clusters: {reads_written:,}")
    logger.info(
        f"  Reads without cell barcode tag {cell_barcode_tag}: {reads_without_barcode:,}"
    )
    logger.info(f"  Reads with unrecognized cluster: {reads_unrecognized_cluster:,}")
    logger.info(f"  Total clusters created: {len(cell_cluster_to_ofh)}")
    logger.info(f"  Peak memory usage: {final_memory:.2f} MB")
    logger.info(f"  Memory increase: {final_memory - initial_memory:.2f} MB")
    logger.info("=" * 80)

    logger.info("Closing output BAM files...")
    for cluster_name, bam_writer in cell_cluster_to_ofh.items():
        bam_writer.close()
    logger.info(f"Closed {len(cell_cluster_to_ofh)} output BAM files")

    logger.info("Done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
