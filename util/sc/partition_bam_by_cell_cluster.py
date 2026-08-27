#!/usr/bin/env python3

import sys, os, re
import itertools
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

    parser.add_argument(
        "--restrict_to_chromosomes",
        type=str,
        default=None,
        help="restrict the EMITTED cluster bams to these references, named as in "
        "the bam header and separated by spaces or commas. WHOLE contigs only, "
        "which is what makes it sound: depth-normalization downstream anchors its "
        "window grid per contig on that contig's first aligned base, so dropping "
        "entire references leaves every retained reference's grid -- and therefore "
        "which reads it thins -- untouched, while a sub-contig filter would move "
        "that anchor. Reads on other references are read and counted, then not "
        "written. Unset means emit every reference, the prior behaviour",
    )

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix
    cell_clusters_filename = args.cell_clusters
    cell_barcode_tag = args.cell_barcode_tag
    num_threads = args.threads

    # Named references to emit, or None for all of them. Accepts either separator
    # because main_chromosomes is space-separated in the WDLs while ChunkedRun's
    # --contigs is comma-separated, and this is called from both conventions.
    restrict_to_chromosomes = None
    if args.restrict_to_chromosomes:
        restrict_to_chromosomes = set(
            re.split(r"[,\s]+", args.restrict_to_chromosomes.strip())
        )
        restrict_to_chromosomes.discard("")

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
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) != 2:
                logger.warn("Skipping line from cell clusters file: {}".format(line))
                continue
            cell_barcode, cluster_name = vals
            # The Seurat step writes a HEADER ROW, and it has exactly two
            # tab-separated fields, so the length check above admits it: the run
            # then gets a cluster literally named "cluster" holding one barcode
            # named "cell_barcode", which matches no read. The partitioner
            # materializes an empty bam for every cluster by design, so that
            # phantom reaches a per-cluster LRAA run as a header-only bam and
            # dies on `count_reads_from_bam`'s "no reads counted" assertion.
            # Observed as LRAA_by_cluster-27 / test_sc_wdl.cluster.bam in the
            # sc_full_pipe cluster-guided run.
            #
            # Matched on the EXACT first line rather than dropped positionally:
            # this argument is documented as headerless and the fixed-cluster
            # fixtures are, so an unconditional skip would discard a real
            # assignment.
            if lineno == 1 and (cell_barcode, cluster_name) == (
                "cell_barcode",
                "cluster",
            ):
                logger.info(
                    "Skipping cell clusters header line: {}".format(line)
                )
                continue
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
    reads_outside_restriction = 0
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

    # Materialize an EMPTY bam for every cluster before dispatching any read, then
    # let the loop below open writers lazily exactly as it always has.
    #
    # The file has to exist even for a cluster that ends up with nothing, because
    # the task globs `*.bam` and `*.cells.txt` into two arrays that are paired BY
    # INDEX downstream: a missing bam shifts every later cluster's cell list onto
    # another cluster's bam, silently. Whole-genome that was unreachable in
    # practice -- every cluster has reads somewhere. With --restrict_to_chromosomes
    # it is ordinary, since a small cluster can legitimately hold nothing on a
    # handful of contigs.
    #
    # Opened and CLOSED here rather than kept live: holding one writer per cluster
    # open across the whole scan would cost a file descriptor and num_threads
    # compression threads per cluster from the outset. The lazy path below reopens
    # (and truncates) this placeholder for any cluster that does have reads, so the
    # resource profile of the scan is unchanged and only genuinely empty clusters
    # keep a header-only bam.
    for cluster_name in sorted(cluster_to_cell_barcodes.keys()):
        placeholder_filename = output_prefix + "." + cluster_name + ".bam"
        pysam.AlignmentFile(
            placeholder_filename, "wb", template=bamfile_reader
        ).close()
    logger.info(
        f"Materialized {len(cluster_to_cell_barcodes)} empty cluster bams; "
        "writers open lazily as reads arrive"
    )

    # WHICH RECORDS GET VISITED, as opposed to which get written.
    #
    # Restricting only the write leaves this scan reading the whole library to use
    # a fraction of it, which is most of what restricting was supposed to avoid. So
    # when an index is available the references are FETCHED BY NAME and nothing else
    # is ever decompressed. Whole references in header order, so the output stays
    # coordinate-sorted exactly as a full scan would leave it.
    #
    # Without an index there is no way to skip records, so the full scan stands and
    # the write-side filter below is what makes it correct. That filter is kept on
    # both paths deliberately: it costs one set-membership test per visited read and
    # it means the two paths cannot disagree about what lands in a cluster bam.
    read_source = bamfile_reader
    if restrict_to_chromosomes is not None:
        # Refuse a name the bam cannot hold rather than silently emitting nothing
        # for it -- the same policy ChunkedRun applies to --contigs.
        absent = sorted(restrict_to_chromosomes - set(bamfile_reader.references))
        if absent:
            logger.error(
                "--restrict_to_chromosomes names {} absent from the bam header: {}".format(
                    len(absent), " ".join(absent)
                )
            )
            sys.exit(1)

        fetch_order = [
            ref for ref in bamfile_reader.references if ref in restrict_to_chromosomes
        ]
        if bamfile_reader.has_index():
            logger.info(
                "indexed fetch: reading only {} of {} references ({})".format(
                    len(fetch_order), len(bamfile_reader.references), " ".join(fetch_order)
                )
            )
            read_source = itertools.chain.from_iterable(
                bamfile_reader.fetch(reference=ref) for ref in fetch_order
            )
        else:
            logger.warning(
                "bam has no index, so every record must be visited and the "
                "restriction can only be applied on write. Index the bam to read "
                "just the %d requested references.",
                len(fetch_order),
            )

    for read in read_source:
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

            # Tested BEFORE a writer is acquired, so a cluster whose reads all sit
            # on unretained references never opens one and keeps the header-only
            # placeholder written above. The read is counted as processed and
            # deliberately not written: it belongs to a real cell of a real
            # cluster, it just sits on a reference this run does not analyze.
            if (
                restrict_to_chromosomes is not None
                and read.reference_name not in restrict_to_chromosomes
            ):
                reads_outside_restriction += 1
                continue

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
    if restrict_to_chromosomes is not None:
        logger.info(
            f"  Reads dropped by --restrict_to_chromosomes: "
            f"{reads_outside_restriction:,} "
            f"(emitting {len(restrict_to_chromosomes)} references)"
        )
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
