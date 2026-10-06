#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict
import intervaltree as itree
import pysam
import csv
import gzip
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Transcript import Transcript, GTF_contig_to_transcripts
from Pretty_alignment import Pretty_alignment
import Util_funcs
from SQANTI_like_annotator import SQANTI_like_annotator

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)

BAM_TSV_FIELDNAMES = [
    "feature_name",
    "sqanti_cat",
    "read_length",
    "alignment_length",
    "num_exon_segments",
    "structure",
    "matching_isoforms",
]


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

    parser.add_argument(
        "--CPU",
        type=int,
        default=1,
        help="with --input_bam: number of contigs to classify in parallel (needs a .bai)",
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

    # The summary barplot is drawn by an Rscript at the very end. Check for R now
    # rather than after classifying every read: lraa-core ships without R, and a
    # run there used to fail only once all the work was done.
    if shutil.which("Rscript") is None:
        exit(
            "Error, Rscript not found on PATH; it is needed for the summary plot. "
            "Run in the lraa-sc image, which includes R."
        )

    if args.CPU < 1:
        exit("Error, --CPU must be at least 1")

    if input_bam is not None and args.CPU > 1:
        feature_category_counter = classify_bam_by_contig(
            ref_annot_gtf, input_bam, output_prefix, args.CPU
        )
        write_summary_and_plot(output_prefix, feature_category_counter)
        sys.exit(0)

    sqanti_classifier = SQANTI_like_annotator(ref_annot_gtf)

    tsv_output_filename = output_prefix + ".iso_cats.tsv"
    tsv_ofh = open(tsv_output_filename, "wt")

    feature_counter = 0
    feature_category_counter = defaultdict(int)

    if input_bam is not None:
        ## Examine aligned reads
        tsv_writer = csv.DictWriter(
            tsv_ofh,
            fieldnames=BAM_TSV_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        tsv_writer.writeheader()

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

            process_bam_record(
                read,
                bamfile_reader,
                sqanti_classifier,
                tsv_writer,
                bamwriter,
                feature_category_counter,
            )

        bamwriter.close()

    else:
        ## Examine isoforms

        tsv_writer = csv.DictWriter(
            tsv_ofh,
            fieldnames=[
                "feature_name",
                "sqanti_cat",
                "cDNA_length",
                "num_exon_segments",
                "structure",
                "matching_isoforms",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        tsv_writer.writeheader()

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
                read_class_info["cDNA_length"] = transcript_obj.get_cdna_len()

                tsv_writer.writerow(read_class_info)
                feature_category_counter[read_class_info["sqanti_cat"]] += 1

    tsv_ofh.close()

    write_summary_and_plot(output_prefix, feature_category_counter)

    sys.exit(0)


def write_summary_and_plot(output_prefix, feature_category_counter):

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


def get_aligned_length(read):
    aligned_length = 0
    for operation, length in read.cigartuples:
        if operation in [0, 7, 8]:  # M, =, or X
            aligned_length += length
    return aligned_length


def process_bam_record(
    read, bamfile_reader, sqanti_classifier, tsv_writer, bamwriter, category_counter
):
    # Shared by the serial loop and the per-contig workers, so the two paths
    # classify, tag and count identically. Every record is written to the bam;
    # only mapped primary alignments are classified.
    if read.is_mapped and (not read.is_secondary) and (not read.is_supplementary):

        read_class_info = classify_read(read, bamfile_reader, sqanti_classifier)
        read_class_info["read_length"] = len(read.query_sequence)
        read_class_info["alignment_length"] = get_aligned_length(read)

        read.set_tag("CL", read_class_info["sqanti_cat"], "Z")
        read.set_tag("CI", read_class_info["matching_isoforms"], "Z")

        tsv_writer.writerow(read_class_info)
        category_counter[read_class_info["sqanti_cat"]] += 1

    bamwriter.write(read)


def classify_bam_by_contig(ref_annot_gtf, input_bam, output_prefix, num_workers):
    """Classify a coordinate-sorted, indexed bam one contig per worker.

    Each worker loads only its contig's reference transcripts, so per-worker memory
    is small and the full reference is never parsed in one process (a whole-genome
    GENCODE annotation takes ~1 min and ~1.7 GB to load). Classification depends only
    on the read's own contig: the interval trees, splice sites and splice patterns are
    all keyed by contig. Parts are concatenated in bam header order, unplaced
    unmapped reads last, which reproduces the serial run's output exactly.
    """

    bam = pysam.AlignmentFile(input_bam, "rb")
    if not bam.has_index():
        exit(
            "Error, --CPU > 1 needs an index for {} (samtools index it first)".format(
                input_bam
            )
        )

    # contigs carrying any record (placed unmapped reads count), in header order
    read_counts = {
        stat.contig: stat.mapped + stat.unmapped for stat in bam.get_index_statistics()
    }
    contigs = [c for c in bam.references if read_counts.get(c, 0) > 0]
    num_unplaced = bam.nocoordinate
    bam.close()

    tmpdir = output_prefix + ".__sqanti_parts"
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    contig_gtfs = split_gtf_by_contig(ref_annot_gtf, contigs, tmpdir)

    # one job per contig, plus one for the unplaced reads (contig None)
    jobs = []
    for i, contig in enumerate(contigs):
        part = os.path.join(tmpdir, "part{:05d}".format(i))
        jobs.append((contig, contig_gtfs[contig], part))
    if num_unplaced > 0:
        jobs.append((None, None, os.path.join(tmpdir, "part_unplaced")))

    logger.info(
        "Classifying {} contigs ({} unplaced reads) from {} with {} workers".format(
            len(contigs), num_unplaced, input_bam, num_workers
        )
    )

    # largest first, so the long contigs are not left until the end
    submit_order = sorted(
        range(len(jobs)),
        key=lambda j: -read_counts.get(jobs[j][0], num_unplaced),
    )
    part_counters = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=min(num_workers, len(jobs))) as executor:
        futures = {
            executor.submit(classify_bam_contig, input_bam, *jobs[j]): j
            for j in submit_order
        }
        for num_done, future in enumerate(as_completed(futures), start=1):
            j = futures[future]
            part_counters[j] = future.result()
            logger.info(
                "[{}/{}] done: {}".format(num_done, len(jobs), jobs[j][0] or "unplaced")
            )

    # merge in header order. Summing the counters in that order keeps the
    # categories in first-seen order, as the serial run writes them.
    feature_category_counter = defaultdict(int)
    for part_counter in part_counters:
        for category, count in part_counter.items():
            feature_category_counter[category] += count

    tsv_output_filename = output_prefix + ".iso_cats.tsv"
    with open(tsv_output_filename, "wt") as ofh:
        csv.DictWriter(
            ofh, fieldnames=BAM_TSV_FIELDNAMES, delimiter="\t", lineterminator="\n"
        ).writeheader()
        for _, _, part in jobs:
            with open(part + ".tsv", "rt") as fh:
                shutil.copyfileobj(fh, ofh)

    bam_output_filename = output_prefix + ".iso_cats.bam"
    pysam.cat(
        "--no-PG", "-o", bam_output_filename, *[part + ".bam" for _, _, part in jobs]
    )

    shutil.rmtree(tmpdir)

    return feature_category_counter


def split_gtf_by_contig(gtf_filename, contigs, tmpdir):
    # One pass over the reference, writing each wanted contig's lines to its own
    # file. Contigs with reads but no annotation get an empty file.
    contig_gtfs = {
        contig: os.path.join(tmpdir, "ref{:05d}.gtf".format(i))
        for i, contig in enumerate(contigs)
    }
    handles = {contig: open(path, "wt") for contig, path in contig_gtfs.items()}
    opener = gzip.open if gtf_filename.endswith(".gz") else open
    with opener(gtf_filename, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            ofh = handles.get(line.split("\t", 1)[0])
            if ofh is not None:
                ofh.write(line)
    for ofh in handles.values():
        ofh.close()
    return contig_gtfs


def classify_bam_contig(input_bam, contig, contig_gtf, part_prefix):
    # Worker: takes only paths and strings, so it runs under any start method.
    if contig_gtf is not None:
        sqanti_classifier = SQANTI_like_annotator(contig_gtf)
    else:
        sqanti_classifier = None  # unplaced reads are unmapped; nothing to classify

    bamfile_reader = pysam.AlignmentFile(input_bam, "rb")
    bamwriter = pysam.AlignmentFile(part_prefix + ".bam", "wb", template=bamfile_reader)
    category_counter = defaultdict(int)

    with open(part_prefix + ".tsv", "wt") as tsv_ofh:
        tsv_writer = csv.DictWriter(
            tsv_ofh,
            fieldnames=BAM_TSV_FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        records = (
            bamfile_reader.fetch(contig)
            if contig is not None
            else bamfile_reader.fetch("*")
        )
        for read in records:
            process_bam_record(
                read,
                bamfile_reader,
                sqanti_classifier,
                tsv_writer,
                bamwriter,
                category_counter,
            )

    bamwriter.close()
    bamfile_reader.close()

    return dict(category_counter)


def classify_read(read, bamfile_reader, sqanti_classifier):

    chrom = bamfile_reader.get_reference_name(read.reference_id)

    read_name = read.query_name
    # TRANSCRIBED strand (ts, fallback flag), so an antisense-sequenced cDNA read
    # is categorized against its transcript's strand, not how it aligned.
    read_strand = Util_funcs.transcribed_strand(read)

    stranded_chrom = "{}:{}".format(chrom, read_strand)

    pretty_alignment = Pretty_alignment.get_pretty_alignment(read)

    read_class_info = sqanti_classifier.classify_alignment_or_isoform(
        chrom, read_strand, read_name, pretty_alignment
    )

    return read_class_info


if __name__ == "__main__":
    main()
