#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess
import collections
import pysam


def main():

    parser = argparse.ArgumentParser(
        description="extract genome contig and bam for target contig of interest",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--genome_fa", type=str, required=True, help="genome fasta file"
    )

    parser.add_argument("--bam", type=str, required=True, help="bam aligned reads")

    parser.add_argument("--gtf", type=str, required=False, help="gtf annotation")

    parser.add_argument(
        "--region",
        type=str,
        required=True,
        help="contig region of interest. Formatted chr\\d[+-]?:\\d+-\\d+ ",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=False,
        default=None,
        help="prefix for output files. Default: region string.",
    )

    parser.add_argument(
        "--max_left_extension",
        type=int,
        default=20000,
        help="how far the left edge may widen to keep reads that start upstream of the region. "
        "The overhang is bounded by alignment span, not read length, so a single spliced read can "
        "reach megabases; anything beyond this is dropped and reported.",
    )

    parser.add_argument(
        "--strict_window",
        action="store_true",
        default=False,
        help="extract exactly the requested region, discarding reads that start before it. "
        "Default is to widen the left edge so those reads are kept: a read overhanging the "
        "start cannot be rebased into the extracted contig, and silently dropping it removes "
        "evidence the region genuinely has.",
    )

    args = parser.parse_args()

    genome_fa_filename = args.genome_fa
    bam_filename = args.bam
    gtf_filename = args.gtf
    region = args.region
    output_prefix = args.output_prefix
    max_left_extension = args.max_left_extension

    region = region.replace(",", "")

    if output_prefix is None:
        output_prefix = region

    m = re.match("([^\\:\\+\\-]+)([\\+\\-]?):(\\d+)-(\\d+)", region)
    if m is None:
        raise RuntimeError("Cannot parse region string: {}".format(region))

    chrom = m.group(1)
    strand = m.group(2)
    lend = int(m.group(3))
    rend = int(m.group(4))

    # Reads starting left of the region cannot be rebased into it -- their
    # alignment would begin before position 1. Dropping them silently removes
    # evidence the region genuinely has: at one locus 22% of the reads spanning
    # the target began 5.3 kb upstream, and without them the models built from
    # them could not form, so the extraction gave an answer the whole genome does
    # not. Widen the left edge to cover them instead.
    #
    # Bounded, because the overhang is limited by ALIGNMENT SPAN and not by read
    # length: a spliced read can reach across megabases, and following the worst
    # one has been observed to turn a 100 kb request into 3.1 Mb and a few
    # thousand reads into 242,000. Reads overhanging beyond the bound are still
    # dropped, but the loss is now reported rather than silent.
    counts = collections.Counter()
    furthest = lend
    furthest_recoverable = lend
    with pysam.AlignmentFile(bam_filename, "rb") as probe:
        for read in probe.fetch(chrom, lend - 1, rend):
            start = read.reference_start + 1
            if start < lend:
                counts["overhanging"] += 1
                furthest = min(furthest, start)
                if start >= lend - max_left_extension:
                    counts["recoverable"] += 1
                    furthest_recoverable = min(furthest_recoverable, start)

    if counts["overhanging"]:
        if args.strict_window:
            print(
                f"WARNING: {counts['overhanging']} read(s) overlapping {chrom}:{lend}-{rend} start up to "
                f"{lend - furthest} bp upstream and are DISCARDED under --strict_window; the extracted "
                f"region under-represents its own left edge.",
                file=sys.stderr,
            )
        else:
            widened = lend - furthest_recoverable
            dropped = counts["overhanging"] - counts["recoverable"]
            if widened:
                print(
                    f"NOTE: widening left edge by {widened} bp ({lend} -> {lend - widened}) to keep "
                    f"{counts['recoverable']} of {counts['overhanging']} read(s) starting upstream of the "
                    f"requested region.",
                    file=sys.stderr,
                )
            if dropped:
                print(
                    f"WARNING: {dropped} read(s) start further upstream than --max_left_extension "
                    f"({max_left_extension} bp), as far as {lend - furthest} bp, and are DISCARDED. "
                    f"Raise --max_left_extension to keep them, at the cost of a larger region.",
                    file=sys.stderr,
                )
            lend -= widened

    # extract contig
    genome_region_fa_filename = f"{output_prefix}.fa"
    cmd = f"samtools faidx {genome_fa_filename} {chrom}:{lend}-{rend}"
    fa_seq = subprocess.check_output(cmd, shell=True).decode()
    # remove header line
    fa_seq = "\n".join([f">{chrom}"] + fa_seq.split("\n")[1:])
    with open(genome_region_fa_filename, "wt") as ofh:
        print(fa_seq, file=ofh)

    cmd = f"samtools faidx {genome_region_fa_filename}"
    subprocess.check_call(cmd, shell=True)

    # extract bam
    bamreader = pysam.AlignmentFile(bam_filename, "rb")

    output_bam_filename = output_prefix + ".bam"
    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)

    for read in bamreader.fetch(chrom, lend - 1, rend):
        if strand != "":
            if (strand == "+" and not read.is_forward) or (
                strand == "-" and not read.is_reverse
            ):
                continue

        # reference_start is 0-based; a read beginning at the 1-based region
        # start has reference_start == lend - 1 and belongs in the extraction
        if read.reference_start >= lend - 1:
            read.reference_start -= lend - 1
            bamwriter.write(read)

    bamwriter.close()

    cmd = f"samtools index {output_bam_filename}"
    subprocess.check_call(cmd, shell=True)

    # extract gtf
    if gtf_filename is not None:
        gtf_output_filename = output_prefix + ".gtf"
        with open(gtf_filename, "rt") as fh, open(gtf_output_filename, "wt") as ofh:

            for line in fh:
                vals = line.split("\t")
                if len(vals) < 8:
                    continue
                gtf_chrom, gtf_lend, gtf_rend, gtf_strand = (
                    vals[0],
                    int(vals[3]),
                    int(vals[4]),
                    vals[6],
                )
                if gtf_chrom != chrom:
                    continue
                if strand != "" and strand != gtf_strand:
                    continue

                if gtf_lend >= lend and gtf_rend <= rend:
                    gtf_lend -= lend - 1
                    gtf_rend -= lend - 1
                    vals[3] = str(gtf_lend)
                    vals[4] = str(gtf_rend)
                    print("\t".join(vals), file=ofh, end="")

    sys.exit(0)


if __name__ == "__main__":
    main()
