#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extract reads from a BAM file that are uniquely assigned to a single isoform
(based on the `XI:Z:` tag) and retain only the N longest reads per isoform.

Example usage:
    ./extract_unique_isoform_reads.py \
        --bam input.bam \
        --output unique_longest.bam \
        --topN 10

Author: Brian Haas (via ChatGPT GPT-5)
"""

import pysam
import argparse
import sys
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract uniquely assigned reads (single XI:Z isoform) "
        "and keep top N longest reads per isoform."
    )
    parser.add_argument(
        "--bam",
        required=True,
        help="Input BAM file (must be indexed if random access is needed)",
    )
    parser.add_argument(
        "--output", required=True, help="Output BAM file for filtered reads"
    )
    parser.add_argument(
        "--topN",
        type=int,
        default=10,
        help="Number of longest reads to retain per isoform (default: 10)",
    )
    parser.add_argument(
        "--min_mapq", type=int, default=0, help="Optional MAPQ filter (default: 0)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        bam_in = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        sys.stderr.write(f"[ERROR] Failed to open input BAM: {e}\n")
        sys.exit(1)

    bam_out = pysam.AlignmentFile(args.output, "wb", header=bam_in.header)

    # Collect reads by isoform
    isoform_reads = defaultdict(list)

    for read in bam_in:
        if read.is_unmapped or read.mapping_quality < args.min_mapq:
            continue

        try:
            xi_tag = read.get_tag("XI")
        except KeyError:
            continue  # skip reads without XI tag

        isoforms = xi_tag.split(",")
        if len(isoforms) != 1:
            continue  # skip multi-mapped reads

        isoform = isoforms[0]
        read_length = read.query_length or 0

        isoform_reads[isoform].append((read_length, read.to_string()))

    bam_in.close()

    # Filter and write top N reads per isoform
    total_written = 0
    for isoform, reads in isoform_reads.items():
        reads.sort(key=lambda x: x[0], reverse=True)
        for _, aln_str in reads[: args.topN]:
            aln = pysam.AlignedSegment.fromstring(aln_str, bam_out.header)
            bam_out.write(aln)
            total_written += 1

    bam_out.close()

    sys.stderr.write(f"[INFO] Wrote {total_written} reads to {args.output}\n")


if __name__ == "__main__":
    main()
