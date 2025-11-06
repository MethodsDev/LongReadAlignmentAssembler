#!/usr/bin/env python3
import argparse
import re
import sys
import gzip


def parse_region(region_str):
    """
    Parse region like: chr12:6532575-6540493 or chr12+:6,532,575-6,540,493
    Returns: (chrom, lend, rend, strand or None)
    """
    m = re.match(
        r"^(?P<chr>(?:chr)?\S+?)(?P<strand>[+-])?:([\d,]+)-([\d,]+)$", region_str
    )
    if not m:
        sys.stderr.write(f"Error: invalid region format: {region_str}\n")
        sys.exit(1)

    chrom = m.group("chr")
    strand = m.group("strand")
    lend = int(m.group(3).replace(",", ""))
    rend = int(m.group(4).replace(",", ""))

    if lend >= rend:
        sys.stderr.write(f"Error: lend >= rend in region: {region_str}\n")
        sys.exit(1)

    return chrom, lend, rend, strand


def line_overlaps_region(chrom, lend, rend, strand, line, fully_contained=False):
    """
    Returns True if GTF record overlaps or is contained in region.
    """
    if line.startswith("#"):
        return False

    fields = line.strip().split("\t")
    if len(fields) < 9:
        return False

    (
        rec_chrom,
        rec_source,
        rec_feature,
        rec_start,
        rec_end,
        rec_score,
        rec_strand,
        rec_frame,
        rec_attr,
    ) = fields
    rec_start, rec_end = int(rec_start), int(rec_end)

    if rec_chrom != chrom:
        return False
    if strand and rec_strand != strand:
        return False

    if fully_contained:
        return rec_start >= lend and rec_end <= rend
    else:
        return not (rec_end < lend or rec_start > rend)


def open_maybe_gzip(fname):
    if fname.endswith(".gz"):
        return gzip.open(fname, "rt")
    return open(fname, "r")


def main():
    parser = argparse.ArgumentParser(
        description="Extract GTF records overlapping a genomic region."
    )
    parser.add_argument("--gtf", required=True, help="Input GTF file (can be gzipped)")
    parser.add_argument(
        "--region",
        required=True,
        help="Region like chr12:6,532,575-6,540,493 or chr12+:6,532,575-6,540,493",
    )
    parser.add_argument(
        "--fully-contained",
        action="store_true",
        help="Only report records fully contained within the region",
    )
    args = parser.parse_args()

    chrom, lend, rend, strand = parse_region(args.region)

    with open_maybe_gzip(args.gtf) as fh:
        for line in fh:
            if line_overlaps_region(
                chrom, lend, rend, strand, line, fully_contained=args.fully_contained
            ):
                sys.stdout.write(line)


if __name__ == "__main__":
    main()
