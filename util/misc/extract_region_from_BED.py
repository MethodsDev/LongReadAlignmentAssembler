#!/usr/bin/env python3
import argparse
import re
import sys
import gzip


def parse_region(region_str):
    """
    Parse region like: chr12:6,532,575-6,540,493 or chr12+:6,532,575-6,540,493
    Returns (chrom, lend, rend, strand or None)
    """
    m = re.match(r"^(?P<chr>chr\S+?)(?P<strand>[+-])?:([\d,]+)-([\d,]+)$", region_str)
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
    Determine whether a BED record overlaps or is contained in region.
    BED: chrom start end [name score strand ...]
    Coordinates are 0-based half-open: start inclusive, end exclusive.
    """
    if line.startswith("#") or not line.strip():
        return False

    fields = line.strip().split("\t")
    if len(fields) < 3:
        return False

    rec_chrom = fields[0]
    try:
        rec_start = int(fields[1])
        rec_end = int(fields[2])
    except ValueError:
        return False

    rec_strand = fields[5] if len(fields) >= 6 else None

    if rec_chrom != chrom:
        return False
    if strand and rec_strand and rec_strand != strand:
        return False

    if fully_contained:
        return rec_start >= lend and rec_end <= rend
    else:
        return not (rec_end <= lend or rec_start >= rend)


def open_maybe_gzip(fname):
    return gzip.open(fname, "rt") if fname.endswith(".gz") else open(fname, "r")


def main():
    parser = argparse.ArgumentParser(
        description="Extract BED records overlapping a genomic region."
    )
    parser.add_argument("--bed", required=True, help="Input BED file (can be gzipped)")
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

    with open_maybe_gzip(args.bed) as fh:
        for line in fh:
            if line_overlaps_region(
                chrom, lend, rend, strand, line, args.fully_contained
            ):
                sys.stdout.write(line)


if __name__ == "__main__":
    main()
