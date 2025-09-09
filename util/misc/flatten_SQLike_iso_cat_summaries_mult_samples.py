#!/usr/bin/env python3

import argparse
import gzip
import os
import csv
from collections import defaultdict


def parse_file(filename):
    """
    Parse a single gzipped tsv file into a dict {category: count}
    """
    basename = os.path.basename(filename)
    sample_id = basename.split(".iso_cats.summary_counts.tsv.gz")[0]

    cat_counts = {}
    with gzip.open(filename, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cat = row["Category"]
            count = int(row["Count"])
            cat_counts[cat] = count

    return sample_id, cat_counts


def main():
    parser = argparse.ArgumentParser(
        description="Convert iso_cats summaries into wide format (counts + fractions)"
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        required=True,
        help="Input .iso_cats.summary_counts.tsv.gz files",
    )
    parser.add_argument(
        "-c", "--counts_out", required=True, help="Output TSV file for counts"
    )
    parser.add_argument(
        "-f", "--fractions_out", required=True, help="Output TSV file for fractions"
    )
    args = parser.parse_args()

    # Collect data
    all_samples = {}
    all_categories = set()

    for filename in args.inputs:
        sample_id, cat_counts = parse_file(filename)
        all_samples[sample_id] = cat_counts
        all_categories.update(cat_counts.keys())

    # Sort categories for consistent column order
    categories = sorted(all_categories)

    # Write counts table
    with open(args.counts_out, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["sample_id"] + categories)
        for sample_id, cat_counts in sorted(all_samples.items()):
            row = [sample_id] + [cat_counts.get(cat, 0) for cat in categories]
            writer.writerow(row)

    # Write fractions table
    with open(args.fractions_out, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["sample_id"] + categories)
        for sample_id, cat_counts in sorted(all_samples.items()):
            total = sum(cat_counts.values())
            row = [sample_id] + [
                (cat_counts.get(cat, 0) / total if total > 0 else 0.0)
                for cat in categories
            ]
            writer.writerow(row)


if __name__ == "__main__":
    main()
