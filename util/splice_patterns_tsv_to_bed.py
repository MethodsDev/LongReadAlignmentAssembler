#!/usr/bin/env python3
"""
Convert introns data to BED12 format representing exons.
Parses intron coordinates, infers exon positions, and creates BED12 entries 
with introns_hashcode as feature names. Adds 25bp terminal exons.
"""

# initial version by claude.io, modified by bhaas

import pandas as pd
import re
import argparse
from collections import defaultdict


def parse_introns(introns_str):
    """
    Parse introns string to extract chromosome, strand, and coordinate pairs.

    Example input: "chr10:(+)[(100913251, 100915998), (100916043, 100916569), ...]"
    Returns: (chromosome, strand, [(start1, end1), (start2, end2), ...])
    """
    # Extract chromosome and strand
    match = re.match(r"(chr\w+):\(([+-])\)\[(.*)\]", introns_str)
    if not match:
        raise ValueError(f"Cannot parse introns string: {introns_str}")

    chromosome = match.group(1)
    strand = match.group(2)
    coords_str = match.group(3)

    # Extract coordinate pairs
    coord_pairs = []
    for coord_match in re.finditer(r"\((\d+),\s*(\d+)\)", coords_str):
        start = int(coord_match.group(1))
        end = int(coord_match.group(2))
        coord_pairs.append((start, end))

    return chromosome, strand, coord_pairs


def introns_to_exons_bed12(introns_str, name, terminal_exon_size=25):
    """
    Convert introns to BED12 format representing exons.
    Introns represent gaps between exons, so we infer exon positions
    and add terminal exons of specified size.

    Args:
        introns_str: String containing intron coordinates
        name: Feature name (introns_hashcode)
        terminal_exon_size: Size of terminal exons to add (default: 25bp)

    Returns:
        Dictionary with BED12 fields
    """
    chromosome, strand, intron_pairs = parse_introns(introns_str)

    if not intron_pairs:
        return None

    # Sort introns by start position
    intron_pairs.sort()

    # Infer exon coordinates from intron gaps
    exon_pairs = []

    # Add first terminal exon (25bp before first intron)
    first_intron_start = intron_pairs[0][0]
    first_exon_start = first_intron_start - terminal_exon_size
    first_exon_end = first_intron_start
    exon_pairs.append((first_exon_start, first_exon_end))

    # Add exons between introns
    for i in range(len(intron_pairs) - 1):
        # Exon starts where current intron ends, ends where next intron starts
        exon_start = intron_pairs[i][1]
        exon_end = intron_pairs[i + 1][0]
        exon_pairs.append((exon_start, exon_end))

    # Add last terminal exon (25bp after last intron)
    last_intron_end = intron_pairs[-1][1]
    last_exon_start = last_intron_end
    last_exon_end = last_intron_end + terminal_exon_size
    exon_pairs.append((last_exon_start, last_exon_end))

    # Calculate overall start and end
    chrom_start = exon_pairs[0][0]
    chrom_end = exon_pairs[-1][1]

    # Calculate block sizes and starts
    block_count = len(exon_pairs)
    block_sizes = []
    block_starts = []

    for start, end in exon_pairs:
        block_sizes.append(end - start)
        block_starts.append(start - chrom_start)

    # Create BED12 entry
    bed_entry = {
        "chrom": chromosome,
        "chromStart": chrom_start,
        "chromEnd": chrom_end,
        "name": name,
        "score": 0,  # Default score
        "strand": strand,
        "thickStart": chrom_start,  # Same as chromStart for now
        "thickEnd": chrom_end,  # Same as chromEnd for now
        "itemRgb": "0,0,0",  # Black color
        "blockCount": block_count,
        "blockSizes": ",".join(map(str, block_sizes)),
        "blockStarts": ",".join(map(str, block_starts)),
    }

    return bed_entry


def convert_to_bed12(input_file, output_file, terminal_exon_size=25):
    """
    Convert the input dataframe to BED12 format representing exons.

    Args:
        input_file: Path to input TSV file (requires column headers: introns, introns_hashcode)
        output_file: Path to output BED12 file
        terminal_exon_size: Size of terminal exons to add (default: 25bp)
    """
    # Read the input file
    df = pd.read_csv(input_file, sep="\t")

    # Create a dictionary to store unique introns -> hashcode mappings
    unique_introns = {}

    # Process each row to collect unique introns
    for _, row in df.iterrows():
        introns = row["introns"]
        hashcode = row["introns_hashcode"]

        # Skip rows with missing data
        if pd.isna(introns) or pd.isna(hashcode):
            continue

        # Store the mapping (introns string maps to hashcode)
        if introns not in unique_introns:
            unique_introns[introns] = hashcode

    print(f"Found {len(unique_introns)} unique intron structures")

    # Convert to BED12 format
    bed_entries = []
    for introns_str, hashcode in unique_introns.items():
        try:
            bed_entry = introns_to_exons_bed12(
                introns_str, hashcode, terminal_exon_size
            )
            if bed_entry:
                bed_entries.append(bed_entry)
        except Exception as e:
            print(f"Error processing introns {hashcode}: {e}")
            continue

    # Create DataFrame and sort by chromosome and position
    bed_df = pd.DataFrame(bed_entries)

    # Sort by chromosome and start position
    # First, extract chromosome number for proper sorting
    def chrom_sort_key(chrom):
        if chrom.startswith("chr"):
            chrom_num = chrom[3:]
            if chrom_num.isdigit():
                return (0, int(chrom_num))
            else:
                return (1, chrom_num)
        return (2, chrom)

    bed_df["sort_key"] = bed_df["chrom"].apply(chrom_sort_key)
    bed_df = bed_df.sort_values(["sort_key", "chromStart"])
    bed_df = bed_df.drop("sort_key", axis=1)

    # Write to output file
    bed_df.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Wrote {len(bed_entries)} BED12 entries to {output_file}")
    print(f"Added {terminal_exon_size}bp terminal exons to each transcript")


def main():
    """
    Main function with command line argument parsing.
    """
    parser = argparse.ArgumentParser(
        description="Convert introns data to BED12 format representing exons",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i input.tsv -o output.bed12
  %(prog)s -i data.tsv -o exons.bed12 --terminal-exon-size 50
        """,
    )

    parser.add_argument(
        "-i", "--input", required=True, help="Input TSV file with introns data"
    )

    parser.add_argument("-o", "--output", required=True, help="Output BED12 file")

    parser.add_argument(
        "--terminal-exon-size",
        type=int,
        default=25,
        help="Size of terminal exons to add in base pairs (default: 25)",
    )

    args = parser.parse_args()

    try:
        convert_to_bed12(args.input, args.output, args.terminal_exon_size)
        print("Conversion completed successfully!")
    except Exception as e:
        print(f"Error during conversion: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

# Example programmatic usage:
# convert_to_bed12("your_input_file.tsv", "your_output_file.bed12", terminal_exon_size=25)
