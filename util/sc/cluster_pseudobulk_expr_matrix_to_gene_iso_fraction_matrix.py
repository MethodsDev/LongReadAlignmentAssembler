#!/usr/bin/env python3

import argparse
import sys
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Convert transcript expression table to per-gene transcript fractions."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input tab-delimited expression table"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output tab-delimited fraction table"
    )
    args = parser.parse_args()

    # Load the table
    df = pd.read_csv(args.input, sep="\t")

    # Check first two columns
    expected_cols = ["gene_id", "transcript_id"]
    for i, expected in enumerate(expected_cols):
        if df.columns[i] != expected:
            sys.exit(
                f"ERROR: Expected column {i+1} to be '{expected}', "
                f"but found '{df.columns[i]}'. Please check your input file."
            )

    # Identify expression/sample columns
    id_cols = df.columns[:2].tolist()
    expr_cols = df.columns[2:].tolist()

    # Group by gene_id, compute gene totals per sample column
    gene_totals = df.groupby("gene_id")[expr_cols].transform("sum")

    # Divide transcript expression by gene total
    df_fraction = df.copy()
    df_fraction[expr_cols] = df[expr_cols] / gene_totals[expr_cols]

    # Save to output file
    df_fraction.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
