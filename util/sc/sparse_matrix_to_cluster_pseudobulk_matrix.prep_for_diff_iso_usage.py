#!/usr/bin/env python3

import pandas as pd
import argparse
import sys


# thanks claude.ai ! min mods by bhaas


def create_lookup_dictionary(lookup_file):
    """
    Create a lookup dictionary from the lookup table file.

    Args:
        lookup_file: Path to the lookup table file

    Returns:
        dict: Dictionary mapping transcript_id to gene_id
    """
    print("Reading lookup table...")
    lookup_df = pd.read_csv(lookup_file, sep="\t")

    print(f"Lookup table shape: {lookup_df.shape}")
    print(f"Lookup table columns: {lookup_df.columns.tolist()}")

    # Create lookup dictionaries for efficient searching
    print("Creating lookup dictionaries...")

    # Dictionary for transcript_id -> gene_id
    transcript_to_gene = {}

    # Add mappings from transcript_id column
    for _, row in lookup_df.iterrows():
        transcript_to_gene[row["transcript_id"]] = row["gene_id"]

        # Add mappings from transcript_splice_hash_code column
        if pd.notna(row["transcript_splice_hash_code"]):
            transcript_to_gene[row["transcript_splice_hash_code"]] = row["gene_id"]

        # Add mappings from new_transcript_id column to new_gene_id
        if pd.notna(row["new_transcript_id"]):
            transcript_to_gene[row["new_transcript_id"]] = row["new_gene_id"]

        # Add mappings from new_transcript_splice_hash_code column to new_gene_id
        if pd.notna(row["new_transcript_splice_hash_code"]):
            transcript_to_gene[row["new_transcript_splice_hash_code"]] = row[
                "new_gene_id"
            ]

    print(f"Created lookup dictionary with {len(transcript_to_gene)} entries")

    return transcript_to_gene


def process_expression_matrix(expression_file, lookup_file, output_file):
    """
    Process expression matrix to add gene_id column and fix header alignment.

    Args:
        expression_file: Path to the expression matrix file
        lookup_file: Path to the lookup table file
        output_file: Path for the output file
    """

    # Read the expression matrix
    print("Reading expression matrix...")

    # Read the data with the first column as index (transcript_id)
    expr_df = pd.read_csv(expression_file, sep="\t", index_col=0)

    # Reset index to make transcript_id a proper column
    expr_df.reset_index(inplace=True)
    expr_df.rename(columns={"index": "transcript_id"}, inplace=True)

    print(f"Expression matrix shape: {expr_df.shape}")
    print(f"First few transcript IDs: {expr_df['transcript_id'].head().tolist()}")

    # Create lookup dictionary
    transcript_to_gene = create_lookup_dictionary(lookup_file)

    # Function to lookup gene_id for a given transcript_id
    def lookup_gene_id(transcript_id):
        return transcript_to_gene.get(transcript_id, "NOT_FOUND")

    # Apply lookup to get gene_id for each transcript
    print("Looking up gene IDs...")
    expr_df["gene_id"] = expr_df["transcript_id"].apply(lookup_gene_id)

    # Reorder columns to put gene_id first, then transcript_id, then expression data
    columns = ["gene_id", "transcript_id"] + [
        col for col in expr_df.columns if col not in ["gene_id", "transcript_id"]
    ]
    expr_df = expr_df[columns]

    # Check for any missing gene_ids
    missing_count = (expr_df["gene_id"] == "NOT_FOUND").sum()
    if missing_count > 0:
        print(
            f"Warning: {missing_count} transcript IDs could not be matched to gene IDs"
        )
        missing_transcripts = expr_df[expr_df["gene_id"] == "NOT_FOUND"][
            "transcript_id"
        ].tolist()
        print(f"First few missing transcript IDs: {missing_transcripts[:5]}")

    # Save the processed matrix
    print(f"Saving processed matrix to {output_file}...")
    expr_df.to_csv(output_file, sep="\t", index=False)

    print(f"Processing complete!")
    print(f"Final matrix shape: {expr_df.shape}")
    print(f"Columns: {expr_df.columns.tolist()}")

    return expr_df


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Process expression matrix to add gene_id column and fix header alignment.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py -e expression_matrix.txt -l lookup_table.txt -o output.txt
  python script.py --expression expression_matrix.txt --lookup lookup_table.txt --output output.txt
        """,
    )

    parser.add_argument(
        "-e", "--expression", required=True, help="Path to the expression matrix file"
    )

    parser.add_argument(
        "-l", "--lookup", required=True, help="Path to the lookup table file"
    )

    parser.add_argument(
        "-o", "--output", required=True, help="Path for the output file"
    )

    return parser.parse_args()


# Example usage
if __name__ == "__main__":
    args = parse_arguments()

    try:
        processed_df = process_expression_matrix(
            args.expression, args.lookup, args.output
        )

        # Show a preview of the results
        print("\nPreview of processed data:")
        print(processed_df.head())

    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        print("Please make sure the input files exist and the paths are correct")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing files: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
