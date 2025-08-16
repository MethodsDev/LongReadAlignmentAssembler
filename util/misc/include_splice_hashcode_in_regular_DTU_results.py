#!/usr/bin/env python3


import pandas as pd
import argparse


def merge_dtu_with_hash_codes(dtu_file, annot_file, output_file):
    # Read the TSV files
    dtu_df = pd.read_csv(dtu_file, sep="\t")
    annot_df = pd.read_csv(annot_file, sep="\t")

    print(f"DTU results shape: {dtu_df.shape}")
    print(f"Annotation mappings shape: {annot_df.shape}")
    print(f"DTU columns: {list(dtu_df.columns)}")
    print(f"Annotation columns: {list(annot_df.columns)}")

    # Create mapping dictionary from new_transcript_id to new_transcript_splice_hash_code
    transcript_to_hash = dict(
        zip(annot_df["new_transcript_id"], annot_df["new_transcript_splice_hash_code"])
    )

    print(f"Created mapping for {len(transcript_to_hash)} transcripts")

    # Function to map transcript IDs to hash codes
    def map_transcript_ids_to_hashes(transcript_ids):
        if pd.isna(transcript_ids):
            return "N/A"

        # Split by comma and strip whitespace
        ids = [tid.strip() for tid in str(transcript_ids).split(",")]

        # Map each ID to its hash code
        hash_codes = [transcript_to_hash.get(tid, "N/A") for tid in ids]

        return ",".join(hash_codes)

    # Add new columns with hash codes
    dtu_df["dominant_transcript_hash_codes"] = dtu_df["dominant_transcript_ids"].apply(
        map_transcript_ids_to_hashes
    )
    dtu_df["alternate_transcript_hash_codes"] = dtu_df[
        "alternate_transcript_ids"
    ].apply(map_transcript_ids_to_hashes)

    print("\nSample of enhanced data:")
    print(
        dtu_df[
            [
                "gene_id",
                "dominant_transcript_ids",
                "dominant_transcript_hash_codes",
                "alternate_transcript_ids",
                "alternate_transcript_hash_codes",
            ]
        ].head()
    )

    # Save the enhanced results
    dtu_df.to_csv(output_file, sep="\t", index=False)

    print(f"\nEnhanced DTU results saved to '{output_file}'")
    print(f"Original columns: {len(dtu_df.columns) - 2}")
    print(f"Enhanced columns: {len(dtu_df.columns)}")

    return dtu_df


# Command line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge DTU results with transcript hash codes"
    )

    parser.add_argument(
        "--dtu_file", required=True, help="Path to DTU results TSV file"
    )
    parser.add_argument(
        "--annot_mappings_file",
        required=True,
        help="Path to annotation mappings TSV file",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="DTU_results_with_hash_codes.tsv",
        help="Output file path (default: DTU_results_with_hash_codes.tsv)",
    )

    args = parser.parse_args()

    # Run the merge
    enhanced_dtu = merge_dtu_with_hash_codes(
        args.dtu_file, args.annot_mappings_file, args.output
    )
