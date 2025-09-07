#!/usr/bin/env python3

import argparse
import logging
import pandas as pd
from scipy.sparse import csr_matrix
import numpy as np
from scipy.io import mmwrite


def main():
    parser = argparse.ArgumentParser(
        description="Build sparse cell x transcript_splice_hash_code matrix from unique reads."
    )
    parser.add_argument(
        "--input", required=True, help="Input TSV file with read assignments"
    )
    parser.add_argument(
        "--annotations", required=True, help="Annotation mapping TSV file"
    )
    parser.add_argument(
        "--out_prefix",
        required=True,
        help="Output prefix (will write .mtx, .barcodes.txt, .features.txt)",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )

    # Read files
    logging.info("Reading input read-assignments file...")
    df = pd.read_csv(args.input, sep="\t")

    logging.info("Reading annotation mapping file...")
    anno = pd.read_csv(args.annotations, sep="\t")

    # Identify unique reads
    logging.info(
        "Identifying unique read_name → transcript_splice_hash_code assignments..."
    )
    read_counts = df.groupby("read_name")["transcript_splice_hash_code"].nunique()
    unique_reads = read_counts[read_counts == 1].index
    df_unique = df[df["read_name"].isin(unique_reads)]

    logging.info(
        f"Retained {len(df_unique)} unique read assignments out of {len(df)} total records."
    )

    # Add cell barcode
    logging.info("Extracting cell barcodes from read_name...")
    df_unique["cell_barcode"] = df_unique["read_name"].str.split("^").str[0]

    # Merge in annotation mapping to swap in new_transcript_splice_hash_code
    logging.info("Merging annotation mapping...")
    df_merged = df_unique.merge(
        anno[["transcript_splice_hash_code", "new_transcript_splice_hash_code"]],
        on="transcript_splice_hash_code",
        how="left",
    )

    if df_merged["new_transcript_splice_hash_code"].isnull().any():
        logging.warning(
            "Some transcript_splice_hash_code values did not have a mapping!"
        )

    # Build sparse matrix
    logging.info("Building sparse matrix...")
    barcodes = df_merged["cell_barcode"].unique()
    transcripts = df_merged["new_transcript_splice_hash_code"].unique()

    barcode_index = {bc: i for i, bc in enumerate(barcodes)}
    transcript_index = {tr: j for j, tr in enumerate(transcripts)}

    row_ind = df_merged["cell_barcode"].map(barcode_index).to_numpy()
    col_ind = (
        df_merged["new_transcript_splice_hash_code"].map(transcript_index).to_numpy()
    )
    data = np.ones(len(df_merged), dtype=np.int32)

    mat = csr_matrix(
        (data, (row_ind, col_ind)), shape=(len(barcodes), len(transcripts))
    )

    # Write outputs
    logging.info("Writing outputs in 10X format...")
    mmwrite(args.out_prefix + ".mtx", mat)

    with open(args.out_prefix + ".barcodes.txt", "w") as f:
        for bc in barcodes:
            f.write(bc + "\n")

    # features.txt with 3-column format: feature_id, gene_name, feature_type
    with open(args.out_prefix + ".features.txt", "w") as f:
        for tr in transcripts:
            # split "GENE^hashcode" → gene, id
            if "^" in tr:
                gene, _ = tr.split("^", 1)
            else:
                gene = tr
            f.write(f"{tr}\t{gene}\tGene\n")

    logging.info("Done.")


if __name__ == "__main__":
    main()
