#!/usr/bin/env python3

import argparse
import logging
import os
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import gzip
import shutil


def main():
    parser = argparse.ArgumentParser(
        description="Build Seurat-compatible 10X-format directory from unique read assignments."
    )
    parser.add_argument(
        "--input", required=True, help="Input TSV file with read assignments"
    )
    parser.add_argument(
        "--annotations", required=True, help="Annotation mapping TSV file"
    )
    parser.add_argument(
        "--outdir", required=True, help="Output directory for Seurat Read10X()"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )

    os.makedirs(args.outdir, exist_ok=True)

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

    # Merge in annotation mapping
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

    # Define barcodes and transcripts
    barcodes = sorted(df_merged["cell_barcode"].unique())
    transcripts = sorted(df_merged["new_transcript_splice_hash_code"].unique())
    barcode_index = {bc: i for i, bc in enumerate(barcodes)}
    transcript_index = {tr: j for j, tr in enumerate(transcripts)}

    # Build sparse matrix with features as rows, cells as columns
    logging.info("Building sparse matrix (features as rows, cells as columns)...")
    row_ind = (
        df_merged["new_transcript_splice_hash_code"].map(transcript_index).to_numpy()
    )
    col_ind = df_merged["cell_barcode"].map(barcode_index).to_numpy()
    data = np.ones(len(df_merged), dtype=np.int32)
    mat = csr_matrix(
        (data, (row_ind, col_ind)), shape=(len(transcripts), len(barcodes))
    )

    # Write matrix (gzipped)
    logging.info("Writing matrix.mtx.gz...")
    tmp_matrix = os.path.join(args.outdir, "matrix.mtx")
    mmwrite(tmp_matrix, mat)
    with open(tmp_matrix, "rb") as f_in, gzip.open(
        os.path.join(args.outdir, "matrix.mtx.gz"), "wb"
    ) as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(tmp_matrix)

    # Write barcodes
    logging.info("Writing barcodes.tsv.gz...")
    with gzip.open(os.path.join(args.outdir, "barcodes.tsv.gz"), "wt") as f:
        for bc in barcodes:
            f.write(bc + "\n")

    # Write features
    logging.info("Writing features.tsv.gz...")
    with gzip.open(os.path.join(args.outdir, "features.tsv.gz"), "wt") as f:
        for tr in transcripts:
            if "^" in tr:
                gene, _ = tr.split("^", 1)
            else:
                gene = tr
            f.write(f"{tr}\t{gene}\tGene\n")

    # Validation
    logging.info("Validating output consistency...")
    n_barcodes = sum(
        1 for _ in gzip.open(os.path.join(args.outdir, "barcodes.tsv.gz"), "rt")
    )
    n_features = sum(
        1 for _ in gzip.open(os.path.join(args.outdir, "features.tsv.gz"), "rt")
    )
    if mat.shape != (n_features, n_barcodes):
        raise ValueError(
            f"Matrix shape {mat.shape} does not match features ({n_features}) x barcodes ({n_barcodes})"
        )

    logging.info("All done! Point Seurat's Read10X() to this directory.")


if __name__ == "__main__":
    main()
