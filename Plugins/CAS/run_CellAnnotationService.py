#!/usr/bin/env python3
"""
Run Cellarium Cell Annotation Service (CAS) on a 10x-style sparse matrix folder.

This script:
  1. Reads a 10x Genomics matrix directory (matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz])
  2. Calls Cellarium CAS for ontology-aware cell type annotation
  3. Inserts the CAS response into the AnnData object
  4. Computes most granular top-k cell type calls per cell
  5. Writes:
       - <output_prefix>.h5ad
       - <output_prefix>.tsv (cell-barcode + all obs columns)
"""

import argparse
import logging
import os
import sys

import scanpy as sc
from cellarium.cas.client import CASClient
from cellarium.cas.postprocessing import insert_cas_ontology_aware_response_into_adata
import cellarium.cas.postprocessing.ontology_aware as pp
from cellarium.cas.postprocessing.cell_ontology import CellOntologyCache


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Cellarium CAS cell type annotation on a 10x sparse matrix folder."
    )
    parser.add_argument(
        "--matrix-dir",
        required=True,
        help="Path to 10x-style matrix directory (matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz]).",
    )
    parser.add_argument(
        "--api-token",
        required=True,
        help="Cellarium CAS API token (REQUIRED).",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Prefix used to name output files: <prefix>.h5ad and <prefix>.tsv",
    )
    parser.add_argument(
        "--cas-model-name",
        default=None,
        help="Optional CAS model name. Leave empty to use default.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500,
        help="Number of cells per CAS chunk (default: 500).",
    )
    parser.add_argument(
        "--min-acceptable-score",
        type=float,
        default=0.2,
        help="Minimum acceptable evidence score (default: 0.2).",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=3,
        help="Top-k CAS cell type calls to store per cell (default: 3).",
    )
    parser.add_argument(
        "--obs-prefix",
        default="cas_cell_type",
        help="Prefix for CAS cell type columns in .obs (default: cas_cell_type).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
    )

    matrix_dir = args.matrix_dir
    token = args.api_token
    output_prefix = args.output_prefix

    out_h5ad = f"{output_prefix}.h5ad"
    out_tsv = f"{output_prefix}.tsv"

    logging.info("Matrix directory: %s", matrix_dir)
    if not os.path.isdir(matrix_dir):
        logging.error("Matrix directory does not exist: %s", matrix_dir)
        sys.exit(1)

    # 1) Read 10x matrix
    logging.info("Reading 10x matrix...")
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=True)
    logging.info("Loaded AnnData: %d cells × %d genes", adata.n_obs, adata.n_vars)

    # 2) CAS client
    logging.info("Initializing CAS client...")
    cas = CASClient(api_token=token)

    # 3) CAS annotation request
    logging.info("Submitting to CAS ...")
    cas_response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        matrix=adata,
        chunk_size=args.chunk_size,
        feature_ids_column_name="gene_ids",
        feature_names_column_name="index",
        cas_model_name=args.cas_model_name,
    )
    logging.info("CAS response received with %d entries.", len(cas_response.data))

    # 4) Insert CAS response into AnnData
    logging.info("Integrating CAS response into AnnData...")
    insert_cas_ontology_aware_response_into_adata(cas_response, adata)

    # 5) Ontology
    logging.info("Loading Cell Ontology Cache...")
    cl = CellOntologyCache()

    logging.info("Computing most granular top-%d calls...", args.top_k)
    pp.compute_most_granular_top_k_calls_single(
        adata=adata,
        cl=cl,
        min_acceptable_score=args.min_acceptable_score,
        top_k=args.top_k,
        obs_prefix=args.obs_prefix,
    )

    # 6) Save outputs
    logging.info("Writing AnnData to %s", out_h5ad)
    adata.write(out_h5ad)

    logging.info("Writing cell-type table to %s", out_tsv)
    adata.obs.to_csv(out_tsv, sep="\t", index=True, index_label="cell_barcode")

    logging.info("Done!")
    logging.info("  H5AD: %s", out_h5ad)
    logging.info("  TSV : %s", out_tsv)


if __name__ == "__main__":
    main()
