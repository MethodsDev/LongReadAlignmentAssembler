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
import gzip
import logging
import os
import re
import shutil
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


def load_symbol_to_ensg_mapping(mapping_file):
    """
    Load SYMBOL -> ENSG mapping from a TSV file.

    Format: ENSG_ID \t SYMBOL

    For symbols that map to multiple ENSG IDs, keeps the first occurrence
    and logs a warning.

    Args:
        mapping_file: Path to mapping file (can be .gz compressed)

    Returns:
        dict: {SYMBOL: ENSG_ID}
    """
    symbol_to_ensg = {}
    duplicate_symbols = {}

    open_func = gzip.open if mapping_file.endswith('.gz') else open
    mode = 'rt' if mapping_file.endswith('.gz') else 'r'

    with open_func(mapping_file, mode) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) != 2:
                continue

            ensg, symbol = parts

            if symbol in symbol_to_ensg:
                # Track duplicates
                if symbol not in duplicate_symbols:
                    duplicate_symbols[symbol] = [symbol_to_ensg[symbol]]
                duplicate_symbols[symbol].append(ensg)
            else:
                symbol_to_ensg[symbol] = ensg

    if duplicate_symbols:
        logging.warning(
            "Found %d symbols mapping to multiple ENSG IDs (using first occurrence)",
            len(duplicate_symbols)
        )
        # Show first 5 examples
        for symbol, ensg_list in list(duplicate_symbols.items())[:5]:
            logging.warning("  %s -> %s", symbol, ', '.join(ensg_list))

    logging.info("Loaded %d unique symbol->ENSG mappings", len(symbol_to_ensg))
    return symbol_to_ensg


def is_cas_compatible_format(matrix_dir):
    """
    Check if the matrix directory is already in CAS-compatible format.
    
    CAS expects:
      - barcodes.tsv (or .tsv.gz)
      - genes.tsv (or .tsv.gz)
      - matrix.mtx (or .mtx.gz)
    
    And genes.tsv should be formatted as:
      ENSG00000000003 TSPAN6
      ENSG00000000419 DPM1
      ...
    
    Returns:
        bool: True if compatible, False if it needs conversion
    """
    # Check for genes.tsv or genes.tsv.gz
    genes_file = None
    if os.path.exists(os.path.join(matrix_dir, "genes.tsv")):
        genes_file = os.path.join(matrix_dir, "genes.tsv")
    elif os.path.exists(os.path.join(matrix_dir, "genes.tsv.gz")):
        genes_file = os.path.join(matrix_dir, "genes.tsv.gz")
    else:
        # No genes.tsv, probably features.tsv format
        return False
    
    # Check first line format
    try:
        if genes_file.endswith(".gz"):
            with gzip.open(genes_file, "rt") as f:
                first_line = f.readline().strip()
        else:
            with open(genes_file, "r") as f:
                first_line = f.readline().strip()
        
        # Expected format: "ENSG00000000003 TSPAN6" or "ENSG00000000003\tTSPAN6"
        parts = first_line.split()
        if len(parts) == 2 and parts[0].startswith("ENSG"):
            return True
    except Exception as e:
        logging.warning("Could not check genes file format: %s", e)
    
    return False


def _find_file(directory, base_name):
    """
    Find a file in directory, checking both compressed and uncompressed versions.
    
    Args:
        directory: Directory to search in
        base_name: Base filename (e.g., "barcodes.tsv")
    
    Returns:
        tuple: (file_path, is_compressed) or raises FileNotFoundError
    """
    compressed_path = os.path.join(directory, f"{base_name}.gz")
    uncompressed_path = os.path.join(directory, base_name)
    
    if os.path.exists(compressed_path):
        return compressed_path, True
    elif os.path.exists(uncompressed_path):
        return uncompressed_path, False
    else:
        raise FileNotFoundError(
            f"Could not find {base_name} or {base_name}.gz in {directory}"
        )


def convert_lraa_matrix_to_cas_format(lraa_matrix_dir, output_dir, symbol_to_ensg_mapping=None):
    """
    Convert LRAA-formatted matrix to CAS-compatible format.

    LRAA format:
      - barcodes.tsv.gz
      - features.tsv.gz (format: SYMBOL^ENSG00000000000.version or SYMBOL^LRAA_gene_id)
      - matrix.mtx.gz

    CAS format:
      - barcodes.tsv
      - genes.tsv (format: ENSG00000000000 SYMBOL)
      - matrix.mtx

    Args:
        lraa_matrix_dir: Path to LRAA-formatted matrix directory
        output_dir: Path to output CAS-compatible directory
        symbol_to_ensg_mapping: Optional dict mapping gene symbols to ENSG IDs

    Returns:
        str: Path to the converted matrix directory
    """
    logging.info("Converting LRAA matrix format to CAS-compatible format...")
    logging.info("Output directory will be: %s", output_dir)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Convert barcodes.tsv[.gz] -> barcodes.tsv
    barcodes_in, is_compressed = _find_file(lraa_matrix_dir, "barcodes.tsv")
    barcodes_out = os.path.join(output_dir, "barcodes.tsv")
    
    logging.info("Converting %s...", os.path.basename(barcodes_in))
    open_func = gzip.open if is_compressed else open
    mode_in = "rt" if is_compressed else "r"
    with open_func(barcodes_in, mode_in) as f_in:
        with open(barcodes_out, "w") as f_out:
            f_out.write(f_in.read())
    
    # 2. Convert features.tsv[.gz] -> genes.tsv
    features_in, is_compressed = _find_file(lraa_matrix_dir, "features.tsv")
    genes_out = os.path.join(output_dir, "genes.tsv")

    logging.info("Converting %s to genes.tsv...", os.path.basename(features_in))

    # Track mapping statistics
    total_features = 0
    mapped_via_ensg = 0  # Already had ENSG
    mapped_via_lookup = 0  # Mapped via symbol lookup
    unmapped = 0  # Could not map

    open_func = gzip.open if is_compressed else open
    mode_in = "rt" if is_compressed else "r"
    with open_func(features_in, mode_in) as f_in:
        with open(genes_out, "w") as f_out:
            for line in f_in:
                line = line.strip()
                if not line:
                    continue

                total_features += 1

                # Parse SYMBOL^GENE_ID format
                if "^" in line:
                    parts = line.split("^")
                    if len(parts) == 2:
                        symbol = parts[0]
                        gene_id = parts[1]

                        # Check if gene_id is already ENSG format
                        if gene_id.startswith("ENSG"):
                            # Remove version number (e.g., ENSG00000000000.15 -> ENSG00000000000)
                            ensg = re.sub(r"\.\d+$", "", gene_id)
                            mapped_via_ensg += 1
                        elif symbol_to_ensg_mapping and symbol in symbol_to_ensg_mapping:
                            # Look up symbol in mapping
                            ensg = symbol_to_ensg_mapping[symbol]
                            mapped_via_lookup += 1
                        else:
                            # No mapping found, use gene_id as-is (for novel/LRAA genes)
                            ensg = gene_id
                            unmapped += 1

                        # Write in CAS format: ENSG SYMBOL
                        f_out.write(f"{ensg}\t{symbol}\n")
                    else:
                        logging.warning("Unexpected features.tsv format: %s", line)
                else:
                    # No '^' separator, treat the full token as both symbol and gene ID
                    token = line

                    # Check if it's an ENSG ID
                    if token.startswith("ENSG"):
                        ensg = re.sub(r"\.\d+$", "", token)
                        mapped_via_ensg += 1
                    elif symbol_to_ensg_mapping and token in symbol_to_ensg_mapping:
                        ensg = symbol_to_ensg_mapping[token]
                        mapped_via_lookup += 1
                    else:
                        ensg = token
                        unmapped += 1

                    f_out.write(f"{ensg}\t{token}\n")

    # Log mapping statistics
    logging.info("Feature mapping statistics:")
    logging.info("  Total features: %d", total_features)
    logging.info("  Already had ENSG ID: %d (%.1f%%)",
                 mapped_via_ensg,
                 100.0 * mapped_via_ensg / total_features if total_features > 0 else 0)
    logging.info("  Mapped via symbol lookup: %d (%.1f%%)",
                 mapped_via_lookup,
                 100.0 * mapped_via_lookup / total_features if total_features > 0 else 0)
    logging.info("  Unmapped (kept original ID): %d (%.1f%%)",
                 unmapped,
                 100.0 * unmapped / total_features if total_features > 0 else 0)
    
    # 3. Convert matrix.mtx[.gz] -> matrix.mtx
    matrix_in, is_compressed = _find_file(lraa_matrix_dir, "matrix.mtx")
    matrix_out = os.path.join(output_dir, "matrix.mtx")
    
    logging.info("Converting %s...", os.path.basename(matrix_in))
    open_func = gzip.open if is_compressed else open
    mode_in = "rt" if is_compressed else "r"
    with open_func(matrix_in, mode_in) as f_in:
        with open(matrix_out, "w") as f_out:
            f_out.write(f_in.read())
    
    logging.info("Conversion complete. CAS-compatible matrix saved to: %s", output_dir)
    return output_dir


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

    # Load symbol to ENSG mapping
    mapping_file = "/usr/local/share/cas/human.refdata-gex-GRCh38-GENCODE47.ENSG_to_SYMBOL.tsv.gz"
    symbol_to_ensg = None
    if os.path.exists(mapping_file):
        logging.info("Loading gene symbol to ENSG mapping from %s", mapping_file)
        symbol_to_ensg = load_symbol_to_ensg_mapping(mapping_file)
    else:
        logging.warning("Symbol to ENSG mapping file not found at %s", mapping_file)
        logging.warning("Will only use ENSG IDs already present in features file")

    # Check if matrix is already CAS-compatible or needs conversion
    cas_matrix_dir = matrix_dir

    if is_cas_compatible_format(matrix_dir):
        logging.info("Matrix directory is already in CAS-compatible format.")
    else:
        logging.info("Matrix directory appears to be in LRAA format. Converting to CAS format...")
        # Create output directory in same location as input with .for_CAS extension
        matrix_basename = os.path.basename(os.path.normpath(matrix_dir))
        matrix_parent = os.path.dirname(os.path.normpath(matrix_dir))
        cas_output_dir = os.path.join(matrix_parent, f"{matrix_basename}.for_CAS")
        cas_matrix_dir = convert_lraa_matrix_to_cas_format(matrix_dir, cas_output_dir, symbol_to_ensg)
        logging.info("Using converted matrix from: %s", cas_matrix_dir)

    # 1) Read 10x matrix
    logging.info("Reading 10x matrix...")
    adata = sc.read_10x_mtx(cas_matrix_dir, var_names="gene_symbols", cache=True)
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
