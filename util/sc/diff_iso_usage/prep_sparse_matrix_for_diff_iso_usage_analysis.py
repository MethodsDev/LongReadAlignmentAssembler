#!/usr/bin/env python3

import sys, os
import argparse
import logging

# Add pylib to path for Pipeliner import
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../../pylib"]))

from Pipeliner import Pipeliner, Command

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s : %(levelname)s : %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare sparse matrix for differential isoform usage analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--sparseM_dir", type=str, required=True, 
                       help="Path to sparse matrix directory")
    parser.add_argument("--clusters", type=str, required=True,
                       help="Path to cell cluster assignments file")
    parser.add_argument("--gene_transcript_splicehash", type=str, required=True,
                       help="Path to gene_transcript_splicehashcode file (with gene symbols)")
    
    args = parser.parse_args()
    
    sparseM_dir = os.path.abspath(args.sparseM_dir)
    clusters_file = os.path.abspath(args.clusters)
    gene_transcript_file = os.path.abspath(args.gene_transcript_splicehash)
    
    # Validate inputs
    if not os.path.exists(sparseM_dir):
        logger.error(f"Sparse matrix directory not found: {sparseM_dir}")
        sys.exit(1)
    if not os.path.exists(clusters_file):
        logger.error(f"Clusters file not found: {clusters_file}")
        sys.exit(1)
    if not os.path.exists(gene_transcript_file):
        logger.error(f"Gene transcript file not found: {gene_transcript_file}")
        sys.exit(1)
    
    # Create output directory
    output_dir = sparseM_dir + ".resources"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logger.info(f"Created output directory: {output_dir}")
    
    # Get the base name from sparseM_dir for output file naming
    sparseM_basename = os.path.basename(sparseM_dir)
    
    # Get script directory for relative path resolution
    scriptdir = os.path.abspath(os.path.dirname(__file__))
    util_sc_dir = os.path.dirname(scriptdir)  # util/sc directory
    
    # Define output file paths
    cluster_pseudobulk_matrix = os.path.join(output_dir, f"{sparseM_basename}.cluster_pseudobulk.matrix")
    cluster_pseudobulk_for_diff = os.path.join(output_dir, f"{sparseM_basename}.cluster_pseudobulk.matrix.for_diff_iso_usage")
    cluster_pseudobulk_iso_fraction = os.path.join(output_dir, f"{sparseM_basename}.cluster_pseudobulk.matrix.for_diff_iso_usage.iso_fraction.matrix")
    cell_fractions_matrix = os.path.join(output_dir, f"{sparseM_basename}.cell_fractions_expressed.matrix")
    cell_fractions_for_diff = os.path.join(output_dir, f"{sparseM_basename}.cell_fractions_expressed.matrix.for_diff_iso_usage")
    
    # Initialize Pipeliner with checkpoint directory in output dir
    checkpoint_dir = os.path.join(output_dir, "__chckpts")
    pipeliner = Pipeliner(checkpoint_dir)
    
    logger.info("Setting up differential isoform usage analysis pipeline...")
    
    # Step 1: Make cell clusters pseudobulk matrix
    cmd = (
        f"{os.path.join(util_sc_dir, 'sparse_matrix_to_cluster_pseudobulk_matrix.R')} "
        f"--sparseM_dir {sparseM_dir} "
        f"--clusters {clusters_file} "
        f"--output_matrix {cluster_pseudobulk_matrix}"
    )
    pipeliner.add_commands([Command(cmd, "cluster_pseudobulk_matrix.ok")])
    
    # Step 2: Prep cluster count matrix for diff iso usage
    cmd = (
        f"{os.path.join(util_sc_dir, 'sparse_matrix_to_cluster_pseudobulk_matrix.prep_for_diff_iso_usage.py')} "
        f"-e {cluster_pseudobulk_matrix} "
        f"-l {gene_transcript_file} "
        f"-o {cluster_pseudobulk_for_diff}"
    )
    pipeliner.add_commands([Command(cmd, "prep_cluster_for_diff_iso.ok")])
    
    # Step 3: Generate isoform fraction matrix from cluster pseudobulk counts
    cmd = (
        f"{os.path.join(util_sc_dir, 'cluster_pseudobulk_expr_matrix_to_gene_iso_fraction_matrix.py')} "
        f"-i {cluster_pseudobulk_for_diff} "
        f"-o {cluster_pseudobulk_iso_fraction}"
    )
    pipeliner.add_commands([Command(cmd, "cluster_iso_fraction.ok")])
    
    # Step 4: Get cell fraction expressed matrix
    cmd = (
        f"{os.path.join(util_sc_dir, 'sparse_matrix_to_cluster_pseudobulk_cell_fractions_expressed.R')} "
        f"--sparseM_dir {sparseM_dir} "
        f"--clusters {clusters_file} "
        f"--output_matrix {cell_fractions_matrix}"
    )
    pipeliner.add_commands([Command(cmd, "cell_fractions_matrix.ok")])
    
    # Step 5: Prep cell fractions matrix for diff iso usage
    cmd = (
        f"{os.path.join(util_sc_dir, 'sparse_matrix_to_cluster_pseudobulk_matrix.prep_for_diff_iso_usage.py')} "
        f"-e {cell_fractions_matrix} "
        f"-l {gene_transcript_file} "
        f"-o {cell_fractions_for_diff}"
    )
    pipeliner.add_commands([Command(cmd, "prep_fractions_for_diff_iso.ok")])
    
    # Run the pipeline
    logger.info("Running pipeline...")
    pipeliner.run()
    
    logger.info(f"Pipeline complete! Output files in: {output_dir}")
    logger.info(f"  - Cluster pseudobulk matrix: {cluster_pseudobulk_matrix}")
    logger.info(f"  - Cluster pseudobulk for diff iso: {cluster_pseudobulk_for_diff}")
    logger.info(f"  - Cluster pseudobulk iso fraction matrix: {cluster_pseudobulk_iso_fraction}")
    logger.info(f"  - Cell fractions matrix: {cell_fractions_matrix}")
    logger.info(f"  - Cell fractions for diff iso: {cell_fractions_for_diff}")
    
    sys.exit(0)


if __name__ == '__main__':
    main()
