version 1.0

import "LRAA.wdl" as LRAA
import "subwdls/LRAA-build_sparse_matrices_from_tracking.wdl" as BuildMatrices
import "subwdls/LRAA-filter_good_cells.wdl" as FilterCells
import "subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl" as Seurat
import "subwdls/Incorporate_gene_symbols.wdl" as GeneSymbols
import "LRAA-cell_cluster_guided.wdl" as ClusterGuided

workflow LRAA_singlecell_wf {
  input {
    # Core inputs
    String sample_id
    File referenceGenome
    File inputBAM

    # Optional: guide initial discovery with an annotation
    File? initial_annot_gtf

    # Single-cell pipeline mode: 'basic' or 'cluster-guided'
    String single_cell_pipe_mode = "basic"
    
    # Quantification-only mode: skip discovery, requires initial_annot_gtf
    Boolean quant_only = false

    # Platform/options
    Boolean HiFi = false
    String oversimplify = "chrM"   # e.g., "chrM" or "chrM,M"
    String main_chromosomes = ""  # if empty, runs without partitioning
    String? region                 # e.g., "chr1:100000-200000"; forces direct mode

    # Optional: reuse outputs from a prior initial discovery run and skip LRAA_init
    File? precomputed_init_quant_tracking
    File? precomputed_init_gtf
    
    # Optional: provide precomputed cluster assignments to skip steps 1-3
    File? precomputed_cluster_assignments_tsv

    # Single-cell barcode and UMI tags
    String cell_barcode_tag = "CB"
    String read_umi_tag = "XM"

    # Resources and docker (propagated to subcalls where applicable)
    Int numThreadsPerWorker = 5
    Int numThreadsPerWorkerScattered = 5
    Int num_parallel_contigs = 3
    Int memoryGB = 64
    Int memoryGBPerWorkerScattered = 32
    Int memoryGBbuildSparseMatrices = 32
    Int memoryGBFilterCells = 32
    Int memoryGBSeurat = 32
    Int memoryGBmergeGTFs = 32
    Int memoryGBquantFinal = 32
    Int memoryGBscSparseMatrices = 32
    Int diskSizeGB = 256
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"

    # Seurat clustering parameters (forwarded to Seurat subworkflow)
    Int min_cells = 10
    Int min_features = 1000
    Float percent_mt_max = 20.0
    String mt_pattern = "^MT-"
    Int npcs = 12
    Float resolution = 0.6
    Int n_variable_features = 2000
    Int seed = 1

    # Filter good cells parameters
    Float fdr_threshold = 0.01
    Int? lower_threshold

    File? ref_annot_gtf_source_gene_symbols # used as source for gene symbol assignment at the end.
  }

  # Only skip the initial discovery call when the downstream-critical artifacts are provided
  Boolean has_precomputed_init = defined(precomputed_init_gtf) && (defined(precomputed_init_quant_tracking) || defined(precomputed_cluster_assignments_tsv))
  Boolean has_precomputed_clusters = defined(precomputed_cluster_assignments_tsv)
  Boolean run_initial_phase = !has_precomputed_init
  Boolean run_clustering_phase = !has_precomputed_clusters
  Boolean run_cluster_guided = single_cell_pipe_mode == "cluster-guided"

  # 1) Initial transcript discovery + quantification on the full BAM (skipped when precomputed inputs are provided)
  if (run_initial_phase) {
    call LRAA.LRAA_wf as LRAA_init {
      input:
        sample_id = sample_id,
        referenceGenome = referenceGenome,
        inputBAM = inputBAM,
        annot_gtf = initial_annot_gtf,
        HiFi = HiFi,
        oversimplify = oversimplify,
        main_chromosomes = main_chromosomes,
        region = region,
        cell_barcode_tag = cell_barcode_tag,
        read_umi_tag = read_umi_tag,
        numThreadsPerWorker = numThreadsPerWorker,
        numThreadsPerWorkerScattered = numThreadsPerWorkerScattered,
        num_parallel_contigs = num_parallel_contigs,
        memoryGB = memoryGB,
        memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
        diskSizeGB = diskSizeGB,
        docker = docker,
        quant_only = quant_only
    }
  }

  File? init_quant_expr_file = LRAA_init.mergedQuantExpr
  File? init_quant_tracking_generated = LRAA_init.mergedQuantTracking
  File? init_quant_tracking_file = if (!run_initial_phase && defined(precomputed_init_quant_tracking)) then select_first([precomputed_init_quant_tracking]) else init_quant_tracking_generated
  File? init_gtf_generated = LRAA_init.mergedGTF
  File? init_gtf_file = if (!run_initial_phase && defined(precomputed_init_gtf)) then select_first([precomputed_init_gtf]) else init_gtf_generated

  # 2) Build single-cell sparse matrices from the initial tracking (skipped when precomputed clusters are provided)
  if (run_clustering_phase) {
    call BuildMatrices.BuildSparseMatricesFromTracking as build_sc_from_init_tracking {
      input:
        sample_id = sample_id,
        tracking_file = select_first([init_quant_tracking_file]),
        docker = docker,
        memoryGB = memoryGBbuildSparseMatrices
    }

    # 2.5) Filter good cells from the gene-level sparse matrix
    call FilterCells.FilterGoodCells as filter_good_cells {
      input:
        sample_id = sample_id,
        gene_sparse_tar_gz = build_sc_from_init_tracking.gene_sparse_dir_tgz,
        isoform_sparse_tar_gz = build_sc_from_init_tracking.isoform_sparse_dir_tgz,
        splice_pattern_sparse_tar_gz = build_sc_from_init_tracking.splice_pattern_sparse_dir_tgz,
        docker = docker,
        memoryGB = memoryGBFilterCells,
        fdr_threshold = fdr_threshold,
        lower_threshold = lower_threshold
    }

    # 3) Cluster cells from the filtered gene-level sparse matrix
    call Seurat.GeneSparseM_To_SeuratClusters as cluster_cells {
      input:
        sample_id = sample_id,
        gene_sparse_tar_gz = filter_good_cells.filtered_gene_sparse_tar_gz,
        docker = docker,
        memoryGB = memoryGBSeurat,
        min_cells = min_cells,
        min_features = min_features,
        percent_mt_max = percent_mt_max,
        mt_pattern = mt_pattern,
        npcs = npcs,
        resolution = resolution,
        n_variable_features = n_variable_features,
        seed = seed
    }
  }

  # Select cluster assignments: use precomputed if provided, otherwise use generated
  File? cluster_assignments_generated = cluster_cells.cluster_assignments_tsv
  File cluster_assignments_file = if (has_precomputed_clusters) then select_first([precomputed_cluster_assignments_tsv]) else select_first([cluster_assignments_generated])

  # 4) Cluster-guided reconstruction + final single-cell matrices (only if mode is 'cluster-guided')
  if (run_cluster_guided) {
    call ClusterGuided.LRAA_cell_cluster_guided as cluster_guided {
      input:
        sample_id = sample_id,
        referenceGenome = referenceGenome,
        inputBAM = inputBAM,
        cell_clusters_info = cluster_assignments_file,
        annot_gtf = if quant_only then initial_annot_gtf else init_gtf_file,
        HiFi = HiFi,
        oversimplify = oversimplify,
        main_chromosomes = main_chromosomes,
        cell_barcode_tag = cell_barcode_tag,
        read_umi_tag = read_umi_tag,
        numThreadsPerWorker = numThreadsPerWorker,
        numThreadsPerWorkerScattered = numThreadsPerWorkerScattered,
        num_parallel_contigs = num_parallel_contigs,
        memoryGB = memoryGB,
        memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
        memoryGBmergeGTFs = memoryGBmergeGTFs,
        memoryGBquantFinal = memoryGBquantFinal,
        memoryGBscSparseMatrices = memoryGBscSparseMatrices,
        diskSizeGB = diskSizeGB,
        docker = docker,
        quant_only_cluster_guided = quant_only
    }
  }

  # 5) Incorporate gene symbols: use cluster-guided outputs if available, otherwise use filtered good cell outputs
  # In quant_only mode, these GTF values will naturally be undefined since discovery doesn't produce them
  File? gtf_for_symbols = if run_cluster_guided then cluster_guided.LRAA_final_gtf else init_gtf_file
  File? gene_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_gene_sparse_tar_gz else filter_good_cells.filtered_gene_sparse_tar_gz
  File? isoform_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_isoform_sparse_tar_gz else filter_good_cells.filtered_isoform_sparse_tar_gz
  File? splice_pattern_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_splice_pattern_sparse_tar_gz else filter_good_cells.filtered_splice_pattern_sparse_tar_gz
  File? mapping_for_symbols = if run_cluster_guided then cluster_guided.sc_gene_transcript_splicehash_mapping else build_sc_from_init_tracking.mapping_file

  if (defined(ref_annot_gtf_source_gene_symbols)) {
    call GeneSymbols.Incorporate_gene_symbols as add_gene_symbols {
      input:
        sample_id = sample_id,
        reference_gtf = select_first([ref_annot_gtf_source_gene_symbols]),
        final_gtf = gtf_for_symbols,
        final_sc_gene_sparse_tar_gz = select_first([gene_sparse_for_symbols]),
        final_sc_isoform_sparse_tar_gz = select_first([isoform_sparse_for_symbols]),
        final_sc_splice_pattern_sparse_tar_gz = select_first([splice_pattern_sparse_for_symbols]),
        final_sc_gene_transcript_splicehash_mapping = select_first([mapping_for_symbols]),
        docker = docker
    }
  }

  output {
    # Initial discovery outputs
    File? init_quant_expr = init_quant_expr_file
    File? init_quant_tracking = init_quant_tracking_generated
    File? init_gtf = init_gtf_generated

    # Initial single-cell matrices and clustering inputs/outputs
    File? init_sc_gene_sparse_tar_gz = build_sc_from_init_tracking.gene_sparse_dir_tgz
    File? init_sc_isoform_sparse_tar_gz = build_sc_from_init_tracking.isoform_sparse_dir_tgz
    File? init_sc_splice_pattern_sparse_tar_gz = build_sc_from_init_tracking.splice_pattern_sparse_dir_tgz
    File? init_sc_gene_transcript_splicehash_mapping = build_sc_from_init_tracking.mapping_file
    
    
    # Filter empty droplets
    File? filtered_gene_sparse_tar_gz = filter_good_cells.filtered_gene_sparse_tar_gz
    File? filtered_isoform_sparse_tar_gz = filter_good_cells.filtered_isoform_sparse_tar_gz
    File? filtered_splice_pattern_sparse_tar_gz = filter_good_cells.filtered_splice_pattern_sparse_tar_gz
    File? good_cell_barcodes = filter_good_cells.good_cell_barcodes
    File? filtering_summary = filter_good_cells.filtering_summary
    
    # Seurat for gene-based filtering, cell clustering, and umap
    File? seurat_umap_pdf = cluster_cells.umap_pdf
    File? seurat_umap_with_clusters_tsv = cluster_cells.umap_with_clusters_tsv
    File? seurat_rds_initial = cluster_cells.seurat_rds_initial
    File? seurat_rds = cluster_cells.seurat_rds
    File? seurat_cluster_assignments = cluster_assignments_generated

    # Final cluster-guided LRAA outputs (main deliverables)
    File? final_gtf = cluster_guided.LRAA_final_gtf
    File? final_gtf_tracking = cluster_guided.LRAA_final_gtf_tracking
    File? final_tracking = cluster_guided.LRAA_final_tracking
    File? final_sc_gene_sparse_tar_gz = cluster_guided.sc_gene_sparse_tar_gz
    File? final_sc_isoform_sparse_tar_gz = cluster_guided.sc_isoform_sparse_tar_gz
    File? final_sc_splice_pattern_sparse_tar_gz = cluster_guided.sc_splice_pattern_sparse_tar_gz
    File? final_sc_gene_transcript_splicehash_mapping = cluster_guided.sc_gene_transcript_splicehash_mapping

    # Convenience tarballs and matrices from the cluster-guided phase
    File? partitioned_cluster_bams_tar = cluster_guided.LRAA_partitioned_cluster_bams_tar
    File? final_cluster_exprs_tar = cluster_guided.LRAA_final_cluster_exprs_tar
    File? final_cluster_trackings_tar = cluster_guided.LRAA_final_cluster_trackings_tar
    
    # Preliminary intermediate outputs from cluster-guided discovery phase
    File? prelim_cluster_gtfs = cluster_guided.LRAA_prelim_cluster_gtfs
    File? prelim_cluster_read_trackings = cluster_guided.LRAA_prelim_cluster_read_trackings
    File? prelim_cluster_pseudobulk_exprs = cluster_guided.LRAA_prelim_cluster_pseudobulk_exprs
    
    File? cluster_gene_counts_matrix = cluster_guided.cluster_gene_counts_matrix
    File? cluster_gene_TPM_matrix = cluster_guided.cluster_gene_TPM_matrix
    File? cluster_isoform_counts_matrix = cluster_guided.cluster_isoform_counts_matrix
    File? cluster_isoform_TPM_matrix = cluster_guided.cluster_isoform_TPM_matrix
    File? cluster_isoform_counts_forDiffIsoUsage = cluster_guided.cluster_isoform_counts_forDiffIsoUsage

    File? incl_gene_symbols_gffcompare_tracking = add_gene_symbols.gffcompare_tracking
    File? incl_gene_symbols_gffcompare_stats = add_gene_symbols.gffcompare_stats
    File? incl_gene_symbols_updated_gtf = add_gene_symbols.updated_gtf_with_gene_symbols
    File? incl_gene_symbols_updated_id_mappings = add_gene_symbols.updated_id_mappings
    File? incl_gene_symbols_gene_sparse_tar_gz = add_gene_symbols.updated_gene_sparse_tar_gz
    File? incl_gene_symbols_isoform_sparse_tar_gz = add_gene_symbols.updated_isoform_sparse_tar_gz
    File? incl_gene_symbols_splice_pattern_sparse_tar_gz = add_gene_symbols.updated_splice_pattern_sparse_tar_gz
  }
}
