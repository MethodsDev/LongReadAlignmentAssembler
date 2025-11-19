version 1.0

import "LRAA.wdl" as LRAA
import "subwdls/LRAA-build_sparse_matrices_from_tracking.wdl" as BuildMatrices
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

    # Platform/options
    Boolean HiFi = false
    String oversimplify = "chrM"   # e.g., "chrM" or "chrM,M"
    String main_chromosomes = ""  # if empty, runs without partitioning
    String? region                 # e.g., "chr1:100000-200000"; forces direct mode

    # Optional: reuse outputs from a prior initial discovery run and skip LRAA_init
    File? precomputed_init_quant_tracking
    File? precomputed_init_gtf

    # Resources and docker (propagated to subcalls where applicable)
    Int numThreadsPerWorker = 5
    Int numThreadsPerWorkerScattered = 5
    Int num_parallel_contigs = 3
    Int memoryGB = 64
    Int memoryGBPerWorkerScattered = 32
    Int memoryGBbuildSparseMatrices = 32
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
  }

  # Only skip the initial discovery call when the downstream-critical artifacts are provided
  Boolean has_precomputed_init = defined(precomputed_init_quant_tracking) && defined(precomputed_init_gtf)
  Boolean run_initial_phase = !has_precomputed_init
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
        numThreadsPerWorker = numThreadsPerWorker,
        numThreadsPerWorkerScattered = numThreadsPerWorkerScattered,
        num_parallel_contigs = num_parallel_contigs,
        memoryGB = memoryGB,
        memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
        diskSizeGB = diskSizeGB,
        docker = docker
    }
  }

  File? init_quant_expr_file = LRAA_init.mergedQuantExpr
  File? init_quant_tracking_generated = LRAA_init.mergedQuantTracking
  File init_quant_tracking_file = if (!run_initial_phase && defined(precomputed_init_quant_tracking)) then select_first([precomputed_init_quant_tracking]) else select_first([init_quant_tracking_generated])
  File? init_gtf_generated = LRAA_init.mergedGTF
  File? init_gtf_file = if (!run_initial_phase && defined(precomputed_init_gtf)) then select_first([precomputed_init_gtf]) else init_gtf_generated

  # 2) Build single-cell sparse matrices from the initial tracking
  call BuildMatrices.BuildSparseMatricesFromTracking as build_sc_from_init_tracking {
    input:
      sample_id = sample_id,
      tracking_file = init_quant_tracking_file,
      docker = docker,
      memoryGB = memoryGBbuildSparseMatrices
  }

  # 3) Cluster cells from the gene-level sparse matrix
  call Seurat.GeneSparseM_To_SeuratClusters as cluster_cells {
    input:
      sample_id = sample_id,
      gene_sparse_tar_gz = build_sc_from_init_tracking.gene_sparse_dir_tgz,
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

  # 4) Cluster-guided reconstruction + final single-cell matrices (only if mode is 'cluster-guided')
  if (run_cluster_guided) {
    call ClusterGuided.LRAA_cell_cluster_guided as cluster_guided {
      input:
        sample_id = sample_id,
        referenceGenome = referenceGenome,
        inputBAM = inputBAM,
        cell_clusters_info = cluster_cells.cluster_assignments_tsv,
        annot_gtf = init_gtf_file,
        HiFi = HiFi,
        oversimplify = oversimplify,
        main_chromosomes = main_chromosomes,
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
        quant_only_cluster_guided = false
    }
  }

  # 5) Incorporate gene symbols: use cluster-guided outputs if available, otherwise use initial outputs
  if (defined(initial_annot_gtf)) {
    File? gtf_for_symbols = if run_cluster_guided then cluster_guided.LRAA_final_gtf else init_gtf_file
    File? gene_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_gene_sparse_tar_gz else build_sc_from_init_tracking.gene_sparse_dir_tgz
    File? isoform_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_isoform_sparse_tar_gz else build_sc_from_init_tracking.isoform_sparse_dir_tgz
    File? splice_pattern_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_splice_pattern_sparse_tar_gz else build_sc_from_init_tracking.splice_pattern_sparse_dir_tgz
    File? mapping_for_symbols = if run_cluster_guided then cluster_guided.sc_gene_transcript_splicehash_mapping else build_sc_from_init_tracking.mapping_file

    if (defined(gtf_for_symbols)) {
      call GeneSymbols.Incorporate_gene_symbols as add_gene_symbols {
        input:
          sample_id = sample_id,
          reference_gtf = select_first([initial_annot_gtf]),
          final_gtf = select_first([gtf_for_symbols]),
          final_sc_gene_sparse_tar_gz = select_first([gene_sparse_for_symbols]),
          final_sc_isoform_sparse_tar_gz = select_first([isoform_sparse_for_symbols]),
          final_sc_splice_pattern_sparse_tar_gz = select_first([splice_pattern_sparse_for_symbols]),
          final_sc_gene_transcript_splicehash_mapping = select_first([mapping_for_symbols]),
          docker = docker
      }
    }
  }

  output {
    # Initial discovery outputs
    File? init_quant_expr = init_quant_expr_file
    File? init_quant_tracking = init_quant_tracking_generated
    File? init_gtf = init_gtf_generated

    # Initial single-cell matrices and clustering inputs/outputs
    File init_sc_gene_sparse_tar_gz = build_sc_from_init_tracking.gene_sparse_dir_tgz
    File init_sc_isoform_sparse_tar_gz = build_sc_from_init_tracking.isoform_sparse_dir_tgz
    File init_sc_splice_pattern_sparse_tar_gz = build_sc_from_init_tracking.splice_pattern_sparse_dir_tgz
    File init_sc_gene_transcript_splicehash_mapping = build_sc_from_init_tracking.mapping_file
    File seurat_umap_pdf = cluster_cells.umap_pdf
    File seurat_cluster_assignments = cluster_cells.cluster_assignments_tsv

    # Final cluster-guided outputs (main deliverables)
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
