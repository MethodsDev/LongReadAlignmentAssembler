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

    # Platform/options
    Boolean HiFi = false
    String oversimplify = "chrM"   # e.g., "chrM" or "chrM,M"
    String main_chromosomes = ""  # if empty, runs without partitioning
    String? region                 # e.g., "chr1:100000-200000"; forces direct mode

    # Resources and docker (propagated to subcalls where applicable)
    Int numThreadsPerWorker = 5
    Int numThreadsPerWorkerScattered = 5
    Int num_parallel_contigs = 3
    Int numThreadsPerLRAA = 4
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

  # 1) Initial transcript discovery + quantification on the full BAM
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

  # 2) Build single-cell sparse matrices from the initial tracking
  call BuildMatrices.BuildSparseMatricesFromTracking as build_sc_from_init_tracking {
    input:
      sample_id = sample_id,
      tracking_file = LRAA_init.mergedQuantTracking,
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

  # 4) Cluster-guided reconstruction + final single-cell matrices
  call ClusterGuided.LRAA_cell_cluster_guided as cluster_guided {
    input:
      sample_id = sample_id,
      referenceGenome = referenceGenome,
      inputBAM = inputBAM,
      cell_clusters_info = cluster_cells.cluster_assignments_tsv,
      annot_gtf = LRAA_init.mergedGTF,
      HiFi = HiFi,
      oversimplify = oversimplify,
      main_chromosomes = main_chromosomes,
      numThreadsPerWorker = numThreadsPerWorker,
      numThreadsPerWorkerScattered = numThreadsPerWorkerScattered,
      num_parallel_contigs = num_parallel_contigs,
      numThreadsPerLRAA = numThreadsPerLRAA,
      memoryGB = memoryGB,
      memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
      memoryGBmergeGTFs = memoryGBmergeGTFs,
      memoryGBquantFinal = memoryGBquantFinal,
      memoryGBscSparseMatrices = memoryGBscSparseMatrices,
      diskSizeGB = diskSizeGB,
      docker = docker,
      quant_only_cluster_guided = false
  }

  if (defined(cluster_guided.LRAA_final_gtf) && defined(initial_annot_gtf)) {
    call GeneSymbols.Incorporate_gene_symbols as add_gene_symbols {
      input:
        sample_id = sample_id,
        reference_gtf = select_first([initial_annot_gtf]),
        final_gtf = select_first([cluster_guided.LRAA_final_gtf]),
        final_sc_gene_sparse_tar_gz = cluster_guided.sc_gene_sparse_tar_gz,
        final_sc_isoform_sparse_tar_gz = cluster_guided.sc_isoform_sparse_tar_gz,
        final_sc_splice_pattern_sparse_tar_gz = cluster_guided.sc_splice_pattern_sparse_tar_gz,
        final_sc_gene_transcript_splicehash_mapping = cluster_guided.sc_gene_transcript_splicehash_mapping,
        docker = docker
    }
  }

  output {
    # Initial discovery outputs
    File init_quant_expr = LRAA_init.mergedQuantExpr
    File init_quant_tracking = LRAA_init.mergedQuantTracking
    File? init_gtf = LRAA_init.mergedGTF

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
    File final_tracking = cluster_guided.LRAA_final_tracking
    File final_sc_gene_sparse_tar_gz = cluster_guided.sc_gene_sparse_tar_gz
    File final_sc_isoform_sparse_tar_gz = cluster_guided.sc_isoform_sparse_tar_gz
    File final_sc_splice_pattern_sparse_tar_gz = cluster_guided.sc_splice_pattern_sparse_tar_gz
    File final_sc_gene_transcript_splicehash_mapping = cluster_guided.sc_gene_transcript_splicehash_mapping

    # Convenience tarballs and matrices from the cluster-guided phase
    File partitioned_cluster_bams_tar = cluster_guided.LRAA_partitioned_cluster_bams_tar
    File final_cluster_exprs_tar = cluster_guided.LRAA_final_cluster_exprs_tar
    File final_cluster_trackings_tar = cluster_guided.LRAA_final_cluster_trackings_tar
    File cluster_gene_counts_matrix = cluster_guided.cluster_gene_counts_matrix
    File cluster_gene_TPM_matrix = cluster_guided.cluster_gene_TPM_matrix
    File cluster_isoform_counts_matrix = cluster_guided.cluster_isoform_counts_matrix
    File cluster_isoform_TPM_matrix = cluster_guided.cluster_isoform_TPM_matrix
    File cluster_isoform_counts_forDiffIsoUsage = cluster_guided.cluster_isoform_counts_forDiffIsoUsage

    File? incl_gene_symbols_gffcompare_tracking = add_gene_symbols.gffcompare_tracking
    File? incl_gene_symbols_gffcompare_stats = add_gene_symbols.gffcompare_stats
    File? incl_gene_symbols_updated_gtf = add_gene_symbols.updated_gtf_with_gene_symbols
    File? incl_gene_symbols_updated_id_mappings = add_gene_symbols.updated_id_mappings
    File? incl_gene_symbols_gene_sparse_tar_gz = add_gene_symbols.updated_gene_sparse_tar_gz
    File? incl_gene_symbols_isoform_sparse_tar_gz = add_gene_symbols.updated_isoform_sparse_tar_gz
    File? incl_gene_symbols_splice_pattern_sparse_tar_gz = add_gene_symbols.updated_splice_pattern_sparse_tar_gz
  }
}
