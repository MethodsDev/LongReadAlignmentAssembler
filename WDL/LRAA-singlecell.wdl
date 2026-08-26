version 1.0

##########################################################################################
# LRAA Single-Cell RNA-seq Workflow
##########################################################################################
#
# This workflow performs isoform discovery and/or quantification on single-cell long-read
# RNA-seq data with two primary execution modes and multiple optimization paths.
#
# ===== EXECUTION MODES =====
#
# 1. BASIC MODE (single_cell_pipe_mode = "basic")
#    - Runs initial discovery on the full BAM to generate a preliminary GTF
#    - Builds single-cell sparse matrices (gene, isoform, splice pattern) from read tracking
#    - Optionally filters good cells from empty droplets
#    - Clusters cells using Seurat on gene-level expression
#    - Outputs: initial GTF, single-cell matrices, cluster assignments, UMAP visualizations
#    - Use when: You want a quick initial analysis without per-cluster refinement
#
# 2. CLUSTER-GUIDED MODE (single_cell_pipe_mode = "cluster-guided")
#    - Performs all basic mode steps PLUS:
#    - Partitions BAM by cell clusters
#    - Re-runs LRAA discovery independently on each cluster's reads
#    - Merges cluster-specific GTFs into a refined final GTF
#    - Re-quantifies all cells against the final GTF
#    - Builds final single-cell sparse matrices with improved isoform definitions
#    - Outputs: refined final GTF, improved single-cell matrices, per-cluster pseudobulk data
#    - Use when: You want the highest quality isoform definitions tailored to cell types
#
# ===== QUANTIFICATION-ONLY MODE =====
#
# Set quant_only = true to skip isoform discovery and quantify against a known annotation:
#    - In BASIC mode with quant_only = true:
#      * Requires: initial_annot_gtf (the annotation to quantify against)
#      * Skips: Initial discovery phase
#      * Runs: Quantification, matrix building, cell filtering, clustering
#      * Outputs: Single-cell matrices and clusters based on the provided annotation
#
#    - In CLUSTER-GUIDED mode with quant_only = true:
#      * Requires: initial_annot_gtf
#      * Skips: Both initial discovery AND per-cluster discovery
#      * Runs: Clustering (if not precomputed), per-cluster quantification only
#      * Outputs: Single-cell matrices from cluster-guided quantification
#
# ===== PRECOMPUTED INPUTS (OPTIMIZATION PATHS) =====
#
# You can skip expensive computation steps by providing precomputed files from prior runs:
#
# 1. Skip Initial Discovery:
#    - Provide: precomputed_init_gtf + precomputed_init_quant_tracking
#    - Effect: Skips LRAA_init, uses your GTF and tracking files directly
#    - Use case: Re-run clustering with different parameters without re-discovering isoforms
#
# 2. Skip Initial Discovery + Matrix Building + Clustering:
#    - Provide: precomputed_init_gtf + precomputed_cluster_assignments_tsv
#    - Effect: Jumps directly to cluster-guided refinement (if mode = "cluster-guided")
#             or just outputs the precomputed results (if mode = "basic")
#    - Use case: Test cluster-guided refinement with different discovery parameters
#
# 3. Kickstart Cluster-Guided Mode:
#    - Provide: precomputed_init_gtf + precomputed_cluster_assignments_tsv
#    - Set: single_cell_pipe_mode = "cluster-guided"
#    - Flow: Uses precomputed GTF as the annotation for cluster-guided discovery
#            (passed as annot_gtf to guide per-cluster LRAA runs)
#    - Use case: You have initial results and want refined cluster-specific isoforms
#               without re-running the expensive initial discovery
#
# 4. Cluster-Guided Quantification-Only:
#    - Provide: precomputed_cluster_assignments_tsv + initial_annot_gtf
#    - Set: single_cell_pipe_mode = "cluster-guided", quant_only = true
#    - Flow: Skips all discovery, quantifies each cluster against initial_annot_gtf
#    - Use case: Quantify known isoforms with cluster-aware quantification
#
# ===== GENE SYMBOL INCORPORATION =====
#
# Optionally provide ref_annot_gtf_source_gene_symbols to incorporate gene symbols:
#    - Runs gffcompare between LRAA GTF and reference annotation
#    - Updates all single-cell matrices with gene symbols from matched loci
#    - Only runs when: All required sparse matrices are available (from either basic or
#                      cluster-guided mode, and not skipped due to missing precomputed inputs)
#
# ===== TYPICAL WORKFLOWS =====
#
# Quick exploration:
#   single_cell_pipe_mode = "basic"
#   → Get initial GTF, matrices, and clusters quickly
#
# High-quality production run:
#   single_cell_pipe_mode = "cluster-guided"
#   → Full pipeline with cluster-specific isoform refinement
#
# Iterative refinement:
#   1st run: mode = "basic" (save init GTF and cluster assignments)
#   2nd run: mode = "cluster-guided" + precomputed files
#   → Skip expensive initial phase, focus on cluster-specific refinement
#
# Quantify known annotation:
#   quant_only = true + initial_annot_gtf
#   → Fast quantification without discovery
#
##########################################################################################

import "LRAA.wdl" as LRAA
import "subwdls/LRAA-build_sparse_matrices_from_tracking.wdl" as BuildMatrices
import "subwdls/LRAA-filter_good_cells.wdl" as FilterCells
import "subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl" as Seurat
import "subwdls/Incorporate_gene_symbols.wdl" as GeneSymbols
import "LRAA-cell_cluster_guided.wdl" as ClusterGuided
# Imported for ONE task, emit_shared_chunk_plan: the run-wide cut plan is selected
# here, before the initial pass, and threaded into every phase. The task lives beside
# LRAA_quant_by_cluster's comparability argument rather than being duplicated here.
import "LRAA_quant_by_cluster.wdl" as QuantByCluster

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
    Boolean rescue_unassigned_reads_via_transcriptome_alignment = true

    # Chunk geometry. Approximate MEGABASES between cut points, and the TOTAL width
    # of the search window centred on each target cut. LRAA's defaults are 10 and 1,
    # so a contig shorter than 10 Mb is never cut and chunking degenerates to one
    # chunk holding the whole thing -- which is what a small test fixture does,
    # leaving each run only its two strand units to parallelise over.
    #
    # The init and the per-cluster runs would each want their own geometry, and no
    # longer CAN have one: ONE cut plan is emitted below, before any phase runs, and
    # every phase applies it, so a run has exactly one set of cut positions. Two
    # geometries would mean a phase whose chunk bounds differ from the plan's, which
    # ChunkedRun refuses by name and value rather than quietly cutting elsewhere.
    #
    # So the *_init pair sets the geometry for the WHOLE run when given, and the plain
    # pair is the fallback -- the same precedence the init call used to apply to itself.
    # The init wins because it is the phase that cannot be parallelised any other way:
    # it is ONE task over every read, and on the 2 Mb test fixture 10 chunks took its
    # slowest unit from 366 s to 36 s. The per-cluster phases already run side by side,
    # so finer cuts cost them wall time and nothing else -- measured on the same
    # fixture, cutting the clusters into 10 took their stage-5 work from 27.2 to 43.7
    # min (discovery) and 14.1 to 24.8 min (final quant), because that phase was
    # already saturating every core. An under-cut init is a serial floor under the
    # whole run; over-cut clusters are a slower run that is still correct.
    #
    # Both pairs unset is the production case and is unchanged: LRAA's own 10 Mb / 1 Mb
    # defaults then apply to every phase, as they did before one plan governed the run.
    Float? approx_MB_per_cut
    Float? approx_MB_per_cut_wiggle_window
    Float? approx_MB_per_cut_init
    Float? approx_MB_per_cut_wiggle_window_init

    # How each LRAA phase divides its work; see LRAA.wdl's `scattering`. The three
    # phases want different shapes. The initial pass is ONE task over every read,
    # so by_chunk is the only way to parallelise it. The per-cluster runs are many
    # small tasks that already run side by side, and cutting them further only
    # multiplies per-unit fixed cost -- measured on the test fixture, giving them
    # 10 chunks took their stage-5 work from 27.2 to 43.7 min for no gain, because
    # that phase already saturated every core.
    String scattering_init = "by_chunk"
    String scattering_per_cluster = "by_chromosome"
    String scattering_final_quant = "by_chromosome"

    # Optional: reuse outputs from a prior initial discovery run and skip LRAA_init
    File? precomputed_init_quant_tracking
    File? precomputed_init_gtf
    
    # Optional: provide precomputed cluster assignments to skip steps 1-3
    File? precomputed_cluster_assignments_tsv

    # Single-cell barcode and UMI tags
    String cell_barcode_tag = "CB"
    String read_umi_tag = "XM"

    # Resources and docker (propagated to subcalls where applicable)
    # Cores per LRAA task: the task's cpu request AND the --cpu_budget it divides across
    # work units. One value each for direct and chromosome-sharded runs; there is no
    # second knob to multiply it by.
    Int cpu = 5
    # Used only when main_chromosomes is non-empty and LRAA runs are chromosome-sharded.
    # OPTIONAL, deliberately. LRAA.wdl resolves this as
    # select_first([cpuScattered, shard_cpu_computed]), so any concrete value here
    # wins unconditionally and its contig-length-derived per-shard sizing never
    # runs -- which is what a default of 5 did on every scattered single-cell run.
    # Unset lets that estimate apply; set it to pin a value.
    Int? cpuScattered
    # The initial whole-library LRAA run, separately from the per-cluster runs. It
    # is ONE task over every read, so its parallelism is bounded by how many chunks
    # it was cut into, while the per-cluster runs are many small tasks whose
    # parallelism comes from running side by side. Those want opposite settings: the
    # init wants a wide budget it can spend across chunks, the cluster tasks want to
    # be narrow enough that many fit at once. Unset means "same as cpu/memoryGB",
    # so nothing changes unless asked.
    Int? cpuInit
    Int? memoryGBInit
    # Optional override for direct LRAA runs. When unset, LRAA.wdl uses max(64 GiB, ceil(1.5 x input BAM GiB)).
    Int? memoryGB
    # Optional override for chromosome-sharded LRAA workers only. Has no effect when main_chromosomes is empty.
    Int? memoryGBPerWorkerScattered
    # Optional by_chunk overrides: make-chunks, chunk leaves, and merge need
    # independently sized resources. Unset preserves LRAA.wdl's defaults.
    Int? chunkMakeChunksCpu
    Int? chunkMakeChunksMemoryGB
    Int? chunkCpu
    Int? chunkMemoryGB
    Int? chunkMergeCpu
    Int? chunkMergeMemoryGB
    # The one-off run-wide chunk-plan task. Cut selection runs one concurrent unit per
    # contig under --cpu_budget, so cpu IS that budget; the memory number is justified
    # where these two inputs are declared in LRAA_quant_by_cluster.wdl. Forwarded to
    # the cluster-guided phase as well, so ONE pair of numbers covers the task
    # wherever a run ends up emitting it.
    Int cpu_chunk_plan = 16
    Int memoryGB_chunk_plan = 8
    Int memoryGBbuildSparseMatrices = 32
    Int memoryGBFilterCells = 32
    Int memoryGBSeurat = 32
    Int memoryGBmergeGTFs = 32
    Int memoryGBquantFinal = 32
    Int memoryGBscSparseMatrices = 32
    String sparseMatrixCsvEngine = "python"
    Int sparseMatrixGzipLevel = 1
    Int diskSizeGB = 256
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    String docker_sc = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"

    # Seurat clustering parameters (forwarded to Seurat subworkflow)
    Int min_cells = 10
    Int min_features = 1000
    Float percent_mt_max = 20.0
    String mt_pattern = "^(MT-|mt-|g:(chrM|MT|M):)"
    Int npcs = 12
    Float resolution = 0.6
    Int n_variable_features = 2000
    Int seed = 1

    # Filter good cells parameters
    Boolean enable_filter_good_cells = true
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

  # ONE geometry for the whole run; see the input comments above for the precedence and
  # for why there cannot be two.
  Float? approx_MB_per_cut_run = if defined(approx_MB_per_cut_init) then approx_MB_per_cut_init else approx_MB_per_cut
  Float? approx_MB_per_cut_wiggle_window_run = if defined(approx_MB_per_cut_wiggle_window_init) then approx_MB_per_cut_wiggle_window_init else approx_MB_per_cut_wiggle_window

  # ONE cut plan for the WHOLE cluster-guided run, emitted BEFORE the initial pass.
  #
  # Emitted HERE because this is the earliest point in the pipeline: a plan has to
  # precede phase 1 to be of any use to it, and all three phases -- initial pass,
  # per-cluster discovery, final quant -- have to apply the same one. Emitted inside
  # the final quant instead, which is where it used to live, it served that phase alone
  # and per-cluster discovery still selected its own cuts: 76 distinct geometries over
  # 270 extractions on the ref-guided smoke run.
  #
  # Selected against the REFERENCE annotation, the only one that exists before any
  # phase has run, and safe BECAUSE it is that early rather than in spite of it: cut
  # selection will not place a boundary inside an annotated model, and discovery then
  # runs INSIDE these chunks, so no novel isoform spanning a cut can be constructed --
  # no chunk ever held the reads that would support one. Measured on the two smoke
  # runs: 0 of 994 consolidated ref-guided models and 0 of 865 de novo models straddle
  # any of the 9 cut positions. De novo has no reference annotation at all and gets an
  # unannotated plan: nothing was protected from severing at selection time, and each
  # consumer checks its own models against the plan's cuts where it applies it.
  #
  # Emitted only when the CLUSTER-GUIDED phases run. Those are the many-sibling phases
  # -- one job per cluster -- and shared geometry is what makes their results
  # comparable. In basic mode the initial pass is the only LRAA run in the pipeline, so
  # there is nothing for its cuts to agree with, and a plan would cost rather than buy:
  # selection would move out of that run's own make_chunks, which OVERLAPS it with
  # extraction (pylib/ChunkedRun.py's prep pool submits a contig's extractions as soon
  # as that contig's selection finishes), into a separate task the run then waits on.
  Boolean want_chunk_plan = run_cluster_guided
      && (scattering_per_cluster != "off" || scattering_final_quant != "off")

  if (want_chunk_plan) {
    call QuantByCluster.emit_shared_chunk_plan as emit_run_chunk_plan {
      input:
        inputBAM = inputBAM,
        referenceGenome = referenceGenome,
        # The most complete annotation that exists BEFORE any phase runs. With the
        # initial pass skipped, that is the precomputed init GTF -- the models the
        # cluster phases then quantify -- rather than the reference, so placement
        # protects those too.
        annot_gtf = if run_initial_phase then initial_annot_gtf else precomputed_init_gtf,
        HiFi = HiFi,
        approx_MB_per_cut = approx_MB_per_cut_run,
        approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window_run,
        cpu = cpu_chunk_plan,
        memoryGB = memoryGB_chunk_plan,
        docker = docker
    }
  }

  # The run's geometry, threaded to every phase below.
  File? run_chunk_plan = emit_run_chunk_plan.chunk_plan

  # The initial pass JOINS that geometry when it runs and is scattered, and it has to:
  # the models it discovers become the annotation the cluster phases quantify, so a
  # model found across one of the plan's cuts would be severed by every phase that
  # applies the plan. With scattering_init=off it cannot -- LRAA.wdl's
  # validate_scattering refuses a plan for a single whole-genome invocation -- so that
  # configuration holds the property only empirically, and a model that does straddle a
  # cut is refused by the consumer that would have severed it rather than quantified
  # quietly.
  if (want_chunk_plan && run_initial_phase && scattering_init != "off") {
    File init_chunk_plan = select_first([run_chunk_plan])
  }

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
        # single-cell: this workflow never surfaces the normalized splice-graph BAM
        retain_normalized_splice_graph_bam = false,
        rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
        cell_barcode_tag = cell_barcode_tag,
        read_umi_tag = read_umi_tag,
        scattering = scattering_init,
        # The run's ONE geometry, selected above. The initial pass discovers inside
        # these chunks, which is what makes the same cuts safe for the cluster phases
        # that quantify what it found.
        internal_chunk_plan = init_chunk_plan,
        approx_MB_per_cut = approx_MB_per_cut_run,
        approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window_run,
        cpu = select_first([cpuInit, cpu]),
        cpuScattered = cpuScattered,
        memoryGB = if defined(memoryGBInit) then memoryGBInit else memoryGB,
        memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
        chunkMakeChunksCpu = chunkMakeChunksCpu,
        chunkMakeChunksMemoryGB = chunkMakeChunksMemoryGB,
        chunkCpu = chunkCpu,
        chunkMemoryGB = chunkMemoryGB,
        chunkMergeCpu = chunkMergeCpu,
        chunkMergeMemoryGB = chunkMergeMemoryGB,
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
        docker = docker_sc,
        memoryGB = memoryGBbuildSparseMatrices,
        csv_engine = sparseMatrixCsvEngine,
        gzip_level = sparseMatrixGzipLevel
    }

    # 2.5) Filter good cells from the gene-level sparse matrix (optional)
    if (enable_filter_good_cells) {
      call FilterCells.FilterGoodCells as filter_good_cells {
        input:
          sample_id = sample_id,
          gene_sparse_tar_gz = build_sc_from_init_tracking.gene_sparse_dir_tgz,
          isoform_sparse_tar_gz = build_sc_from_init_tracking.isoform_sparse_dir_tgz,
          splice_pattern_sparse_tar_gz = build_sc_from_init_tracking.splice_pattern_sparse_dir_tgz,
          docker = docker_sc,
          memoryGB = memoryGBFilterCells,
          fdr_threshold = fdr_threshold,
          lower_threshold = lower_threshold
      }
    }

    # 3) Cluster cells from the gene-level sparse matrix (filtered or unfiltered)
    File gene_sparse_for_clustering = if enable_filter_good_cells then select_first([filter_good_cells.filtered_gene_sparse_tar_gz]) else build_sc_from_init_tracking.gene_sparse_dir_tgz
    
    call Seurat.GeneSparseM_To_SeuratClusters as cluster_cells {
      input:
        sample_id = sample_id,
        gene_sparse_tar_gz = gene_sparse_for_clustering,
        docker = docker_sc,
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
        scattering = scattering_per_cluster,
        scattering_final_quant = scattering_final_quant,
        # The run's ONE geometry, emitted above before the initial pass. Passed
        # whenever it exists; that workflow gates it per phase, since a phase running
        # scattering=off takes no plan.
        internal_chunk_plan = run_chunk_plan,
        cpu_chunk_plan = cpu_chunk_plan,
        memoryGB_chunk_plan = memoryGB_chunk_plan,
        approx_MB_per_cut = approx_MB_per_cut_run,
        approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window_run,
        cpu = cpu,
        cpuScattered = cpuScattered,
        memoryGB = memoryGB,
        memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
        chunkMakeChunksCpu = chunkMakeChunksCpu,
        chunkMakeChunksMemoryGB = chunkMakeChunksMemoryGB,
        chunkCpu = chunkCpu,
        chunkMemoryGB = chunkMemoryGB,
        chunkMergeCpu = chunkMergeCpu,
        chunkMergeMemoryGB = chunkMergeMemoryGB,
        memoryGBmergeGTFs = memoryGBmergeGTFs,
        memoryGBquantFinal = memoryGBquantFinal,
        memoryGBscSparseMatrices = memoryGBscSparseMatrices,
        sparseMatrixCsvEngine = sparseMatrixCsvEngine,
        sparseMatrixGzipLevel = sparseMatrixGzipLevel,
        diskSizeGB = diskSizeGB,
        docker = docker,
        docker_sc = docker_sc,
        quant_only_cluster_guided = quant_only
    }
  }

  # 5) Incorporate gene symbols: use cluster-guided outputs if available, otherwise use filtered (if enabled) or unfiltered good cell outputs
  # In quant_only mode, these GTF values will naturally be undefined since discovery doesn't produce them
  File? gtf_for_symbols = if run_cluster_guided then cluster_guided.LRAA_final_gtf else init_gtf_file
  File? gene_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_gene_sparse_tar_gz else (if enable_filter_good_cells then filter_good_cells.filtered_gene_sparse_tar_gz else build_sc_from_init_tracking.gene_sparse_dir_tgz)
  File? isoform_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_isoform_sparse_tar_gz else (if enable_filter_good_cells then filter_good_cells.filtered_isoform_sparse_tar_gz else build_sc_from_init_tracking.isoform_sparse_dir_tgz)
  File? splice_pattern_sparse_for_symbols = if run_cluster_guided then cluster_guided.sc_splice_pattern_sparse_tar_gz else (if enable_filter_good_cells then filter_good_cells.filtered_splice_pattern_sparse_tar_gz else build_sc_from_init_tracking.splice_pattern_sparse_dir_tgz)
  File? mapping_for_symbols = if run_cluster_guided then cluster_guided.sc_gene_transcript_splicehash_mapping else build_sc_from_init_tracking.mapping_file

  if (defined(ref_annot_gtf_source_gene_symbols) && defined(gene_sparse_for_symbols) && defined(isoform_sparse_for_symbols) && defined(splice_pattern_sparse_for_symbols) && defined(mapping_for_symbols)) {
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
    # The ONE cut geometry every phase of a cluster-guided run applied. Absent in basic
    # mode, which has a single LRAA run and so nothing to be consistent with, and when
    # every cluster-guided phase is scattering=off.
    File? shared_chunk_plan = run_chunk_plan

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
