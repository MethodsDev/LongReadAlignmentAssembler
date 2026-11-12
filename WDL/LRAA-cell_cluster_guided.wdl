version 1.0

import "LRAA.wdl" as LRAA
import "LRAA_quant_by_cluster.wdl" as LRAA_quant_by_cluster
import "subwdls/partition_bam_by_cell_cluster.wdl" as PartitionBam
import "subwdls/LRAA-build_sparse_matrices_from_tracking.wdl" as BuildMatrices



workflow LRAA_cell_cluster_guided {

    input {
        String sample_id
        File referenceGenome

        File inputBAM
        File cell_clusters_info

        File? annot_gtf

        # The pooled (all-cells) LRAA gtf -- in LRAA-singlecell.wdl the initial pass's
        # catalog, the same file passed as annot_gtf there. Merged in alongside the
        # per-cluster gtfs (merge_LRAA_GTFs.py --pooled_gtf) so structures that need the
        # pooled read depth are kept, its TSS/PolyA support counted as a residual over the
        # cluster support (about max, not sum) so boundaries are not inflated. Unset =
        # merge the per-cluster gtfs only. Not the reference annotation: every model in it
        # is carried into the final catalog.
        File? pooled_gtf

        # --polyA_known: trusted cleavage-site BED/GTF (v0.43.0) forwarded to the
        # per-cluster discovery runs so they leverage the same trusted 3' ends as the
        # init pass (see LRAA-singlecell.wdl, which threads it here from init). Default
        # unset => the internal-priming veto applies to every A-rich read candidate.
        File? polyA_known

        # INTERNAL PLUMBING, not a user knob: the ONE cut plan for the WHOLE run,
        # emitted by a caller that starts even earlier than this workflow does.
        # LRAA-singlecell.wdl emits it before ITS initial pass and passes it here, so
        # the initial pass, the per-cluster discovery below and the final quant all cut
        # at identical positions. Unset -- this workflow driven directly -- and the
        # emission below produces it instead.
        File? internal_chunk_plan

        Boolean HiFi = false
        # Merge-only: ignore the TSS/PolyA annotations on the per-cluster gtfs when
        # merging them (the former non-HiFi merge). Default false = respect them, since
        # LRAA models TSS/PolyA in both modes; independent of the HiFi run flag above.
        Boolean ignore_TSS_POLYA = false
        # Chunk geometry, forwarded to every per-cluster LRAA run. Unset means
        # LRAA's own defaults (10 Mb spacing, 1 Mb window), under which a contig
        # shorter than 10 Mb is never cut.
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        # `scattering` is per-cluster DISCOVERY, `scattering_final_quant` the final
        # pooled quant. These MUST agree with LRAA-singlecell.wdl's defaults for
        # the same two knobs: that workflow passes both through when it is the
        # entrypoint, so a divergence here would only show when this workflow is
        # called DIRECTLY -- the same run shape silently differing by entry path.
        # Rationale for the split, including the measurement behind keeping discovery
        # on by_chromosome, is at LRAA-singlecell.wdl:167-186.
        String scattering = "by_chromosome"
        String scattering_final_quant = "by_chromosome"
        String? oversimplify # comma-separated contig names to simplify (e.g., "chrM" or "chrM,MT")
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true

        String main_chromosomes = "" # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

        # Contigs each per-cluster partition extracts concurrently. 1 locally, where
        # these jobs share one host; 4 on Terra, where each is its own VM and the
        # 2.19x measured for the fan-out is available on every one of them.
        Int cluster_partition_workers = 1
        # Contigs the final pooled quant's ONE full-library partition extracts
        # concurrently. LRAA-singlecell.wdl passes its initial_partition_workers, since
        # that partition has the same shape as the init pass's; see the measurement there.
        Int final_quant_partition_workers = 1
        
        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
        
        # Cores per LRAA task: the task's cpu request AND the --cpu_budget it divides
        # across work units. There is no second knob to multiply it by.
        Int cpu = 2
        # Used only when main_chromosomes is non-empty and inner LRAA.wdl calls are chromosome-sharded.
        # OPTIONAL: LRAA.wdl resolves select_first([cpuScattered, shard_cpu_computed]),
        # so a concrete value here silently disables its per-shard estimate.
        Int? cpuScattered
        # Optional override for direct inner LRAA.wdl runs.
        Int? memoryGB
        # Optional override for chromosome-sharded inner LRAA.wdl workers only.
        # The outer per-cluster scatter does not use this unless main_chromosomes triggers chromosome sharding inside each cluster run.
        Int? memoryGBPerWorkerScattered
        # Optional by_chunk overrides: make-chunks, chunk leaves, and merge need
        # independently sized resources. Unset preserves LRAA.wdl's defaults.
        Int? chunkMakeChunksCpu
        Int? chunkMakeChunksMemoryGB
        Int? chunkCpu
        Int? chunkMemoryGB
        Int? chunkMergeCpu
        Int? chunkMergeMemoryGB
        Int memoryGBmergeGTFs = 64
        # merge_LRAA_GTFs.py merges each contig/strand in its own worker process.
        # --CPU auto uses the cores the task is granted, capped by its memory:
        # peak ~ parsed inputs + workers x ~5 GB on human chromosomes, and auto
        # budgets 6 GB per worker within 90% of the granted memory. MEASURED, 14
        # whole-genome PBMC cluster gtfs: 4.9 GB parsed, 25.6 GB peak at 4 workers.
        # At 64 GB and 8 cpus that is ~8 workers. Output is identical to serial.
        Int cpu_mergeGTFs = 8
        # The final pooled quant: one LRAA run over the whole library, like the init
        # pass. memoryGB applies only to scattering_final_quant = "off" (one whole-genome
        # invocation); unset, LRAA.wdl sizes it from the bam as for the init pass. Under
        # by_chromosome each shard is sized by LRAA.wdl or memoryGBPerWorkerScattered.
        # cpu unset = `cpu` above.
        Int? memoryGBquantFinal
        Int? cpuFinalQuant
        # The final pooled quant's TPM denominator (RPM_total_reads). Unset -- the
        # normal case -- LRAA counts inputBAM itself, with its -F 0x904 policy. The
        # per-cluster discovery runs never get it: each counts its own cluster bam.
        Int? num_total_reads
        # The one-off shared-chunk-plan task below. Cut selection runs one concurrent
        # unit per contig under --cpu_budget, so cpu IS that budget; the memory number
        # is justified where these same two inputs are declared in
        # LRAA_quant_by_cluster.wdl.
        Int cpu_chunk_plan = 16
        Int memoryGB_chunk_plan = 8

        # Turn off 3' end weighting in the EM. Forwarded to BOTH phases below --
        # per-cluster discovery and the final quant -- because a weighting that
        # differed between the phases would make their numbers incomparable.
        Boolean no_weight_reads_by_3prime_agreement = false
        # LRAA_Globals.config overrides (JSON object) for the per-cluster DISCOVERY runs
        # only; the final quant gets LRAA's defaults. See LRAA-singlecell.wdl.
        String? cluster_config_update_json
        # LRAA_Globals.config overrides (JSON object) for the final pooled quant. Unset by
        # default: it runs on LRAA's defaults. Not set from LRAA-singlecell.wdl.
        String? final_quant_config_update_json
        # LRAA_Globals.config overrides (JSON object) for the cross-cluster merge
        # (merge_LRAA_GTFs.py --config_update). Unset by default. E.g.
        # '{"allow_unflagged_LRAA_guide": true}' to merge GTFs written before v0.50.0,
        # which lack ref_splice_pattern / ref_supported and are otherwise refused.
        String? merge_config_update_json
        # The library's TSS cap calibration (LRAA.wdl capCalibrationJson of the init
        # run), handed to every per-cluster discovery run so none decides from its own
        # sparse reads. Unset, each cluster's runs calibrate on the cluster's bam.
        String? cap_calibration_json
        Int memoryGBscSparseMatrices = 16
        # Build the single-cell matrices per contig shard during the final
        # quantification and merge them, instead of one single-threaded pass
        # over its tracking file. On by
        # default; see LRAA-singlecell.wdl's input for the measured basis.
        # Requires scattering_final_quant = "by_chromosome", the default since
        # v0.32.0; "off" produces no per-contig shards to merge and forces this
        # false below.
        Boolean sc_sparse_from_shards = true
        String sparseMatrixCsvEngine = "direct"
        Int sparseMatrixGzipLevel = 1
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
        String docker_sc = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"
        Boolean quant_only_cluster_guided = false

        # INTERNAL PLUMBING, not a user knob: a whole-library quant the CALLER already
        # ran that is the same run the final quant below would be -- same bam, same gtf,
        # same chunk plan, same settings. LRAA-singlecell.wdl passes its init pass's
        # outputs here in quant_only mode, where the init pass IS a quant-only run of
        # initial_annot_gtf over the full bam. When the tracking and expr are given,
        # LRAA_final_quant does not run; when the four matrices are given too, they are
        # used as the final single-cell matrices instead of being rebuilt.
        File? internal_final_quant_tracking
        File? internal_final_quant_expr
        File? internal_final_quant_read_assignment_summary
        File? internal_final_sc_gene_tgz
        File? internal_final_sc_isoform_tgz
        File? internal_final_sc_splice_pattern_tgz
        File? internal_final_sc_mapping

     }

    # The shard path needs shards; scattering_final_quant = "off" produces none,
    # so the request is forced off rather than left to fail on an empty list.
    Boolean reuse_final_quant = defined(internal_final_quant_tracking) && defined(internal_final_quant_expr)
    Boolean reuse_final_matrices = reuse_final_quant
        && defined(internal_final_sc_gene_tgz) && defined(internal_final_sc_isoform_tgz)
        && defined(internal_final_sc_splice_pattern_tgz) && defined(internal_final_sc_mapping)
    # A reused quant left no shards here, so without its matrices they are rebuilt
    # from its tracking.
    Boolean use_sc_sparse_from_shards =
        sc_sparse_from_shards && scattering_final_quant == "by_chromosome" && !reuse_final_quant

      
     call PartitionBam.partition_bam_by_cell_cluster {
         input:
            sample_id = sample_id,
            cell_clusters_info = cell_clusters_info,
            inputBAM = inputBAM,
            cell_barcode_tag = cell_barcode_tag,
            # Restrict what the cluster bams CONTAIN. Everything downstream of this
            # task on the bam-preparation side -- per-cluster normalization, the
            # merge, and the merged normalization -- consumes these bams, so this one
            # input is what keeps them from doing whole-library work for a run that
            # analyzes a few contigs. inputBAM itself is untouched.
            main_chromosomes = main_chromosomes,
            docker = docker,
     }

    # Package the partitioned cluster BAMs for convenient retrieval
    call LRAA_tar_outputs as tar_partitioned_cluster_bams {
        input:
            tar_directory_name = sample_id + ".LRAA.cluster_partitioned_bams",
            input_files = partition_bam_by_cell_cluster.partitioned_bams,
            docker = docker
    }

    # If running in quant-only mode, require annot_gtf to be provided
    if (quant_only_cluster_guided) {
        call require_annot_gtf {
            input:
                annot_gtf = select_first([annot_gtf]),
                docker = docker
        }
    }

    # ONE cut plan for the WHOLE run, selected BEFORE the phase-1 discovery below.
    #
    # Emitted here rather than inside the final quant, because a plan emitted there
    # serves that phase and nothing earlier: per-cluster discovery then still selects
    # its own cuts, measured at 76 distinct geometries over 270 extractions on the
    # ref-guided smoke run (when the final quant was per cluster). Cut positions are
    # read off a bam, so per-phase selection means per-phase boundaries: the clusters'
    # discovery runs would each cut differently, and the final quant would cut through
    # places no discovery run cut, where a model found inside one chunk can straddle
    # the quant's boundary.
    #
    # Selected against the REFERENCE annotation, the only one that exists this early,
    # and that is correct rather than a compromise. Cut selection will not place a
    # boundary inside an annotated model, so reference genes stay intact; and discovery
    # below runs INSIDE these chunks, so no novel isoform can be constructed spanning a
    # cut -- no chunk ever held the reads that would support one. The models phase 1
    # finds are therefore confined to chunk interiors by construction, and the
    # consolidated gtf the final quant targets cannot be severed by this geometry.
    # Measured on the two smoke runs: 0 of 994 consolidated ref-guided models and 0 of
    # 865 de novo models straddle any of the 9 cut positions -- and those models were
    # discovered under SEPARATE per-phase geometries, which is the harsher test.
    #
    # Do NOT move this back after discovery to select against the consolidated gtf
    # instead. That is what it used to do, and it is exactly what leaves phase 1
    # unplanned: the consolidated gtf does not exist until phase 1 has already chosen
    # its own cuts.
    #
    # Selected on the WHOLE pre-partition bam, the unthinned superset of every cluster
    # bam: a position no read spans in the superset is spanned in no subset, so one
    # plan is safe for all of them.
    Boolean want_chunk_plan = (scattering != "off" || scattering_final_quant != "off")

    if (want_chunk_plan && !defined(internal_chunk_plan)) {
        call LRAA_quant_by_cluster.emit_shared_chunk_plan as emit_run_chunk_plan {
            input:
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                annot_gtf = annot_gtf,
                # Restrict selection to what this run processes; see the task input.
                main_chromosomes = main_chromosomes,
                HiFi = HiFi,
                approx_MB_per_cut = approx_MB_per_cut,
                approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                cpu = cpu_chunk_plan,
                memoryGB = memoryGB_chunk_plan,
                docker = docker
        }
    }

    # The run's geometry: the caller's when it emitted one, else the emission above.
    File? shared_chunk_plan_resolved = if defined(internal_chunk_plan) then internal_chunk_plan else emit_run_chunk_plan.chunk_plan

    # Gated PER PHASE, because a phase running scattering=off takes no plan at all: it
    # is one whole-genome invocation with no sibling run for its geometry to match, and
    # LRAA.wdl's validate_scattering refuses a plan there rather than dropping it. Both
    # branches imply want_chunk_plan, so the plan is defined wherever these are.
    if (scattering != "off") {
        File discovery_chunk_plan = select_first([shared_chunk_plan_resolved])
    }
    if (scattering_final_quant != "off") {
        File final_quant_chunk_plan = select_first([shared_chunk_plan_resolved])
    }


    if (quant_only_cluster_guided == false) {

        scatter (i in range(length(partition_bam_by_cell_cluster.partitioned_bams))) {
            String cluster_sample_id = sub(basename(partition_bam_by_cell_cluster.partitioned_bams[i]), ".bam$", "")
            call LRAA.LRAA_wf as LRAA_by_cluster {
                input:
                    sample_id = cluster_sample_id,
                    referenceGenome = referenceGenome,
                    annot_gtf = annot_gtf,
                    # The run's ONE geometry, selected above BEFORE this discovery ran,
                    # so every model it finds sits inside these chunks and the final
                    # quant can reuse the same cuts without severing one.
                    internal_chunk_plan = discovery_chunk_plan,
                    inputBAM = partition_bam_by_cell_cluster.partitioned_bams[i],
                    cell_list = partition_bam_by_cell_cluster.partitioned_cell_lists[i],
                    HiFi = HiFi,
                    oversimplify = oversimplify,
                    polyA_known = polyA_known,
                    # single-cell: this workflow never surfaces the normalized splice-graph BAM
                    retain_normalized_splice_graph_bam = false,
                    # These per-cluster gtfs are inputs to the cross-cluster merge
                    # below; the collapse runs once on its result instead.
                    emit_splice_pattern_collapse = false,
                    # Discovery only: the single-cell matrices come from the final
                    # pooled quant below.
                    build_sc_sparse_shards = false,
                    docker_sc = docker_sc,
                    rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                    no_weight_reads_by_3prime_agreement = no_weight_reads_by_3prime_agreement,
                    config_update_json = cluster_config_update_json,
                    internal_cap_calibration_json = cap_calibration_json,
                    main_chromosomes = main_chromosomes,
                    partition_workers = cluster_partition_workers,
                    quant_only = false,
                    cell_barcode_tag = cell_barcode_tag,
                    read_umi_tag = read_umi_tag,
                    scattering = scattering,
                    approx_MB_per_cut = approx_MB_per_cut,
                    approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
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
                    # Was DECLARED at :86 and passed nowhere, so raising it bought
                    # nothing and a run could die on a full disk while reporting that
                    # its request had been honoured. miniwdl said so the whole time:
                    # "UnusedDeclaration, nothing references Int diskSizeGB".
                    # LRAA-singlecell.wdl and LRAA.wdl both thread it; this workflow
                    # was the one that did not.
                    #
                    # Reaches BOTH phases: the per-cluster discovery runs here, and
                    # the pooled final quant (LRAA_final_quant) below.
                    diskSizeGB = diskSizeGB,
                    docker = docker
            }

            # Same pairing discipline as LRAA_quant_by_cluster.wdl applies:
            # mergedReadAssignmentSummary is File?, and
            # select_all() compacts, so the cluster id has to be filtered by the SAME
            # predicate to stay aligned with the surviving files.
            if (defined(LRAA_by_cluster.mergedReadAssignmentSummary)) {
                String prelim_cluster_id_for_summary = cluster_sample_id
            }
        }

        # The phase-1 discovery pass's read-assignment accounting, collated into its
        # OWN table. Deliberately not merged with the final quant's: the two passes
        # count the same reads under different normalization and against different
        # annotations, so their totals legitimately disagree and one combined table
        # would invite a false reconciliation.
        Array[File] prelim_read_assignment_summaries = select_all(LRAA_by_cluster.mergedReadAssignmentSummary)
        Array[String] prelim_read_assignment_ids = select_all(prelim_cluster_id_for_summary)

        if (length(prelim_read_assignment_summaries) > 0) {
            call LRAA_quant_by_cluster.collate_read_assignment_summaries as collate_prelim_read_assignment_summaries {
                input:
                    summaryFiles = prelim_read_assignment_summaries,
                    clusterIds = prelim_read_assignment_ids,
                    emitAllClusters = true,
                    outputFilePrefix = sample_id + ".prelim_clusters",
                    population = "cluster_thinned",
                    docker = docker
            }
        }

        # package them up
        call LRAA_tar_outputs as tar_cluster_prelim_gtf_files {
            input:
                tar_directory_name = sample_id + ".LRAA.prelim.cluster_gtfs",
                input_files = select_all(LRAA_by_cluster.mergedGTF),
                docker = docker
        }

        call LRAA_tar_outputs as tar_cluster_prelim_tracking_files {
            input:
                tar_directory_name = sample_id + ".LRAA.prelim.cluster_read_trackings",
                input_files = LRAA_by_cluster.mergedQuantTracking,
                docker = docker
        }

        call LRAA_tar_outputs as tar_cluster_prelim_pseudobulk_expr_files {
            input:
                tar_directory_name = sample_id + ".LRAA.prelim.cluster_pseudobulk.EXPRs",
                input_files = LRAA_by_cluster.mergedQuantExpr,
                docker = docker
        }

        # merge gtfs from the per-cluster runs into a single final gtf.
        call lraa_merge_gtf_task {
            input:
                sample_id = sample_id,
                LRAA_cell_cluster_gtfs = select_all(LRAA_by_cluster.mergedGTF),
                pooled_gtf = pooled_gtf,
                referenceGenome = referenceGenome,
                ignore_TSS_POLYA = ignore_TSS_POLYA,
                # Same list the per-cluster runs were given, so the contigs the
                # merge carries forward are exactly the ones those runs collapsed.
                # select_first because the workflow input is String? (:43) while the
                # task takes String; "" is the task's own "no contigs" value.
                oversimplify = select_first([oversimplify, ""]),
                config_update_json = merge_config_update_json,
                docker = docker,
                memoryGB = memoryGBmergeGTFs,
                cpu = cpu_mergeGTFs,
        }

        # ONE collapse for the whole sample, on the reconciled cross-cluster catalog.
        # Collapsing per cluster instead would group by (gene_id, intron chain) within
        # a cluster, before merge_LRAA_GTFs.py has unified gene assignments across
        # them, so the same splice pattern could land under different gene_ids in
        # different clusters. Same task the bulk workflow uses.
        call LRAA.splice_pattern_collapse as sc_splice_pattern_collapse {
            input:
                lraaGtf = lraa_merge_gtf_task.mergedGTF,
                outputFilePrefix = sample_id + ".LRAA.sc_merged",
                docker = docker
        }
    }


    # Final quantification: ONE quant-only pass over every cell, against the final
    # catalog, in place of the former per-cluster quant (LRAA_quant_by_cluster.wdl).
    #
    # Why pooled. Quantified per cluster, each cluster's EM divides its ambiguous reads
    # by that cluster's own isoform proportions, so every pseudo-replicate of a cluster
    # shares that cluster's prior and differential-usage tests between clusters see
    # differences the prior put there. One EM over all cells gives every cell the same
    # prior: a cell's ambiguous reads are divided by the library-wide proportions, and
    # a difference between clusters has to come from the reads themselves.
    # Cluster-guided DISCOVERY above is unchanged; only the quantification is pooled.
    #
    # The call is the init pass's shape (LRAA-singlecell.wdl LRAA_init) in quant-only
    # mode: the full pre-partition bam, LRAA's own coverage normalization for the
    # splice graph, the run's ONE chunk plan, and the single-cell sparse shards built
    # alongside quantification. Reads from barcodes in no cluster are quantified too:
    # they are in the single-cell matrices and in no cluster's pseudobulk below.
    # Skipped when the caller passed an identical quant to reuse (internal_final_quant_*).
    if (!reuse_final_quant) {
        call LRAA.LRAA_wf as LRAA_final_quant {
            input:
                sample_id = sample_id,
                referenceGenome = referenceGenome,
                inputBAM = inputBAM,
                annot_gtf = select_first([lraa_merge_gtf_task.mergedGTF, annot_gtf]),
                quant_only = true,
                HiFi = HiFi,
                oversimplify = oversimplify,
                # As the init pass gets it, in quant-only mode too.
                polyA_known = polyA_known,
                config_update_json = final_quant_config_update_json,
                main_chromosomes = main_chromosomes,
                # The one full-library partition of this pass, like the init pass's; see
                # final_quant_partition_workers.
                partition_workers = final_quant_partition_workers,
                # Whole library, so the whole library's denominator; unset = LRAA counts
                # inputBAM itself.
                num_total_reads = num_total_reads,
                # single-cell: this workflow never surfaces the normalized splice-graph BAM
                retain_normalized_splice_graph_bam = false,
                # quant-only against the already-collapsed catalog; nothing to collapse
                emit_splice_pattern_collapse = false,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                no_weight_reads_by_3prime_agreement = no_weight_reads_by_3prime_agreement,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                scattering = scattering_final_quant,
                build_sc_sparse_shards = use_sc_sparse_from_shards,
                docker_sc = docker_sc,
                # The SAME geometry phase-1 discovery above ran under, so no model it found
                # can straddle one of this pass's cuts.
                internal_chunk_plan = final_quant_chunk_plan,
                approx_MB_per_cut = approx_MB_per_cut,
                approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                cpu = select_first([cpuFinalQuant, cpu]),
                cpuScattered = cpuScattered,
                memoryGB = memoryGBquantFinal,
                memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
                chunkMakeChunksCpu = chunkMakeChunksCpu,
                chunkMakeChunksMemoryGB = chunkMakeChunksMemoryGB,
                chunkCpu = chunkCpu,
                chunkMemoryGB = chunkMemoryGB,
                chunkMergeCpu = chunkMergeCpu,
                chunkMergeMemoryGB = chunkMergeMemoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker
        }
    }

    File final_quant_tracking = select_first([internal_final_quant_tracking, LRAA_final_quant.mergedQuantTracking])
    File final_quant_expr = select_first([internal_final_quant_expr, LRAA_final_quant.mergedQuantExpr])
    File? final_quant_read_assignment_summary = if reuse_final_quant then internal_final_quant_read_assignment_summary else LRAA_final_quant.mergedReadAssignmentSummary

    # This pass's read-assignment accounting: ONE invocation over the whole bam,
    # normalized once, so the same table shape as the init pass's (LRAA-singlecell.wdl
    # collate_init_read_assignment_summary), not the per-cluster phase's.
    if (defined(final_quant_read_assignment_summary)) {
        call LRAA_quant_by_cluster.collate_read_assignment_summaries as collate_final_read_assignment_summary {
            input:
                summaryFiles = [select_first([final_quant_read_assignment_summary])],
                clusterIds = ["pseudobulk"],
                emitAllClusters = false,
                outputFilePrefix = sample_id + ".final",
                population = "pooled_thinned",
                docker = docker
        }
    }

    # Per-cluster pseudobulk: the pooled assignment summed over each cluster's cells.
    # No EM per cluster; see util/sc/tracking_to_cluster_pseudobulk_exprs.py. File names
    # are the ones the former per-cluster quant runs wrote,
    # <sample_id>.<cluster>.LRAA.quant-only.quant.expr.
    call split_pooled_quant_by_cluster {
        input:
            sample_id = sample_id,
            quant_tracking = final_quant_tracking,
            quant_expr = final_quant_expr,
            cell_clusters_info = cell_clusters_info,
            docker = docker
    }

    # Tar the per-cluster expression files for convenient retrieval
    call LRAA_tar_outputs as tar_final_cluster_expr_files {
        input:
            tar_directory_name = sample_id + ".LRAA.final.cluster_quant.EXPRs",
            input_files = split_pooled_quant_by_cluster.cluster_quant_exprs,
            docker = docker
    }

    # Build cluster-level pseudobulk matrices from the per-cluster expr files
    call build_cluster_pseudobulk_matrices as build_cluster_pseudobulk {
        input:
            sample_id = sample_id,
            quant_expr_files = split_pooled_quant_by_cluster.cluster_quant_exprs,
            docker = docker
    }

     # Build sparse matrices (Seurat-compatible) from the pooled tracking
     if (!reuse_final_matrices && !use_sc_sparse_from_shards) {
         call sc_build_sparse_matrices as build_sc_sparse_matrices {
             input:
                     sample_id = sample_id,
                     tracking_file = final_quant_tracking,
                     docker = docker_sc,
                     memoryGB = memoryGBscSparseMatrices,
                     csv_engine = sparseMatrixCsvEngine,
                     gzip_level = sparseMatrixGzipLevel
         }
     }

     # shared_features stays false, as in basic mode: one LRAA run, so features
     # partition across contigs and a feature turning up in two shards is an error.
     if (!reuse_final_matrices && use_sc_sparse_from_shards) {
         call BuildMatrices.merge_sc_shard_sparse as merge_sc_cluster_shards {
             input:
                 sample_id = sample_id,
                 shard_sparse_tars = select_first([LRAA_final_quant.scShardSparse]),
                 docker = docker_sc,
                 gzip_level = sparseMatrixGzipLevel,
                 shared_features = false
         }
     }

     File sc_gene_tgz = select_first([if reuse_final_matrices then internal_final_sc_gene_tgz else merge_sc_cluster_shards.gene_sparse_dir_tgz,
                                      build_sc_sparse_matrices.gene_sparse_dir_tgz])
     File sc_isoform_tgz = select_first([if reuse_final_matrices then internal_final_sc_isoform_tgz else merge_sc_cluster_shards.isoform_sparse_dir_tgz,
                                         build_sc_sparse_matrices.isoform_sparse_dir_tgz])
     File sc_splice_tgz = select_first([if reuse_final_matrices then internal_final_sc_splice_pattern_tgz else merge_sc_cluster_shards.splice_pattern_sparse_dir_tgz,
                                        build_sc_sparse_matrices.splice_pattern_sparse_dir_tgz])
     File sc_mapping_tsv = select_first([if reuse_final_matrices then internal_final_sc_mapping else merge_sc_cluster_shards.mapping_file,
                                         build_sc_sparse_matrices.mapping_file])

     output {
         # final outputs
         File? LRAA_final_gtf = lraa_merge_gtf_task.mergedGTF
         File? LRAA_final_gtf_tracking = lraa_merge_gtf_task.mergedTracking
         # Derived from the final gtf above. Null in quant-only-cluster-guided mode,
         # which reuses a supplied annotation instead of discovering a catalog to
         # collapse. The two beds are always written when the collapse runs; the
         # collapsed gtf and merge report are absent if a splice pattern arrived under
         # two gene_ids, which reports itself in the conflicts file rather than
         # failing the task and publishing nothing.
         File? LRAA_final_splice_pattern_collapsed_gtf = sc_splice_pattern_collapse.collapsedGtf
         File? LRAA_final_splice_pattern_collapsed_merge_report = sc_splice_pattern_collapse.mergeReport
         File? LRAA_final_splice_pattern_collapsed_gene_conflicts = sc_splice_pattern_collapse.geneConflictsReport
         File? LRAA_final_TSS_bed = sc_splice_pattern_collapse.tssBed
         File? LRAA_final_PolyA_bed = sc_splice_pattern_collapse.polyaBed
         # partitioned cluster BAMs (always produced)
         File LRAA_partitioned_cluster_bams_tar = tar_partitioned_cluster_bams.tar_gz
         # Per-cluster quant.expr files, one per cluster, each the pooled final
         # quantification summed over that cluster's cells (no per-cluster EM).
         File LRAA_final_cluster_exprs_tar = tar_final_cluster_expr_files.tar_gz
         # The pooled final quantification: its quant.expr over all cells, and its read
         # tracking, from which the single-cell matrices and the per-cluster files are
         # built. The tracking includes reads from barcodes in no cluster.
         File LRAA_final_quant_expr = final_quant_expr
         File LRAA_final_tracking = final_quant_tracking
         # Assigned reads per cluster, and of barcodes in no cluster, in the pooled
         # quantification.
         File LRAA_final_cluster_read_partition = split_pooled_quant_by_cluster.read_partition_summary

         # single-cell sparse matrices and intermediate counts
         File sc_gene_transcript_splicehash_mapping = sc_mapping_tsv
         File? sc_gene_counts = build_sc_sparse_matrices.gene_counts
         File? sc_isoform_counts = build_sc_sparse_matrices.isoform_counts
         File? sc_splice_pattern_counts = build_sc_sparse_matrices.splice_pattern_counts
         # gene/isoform/splice-pattern sparse directories as tar.gz
         File sc_gene_sparse_tar_gz = sc_gene_tgz
         File sc_isoform_sparse_tar_gz = sc_isoform_tgz
         File sc_splice_pattern_sparse_tar_gz = sc_splice_tgz

         # preliminary intermediate outputs (only in discovery mode)
         File? LRAA_prelim_cluster_gtfs = tar_cluster_prelim_gtf_files.tar_gz
         File? LRAA_prelim_cluster_read_trackings = tar_cluster_prelim_tracking_files.tar_gz
         File? LRAA_prelim_cluster_pseudobulk_exprs = tar_cluster_prelim_pseudobulk_expr_files.tar_gz

         # cluster pseudobulk matrices
         File cluster_gene_counts_matrix = build_cluster_pseudobulk.cluster_gene_counts_matrix
         File cluster_gene_TPM_matrix = build_cluster_pseudobulk.cluster_gene_TPM_matrix
         File cluster_isoform_counts_matrix = build_cluster_pseudobulk.cluster_isoform_counts_matrix
         File cluster_isoform_TPM_matrix = build_cluster_pseudobulk.cluster_isoform_TPM_matrix
         File cluster_isoform_counts_forDiffIsoUsage = build_cluster_pseudobulk.cluster_isoform_counts_forDiffIsoUsage

         # The ONE geometry every phase of this run cut on. Absent only when both
         # scattering knobs are off, the mode with nothing to share.
         File? shared_chunk_plan = shared_chunk_plan_resolved

         # Read-assignment accounting, one table per PHASE and never pooled across
         # them; rows of different level must not be summed, which the file states on
         # its first line. The prelim table (per-cluster discovery) carries level=chunk
         # rows, one level=cluster row per cluster and one level=all_clusters row, and
         # is absent in quant_only mode, which runs no phase-1 discovery. The final
         # table is the one pooled quantification: level=chunk rows and one row
         # labelled "pseudobulk".
         File? prelim_collated_read_assignment_summary = collate_prelim_read_assignment_summaries.collatedSummary
         File? collated_read_assignment_summary = collate_final_read_assignment_summary.collatedSummary
    
     }
     
}


task LRAA_tar_outputs {
    input {
        Int preemptible_tries = 3
        String tar_directory_name
        Array[File] input_files
        String docker
    }

    Int memoryGB = 8
    Int disksize = 20 + ceil(5 * size(input_files, "GiB"))
    # Machine shape (v0.47.1, generated text, identical in every N2D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). N2D rather than C3D because the task uses at most
    # 2 cores, and C3D has no 2-core shape, or stays on an HDD disk, and Compute Engine docs
    # appear not to offer pd-standard on C3D (unprobed).
    Int shape_cpu = 1
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 2 then 2
        else if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 32 then 32
        else if shape_cpu <= 48 then 48
        else if shape_cpu <= 64 then 64
        else if shape_cpu <= 80 then 80
        else 96
    Int shape_mem_tier = if shape_memGB <= 16 then 2
        else if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 256 then 32
        else if shape_memGB <= 384 then 48
        else if shape_memGB <= 512 then 64
        else if shape_memGB <= 640 then 80
        else 96
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    String machine_type = if shape_memGB <= shape_cores then "n2d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "n2d-standard-~{shape_cores}"
        else "n2d-highmem-~{shape_cores}"

    
    command <<<

        set -ex

        mkdir ~{tar_directory_name}

        for file in "~{sep=' ' input_files}"; do
           cp $file ~{tar_directory_name}/
        done
        
        tar -zcvf ~{tar_directory_name}.tar.gz ~{tar_directory_name}/

    >>>


    output {
        File tar_gz = "~{tar_directory_name}.tar.gz"
    }

    runtime {

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} SSD"
    }
}


task split_pooled_quant_by_cluster {
    # Per-cluster quant.expr files from the ONE pooled final quantification: every read
    # row of the tracking is credited to the cluster of its cell barcode and summed.
    # No EM; util/sc/tracking_to_cluster_pseudobulk_exprs.py documents each column.
    # Reads of barcodes in no cluster are counted in no file, and reported in
    # read_partition_summary.
    input {
        Int preemptible_tries = 3
        String sample_id
        File quant_tracking
        File quant_expr
        File cell_clusters_info
        String docker
        # One pass over the tracking, holding per (cluster, isoform) totals and the
        # barcodes seen; memory scales with isoforms x clusters plus barcodes, not with
        # the tracking's rows. MEASURED on the PBMC v0.50.0 ref-guided pooled tracking
        # (2.8 GB gz, 131 M rows, 14 clusters): 10.9 min wall, 1.3 GB peak RSS.
        Int memoryGB = 8
    }

    Int disksize = 20 + ceil(3 * size(quant_tracking, "GB"))
    # Machine shape (v0.47.1, generated text, identical in every N2D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). N2D rather than C3D because the task uses at most
    # 2 cores, and C3D has no 2-core shape, or stays on an HDD disk, and Compute Engine docs
    # appear not to offer pd-standard on C3D (unprobed).
    Int shape_cpu = 2
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 2 then 2
        else if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 32 then 32
        else if shape_cpu <= 48 then 48
        else if shape_cpu <= 64 then 64
        else if shape_cpu <= 80 then 80
        else 96
    Int shape_mem_tier = if shape_memGB <= 16 then 2
        else if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 256 then 32
        else if shape_memGB <= 384 then 48
        else if shape_memGB <= 512 then 64
        else if shape_memGB <= 640 then 80
        else 96
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    String machine_type = if shape_memGB <= shape_cores then "n2d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "n2d-standard-~{shape_cores}"
        else "n2d-highmem-~{shape_cores}"


    command <<<
        set -euo pipefail

        # The helper is INLINED, like merge_bams' header check in LRAA_quant_by_cluster.wdl:
        # this task runs in a separately pinned image while the .wdl comes from the
        # checkout, so calling a newly added script by name breaks whenever the two are
        # out of step. The text between the markers is util/sc/tracking_to_cluster_pseudobulk_exprs.py
        # verbatim (its unit tests are the tested source); pylib/test_inlined_wdl_scripts.py
        # fails when they differ. Edit the util script, then paste it here unchanged.
        cat > tracking_to_cluster_pseudobulk_exprs.py <<'LRAA_INLINED_SCRIPT_EOF'
#!/usr/bin/env python3
"""Per-cluster pseudobulk quant.expr files from ONE pooled single-cell quantification.

The cluster-guided single-cell workflow quantifies every cell against the final catalog
in a single LRAA quant-only pass (the pooled, all-cells EM). This script splits that one
assignment by cell cluster: each read row of the pooled tracking file is credited to the
cluster of its cell barcode, and the rows are summed per (cluster, isoform). No EM is run
here, so every cluster's ambiguous reads keep the fractions the pooled EM gave them.

Columns are computed exactly as LRAA's streaming quant.expr writer computes them
(pylib/StreamingQuant.py write_expr), restricted to the cluster's reads:

    uniq_reads                 read rows with is_unique == 1
    all_reads                  sum of frac_assigned
    isoform_fraction           all_reads / the gene's all_reads in the cluster
    unique_gene_read_fraction  uniq_reads / the gene's all_reads in the cluster
    TPM                        all_reads / the cluster's all_reads, x 1e6
    RPM_total_reads            see below
    uniq_FSM_reads             read rows with is_unique == 1 and is_FSM == 1
    has_FSM_read               1 if any read row has is_FSM == 1

The remaining columns (exons, introns, splice_hash_code and any others) are copied from
the pooled quant.expr, which also fixes the isoform list: every cluster file lists every
isoform of the pooled file, zeros included, sorted by transcript_id as LRAA writes them.

RPM_total_reads: LRAA divides by the mapped reads in the bam it quantified. A cluster's
own mapped-read total is not known here, so the denominator is the cluster's assigned
reads (its all_reads total), which makes this column equal TPM. A comment line in each
file says so.

Reads from barcodes that are in no cluster are counted in the pooled quantification (and
in the single-cell matrices built from it) but in no cluster file. Their totals are
reported in the --partition_summary table.

Before writing anything the script checks that the tracking reproduces the pooled
quant.expr: per isoform, the per-cluster sums plus the unclustered remainder must equal
the pooled uniq_reads, uniq_FSM_reads and has_FSM_read exactly; a mismatch means the two
files are not from the same run and stops the script. all_reads is compared to the pooled
file's one-decimal rounding and only warned about: LRAA's default streaming writer reports
the tracking's sum, while the EM report path (stream_reads off) reports EM's own counts.

Output files are named <output_prefix>.<cluster>.LRAA.quant-only.quant.expr, the names
the former per-cluster quant runs wrote, so readers that take the cluster from the file
name (util/sc/site_usage/build_site_event_read_tracks.py) keep working.
"""

import argparse
import collections
import gzip
import logging
import os
import shutil
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

UNCLUSTERED = "__unclustered__"
# all_reads is written with one decimal, so a correct sum can sit up to 0.05 away
ALL_READS_TOLERANCE = 0.051


def open_text(path):
    """Text handle on a plain or gzipped file; gzip goes through pigz when available."""
    if not path.endswith(".gz"):
        return open(path, "rt")
    if shutil.which("pigz"):
        proc = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE, text=True)
        return _ProcReader(proc, path)
    return gzip.open(path, "rt")


class _ProcReader:
    """File-like wrapper that checks the decompressor's exit status on close."""

    def __init__(self, proc, path):
        self.proc = proc
        self.path = path

    def __iter__(self):
        return iter(self.proc.stdout)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.proc.stdout.close()
        rc = self.proc.wait()
        if rc != 0 and exc[0] is None:
            raise RuntimeError(f"pigz -dc {self.path} exited {rc}")
        return False


def read_cell_clusters(path):
    """cell_barcode -> cluster, parsed as util/sc/partition_bam_by_cell_cluster.py does:
    two tab-separated fields per line; the exact header row 'cell_barcode<TAB>cluster'
    on line 1 is skipped; other lines without two fields are skipped."""
    cell_to_cluster = dict()
    with open(path, "rt") as fh:
        for lineno, line in enumerate(fh, start=1):
            vals = line.rstrip().split("\t")
            if len(vals) != 2:
                continue
            barcode, cluster = vals
            if lineno == 1 and (barcode, cluster) == ("cell_barcode", "cluster"):
                continue
            cell_to_cluster[barcode] = cluster
    if not cell_to_cluster:
        raise SystemExit(f"no cell -> cluster assignments read from {path}")
    return cell_to_cluster


def read_pooled_expr(path):
    """(comment lines, header fields, rows as lists) of the pooled quant.expr."""
    comments, header, rows = [], None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
                continue
            vals = line.rstrip("\n").split("\t")
            if header is None:
                header = vals
                continue
            rows.append(vals)
    if header is None:
        raise SystemExit(f"no header line in {path}")
    for col in ("gene_id", "transcript_id", "uniq_reads", "all_reads", "isoform_fraction",
                "unique_gene_read_fraction", "TPM", "RPM_total_reads", "uniq_FSM_reads",
                "has_FSM_read"):
        if col not in header:
            raise SystemExit(f"{path}: column {col} missing")
    return comments, header, rows


class Totals:
    __slots__ = ("frac", "uniq", "uniq_fsm", "any_fsm")

    def __init__(self):
        self.frac = 0.0
        self.uniq = 0
        self.uniq_fsm = 0
        self.any_fsm = False


def tally_tracking(path, cell_to_cluster, known_transcripts):
    """(cluster, transcript_id) -> Totals over the tracking rows, plus per-group counts of
    read rows and distinct barcodes. Barcodes in no cluster go to UNCLUSTERED."""
    totals = collections.defaultdict(Totals)
    barcodes_seen = collections.defaultdict(set)
    col = None
    n_rows = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            vals = line.rstrip("\n").split("\t")
            if vals[0] == "gene_id":
                if col is None:
                    col = {name: i for i, name in enumerate(vals)}
                    for name in ("transcript_id", "read_name", "frac_assigned", "is_unique", "is_FSM"):
                        if name not in col:
                            raise SystemExit(f"{path}: tracking column {name} missing")
                    i_tx, i_read = col["transcript_id"], col["read_name"]
                    i_frac, i_uniq, i_fsm = col["frac_assigned"], col["is_unique"], col["is_FSM"]
                continue
            if col is None:
                raise SystemExit(f"{path}: data row before the column header")
            n_rows += 1
            transcript_id = vals[i_tx]
            if transcript_id not in known_transcripts:
                raise SystemExit(
                    f"{path}: transcript {transcript_id} is not in the pooled quant.expr; "
                    "the tracking and quant.expr must come from the same quantification")
            barcode = vals[i_read].split("^", 1)[0]
            cluster = cell_to_cluster.get(barcode, UNCLUSTERED)
            barcodes_seen[cluster].add(barcode)
            t = totals[(cluster, transcript_id)]
            t.frac += float(vals[i_frac])
            if vals[i_fsm] == "1":
                t.any_fsm = True
            if vals[i_uniq] == "1":
                t.uniq += 1
                if vals[i_fsm] == "1":
                    t.uniq_fsm += 1
    if col is None:
        raise SystemExit(f"{path}: no tracking column header found")
    logger.info("read %d tracking rows", n_rows)
    return totals, barcodes_seen


def check_against_pooled(header, rows, totals):
    """Per isoform, clusters + unclustered must reproduce the pooled quant.expr."""
    ix = {name: i for i, name in enumerate(header)}
    by_tx = collections.defaultdict(Totals)
    for (cluster, tx), t in totals.items():
        s = by_tx[tx]
        s.frac += t.frac
        s.uniq += t.uniq
        s.uniq_fsm += t.uniq_fsm
        s.any_fsm = s.any_fsm or t.any_fsm
    bad, off = [], []
    for r in rows:
        tx = r[ix["transcript_id"]]
        s = by_tx.get(tx, Totals())
        if (int(r[ix["uniq_reads"]]) != s.uniq
                or int(r[ix["uniq_FSM_reads"]]) != s.uniq_fsm
                or (r[ix["has_FSM_read"]] == "1") != s.any_fsm):
            bad.append((tx, r[ix["uniq_reads"]], s.uniq, r[ix["uniq_FSM_reads"]], s.uniq_fsm,
                        r[ix["has_FSM_read"]], int(s.any_fsm)))
        pooled_all = float(r[ix["all_reads"]])
        if abs(pooled_all - s.frac) > ALL_READS_TOLERANCE + 1e-4 * pooled_all:
            off.append((tx, r[ix["all_reads"]], f"{s.frac:.3f}"))
    if off:
        # The streaming writer (LRAA's default, stream_reads) reports all_reads as the sum
        # of the tracking's frac_assigned, so this should not happen there. The EM report
        # path (stream_reads off) reports EM's own counts, which can differ slightly.
        for o in off[:10]:
            logger.warning("all_reads differs, pooled quant.expr vs tracking sum: %s %s / %s", *o)
        logger.warning("%d of %d isoforms: pooled all_reads differs from the tracking's sum of "
                       "frac_assigned; cluster values are tracking sums", len(off), len(rows))
    if bad:
        for b in bad[:10]:
            logger.error("pooled vs tracking mismatch: %s uniq %s/%s uniqFSM %s/%s hasFSM %s/%s", *b)
        raise SystemExit(
            f"{len(bad)} isoforms: the tracking does not reproduce the pooled quant.expr's "
            "uniq_reads, uniq_FSM_reads or has_FSM_read; the two files must come from the "
            "same quantification")


def write_cluster_expr(path, comments, header, rows, cluster, totals):
    ix = {name: i for i, name in enumerate(header)}
    gene_of = {r[ix["transcript_id"]]: r[ix["gene_id"]] for r in rows}
    gene_frac = collections.defaultdict(float)
    cluster_total = 0.0
    for r in rows:
        t = totals.get((cluster, r[ix["transcript_id"]]))
        if t is not None:
            gene_frac[gene_of[r[ix["transcript_id"]]]] += t.frac
            cluster_total += t.frac
    zero = Totals()
    with open(path, "wt") as ofh:
        for c in comments:
            print(c, file=ofh)
        print(f"# cluster {cluster}: derived from the pooled (all cells) quantification by "
              "summing the cluster's cells; RPM_total_reads is per million of the cluster's "
              "assigned reads (equal to TPM)", file=ofh)
        print("\t".join(header), file=ofh)
        for r in sorted(rows, key=lambda r: r[ix["transcript_id"]]):
            tx = r[ix["transcript_id"]]
            t = totals.get((cluster, tx), zero)
            gene_total = gene_frac.get(gene_of[tx], 0.0)
            out = list(r)
            out[ix["uniq_reads"]] = f"{t.uniq}"
            out[ix["all_reads"]] = f"{t.frac:.1f}"
            out[ix["isoform_fraction"]] = f"{(t.frac / gene_total if gene_total > 0 else 0.0):.3f}"
            out[ix["unique_gene_read_fraction"]] = f"{(t.uniq / gene_total if gene_total > 0 else 0.0):.3f}"
            tpm = t.frac / cluster_total * 1e6 if cluster_total > 0 else 0.0
            out[ix["TPM"]] = f"{tpm:.3f}"
            out[ix["RPM_total_reads"]] = f"{tpm:.3f}"
            out[ix["uniq_FSM_reads"]] = f"{t.uniq_fsm}"
            out[ix["has_FSM_read"]] = "1" if t.any_fsm else "0"
            print("\t".join(out), file=ofh)
    return cluster_total


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tracking", required=True,
                        help="pooled LRAA quant.tracking (.gz ok)")
    parser.add_argument("--quant_expr", required=True,
                        help="pooled LRAA quant.expr from the same quantification")
    parser.add_argument("--cell_clusters", required=True,
                        help="cell_barcode <tab> cluster (as given to partition_bam_by_cell_cluster.py)")
    parser.add_argument("--output_prefix", required=True,
                        help="files are written as <output_prefix>.<cluster>.LRAA.quant-only.quant.expr")
    parser.add_argument("--output_dir", default=".", help="directory for the per-cluster files")
    parser.add_argument("--partition_summary", default=None,
                        help="optional TSV: per cluster (and unclustered) barcodes with reads, "
                             "read rows' assigned reads, unique reads")
    args = parser.parse_args()

    cell_to_cluster = read_cell_clusters(args.cell_clusters)
    clusters = sorted(set(cell_to_cluster.values()))
    if UNCLUSTERED in clusters:
        raise SystemExit(f"a cluster is named {UNCLUSTERED}, which this script reserves")
    logger.info("%d barcodes in %d clusters", len(cell_to_cluster), len(clusters))

    comments, header, rows = read_pooled_expr(args.quant_expr)
    known = {r[header.index("transcript_id")] for r in rows}
    logger.info("%d isoforms in the pooled quant.expr", len(rows))

    totals, barcodes_seen = tally_tracking(args.tracking, cell_to_cluster, known)
    check_against_pooled(header, rows, totals)

    os.makedirs(args.output_dir, exist_ok=True)
    assigned = dict()
    for cluster in clusters:
        path = os.path.join(args.output_dir, f"{args.output_prefix}.{cluster}.LRAA.quant-only.quant.expr")
        assigned[cluster] = write_cluster_expr(path, comments, header, rows, cluster, totals)
        logger.info("cluster %s: %.1f assigned reads -> %s", cluster, assigned[cluster], path)

    group_frac = collections.defaultdict(float)
    group_uniq = collections.defaultdict(int)
    for (c, _), t in totals.items():
        group_frac[c] += t.frac
        group_uniq[c] += t.uniq
    unclustered = group_frac[UNCLUSTERED]
    pooled = sum(group_frac.values())
    logger.info("unclustered barcodes: %d, %.1f assigned reads (%.2f%% of %.1f)",
                len(barcodes_seen.get(UNCLUSTERED, ())), unclustered,
                100.0 * unclustered / pooled if pooled else 0.0, pooled)

    if args.partition_summary:
        with open(args.partition_summary, "wt") as ofh:
            print("# Reads of the pooled quantification by cell cluster. assigned_reads: sum of "
                  "frac_assigned; unique_reads: read rows assigned to one isoform only. The "
                  "__unclustered__ row holds barcodes in no cluster: counted in the single-cell "
                  "matrices, in no cluster's pseudobulk.", file=ofh)
            print("\t".join(["cluster", "barcodes_with_reads", "assigned_reads", "unique_reads",
                             "fraction_of_assigned_reads"]), file=ofh)
            for cluster in clusters + [UNCLUSTERED]:
                a, u = group_frac[cluster], group_uniq[cluster]
                print("\t".join([cluster, str(len(barcodes_seen.get(cluster, ()))), f"{a:.1f}", str(u),
                                 f"{(a / pooled if pooled else 0.0):.6f}"]), file=ofh)

    return 0


if __name__ == "__main__":
    sys.exit(main())
LRAA_INLINED_SCRIPT_EOF

        python3 tracking_to_cluster_pseudobulk_exprs.py \
            --tracking ~{quant_tracking} \
            --quant_expr ~{quant_expr} \
            --cell_clusters ~{cell_clusters_info} \
            --output_prefix ~{sample_id} \
            --output_dir cluster_exprs \
            --partition_summary ~{sample_id}.LRAA.final.cluster_read_partition.tsv
    >>>

    output {
        Array[File] cluster_quant_exprs = glob("cluster_exprs/*.quant.expr")
        File read_partition_summary = "~{sample_id}.LRAA.final.cluster_read_partition.tsv"
    }

    runtime {

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        # 1 core for the parse, 1 for pigz decompressing the tracking
        cpu: 2
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} SSD"
    }
}




task lraa_merge_gtf_task {
    input {
        Int preemptible_tries = 1
        String sample_id
        Array[File] LRAA_cell_cluster_gtfs
        File? pooled_gtf
        File referenceGenome
        # TSS/PolyA annotations on the input gtfs are respected by default; set true to
        # ignore them (the former non-HiFi merge). No longer gated on --HiFi, since LRAA
        # models TSS/PolyA in both modes.
        Boolean ignore_TSS_POLYA = false
        # An oversimplified contig was never assembled, so there is nothing across
        # the per-cluster gtfs to reconcile. Without this, the merge rebuilds a
        # splice graph over its aggregate models and reconstructs them, which moved
        # a measured chrM from 1-16569 to 1-16570 -- one base past the contig end --
        # and replaced the OVSIMP ids with comp-N:iso-N. The merge copies these
        # contigs' records verbatim from the inputs instead, and refuses if the
        # inputs disagree about them.
        String oversimplify = ""
        # LRAA_Globals.config overrides (JSON object), handed to merge_LRAA_GTFs.py as
        # --config_update. Unset: the merge runs on LRAA's defaults.
        String? config_update_json
        String docker
        Int memoryGB
        # Contig/strand units are merged in parallel worker processes. The script
        # runs --CPU auto: the cores the backend actually grants, capped by the
        # memory granted (6 GB budgeted per worker).
        Int cpu = 8
    }
    # Machine shape (v0.47.1, generated text, identical in every C3D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). The command's thread count stays on the
    # requested cpu even when the shape has more cores. Batch rejects 8 GiB on c3d-highcpu-4
    # (7.5 GiB passes), hence 7 for that tier: such a request steps up to c3d-standard-4.
    Int shape_cpu = cpu
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 30 then 30
        else if shape_cpu <= 60 then 60
        else if shape_cpu <= 90 then 90
        else 180
    Int shape_mem_tier = if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 240 then 30
        else if shape_memGB <= 480 then 60
        else if shape_memGB <= 720 then 90
        else 180
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    Int shape_highcpu_memGB = if shape_cores == 4 then 7
        else if shape_cores == 30 then 59
        else if shape_cores == 60 then 118
        else if shape_cores == 90 then 177
        else if shape_cores == 180 then 354
        else shape_cores * 2
    String machine_type = if shape_memGB <= shape_highcpu_memGB then "c3d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "c3d-standard-~{shape_cores}"
        else "c3d-highmem-~{shape_cores}"


    
    command <<<
        set -ex

        lraa_version="$(LRAA --version | awk '{print $NF}')"
        version_comment="# LRAA version ${lraa_version}"


        # Written to a file with a quoted heredoc, so the shell never sees the JSON.
        config_update_file="$(pwd)/merge_config_update.json"
        cat > "$config_update_file" <<'LRAA_CONFIG_UPDATE_EOF'
~{default="" config_update_json}
LRAA_CONFIG_UPDATE_EOF

      (
      
        merge_LRAA_GTFs.py --genome ~{referenceGenome} \
                           ~{if ignore_TSS_POLYA then "--ignore_TSS_POLYA" else ""} \
                           ~{if select_first([config_update_json, ""]) != "" then "--config_update \"$config_update_file\"" else ""} \
                           ~{if (oversimplify != "") then "--oversimplify '" + oversimplify + "'" else ""} \
                           --CPU auto \
                           --gtf ~{sep=' ' LRAA_cell_cluster_gtfs } \
                           ~{"--pooled_gtf " + pooled_gtf} \
                           --output_gtf ~{sample_id}.LRAA.sc_merged.gtf  > command_output.log 2>&1
      ) || {
             echo "Command failed with exit code $?" >&2
             echo "Last 100 lines of output:" >&2
             tail -n 100 command_output.log >&2
             exit 1
      }

      lraa_merge_header.py \
          --version_comment "$version_comment" \
          --inputs ~{sep=' ' LRAA_cell_cluster_gtfs} \
          --output merge_header.txt

      prepend_lraa_merge_header() {
          local path="$1"
          local tmp="${path}.with_lraa_header"
          if [[ -s "$path" ]] && head -n 1 "$path" | grep -qxF "$version_comment"; then
              return
          fi
          cat merge_header.txt "$path" > "$tmp"
          mv "$tmp" "$path"
      }

      prepend_lraa_merge_header ~{sample_id}.LRAA.sc_merged.gtf

      if [[ -f ~{sample_id}.LRAA.sc_merged.gtf.tracking.tsv ]]; then
          prepend_lraa_merge_header ~{sample_id}.LRAA.sc_merged.gtf.tracking.tsv
      fi

    >>>


    output {
        File mergedGTF = "~{sample_id}.LRAA.sc_merged.gtf"
        File? mergedTracking = "~{sample_id}.LRAA.sc_merged.gtf.tracking.tsv"
    }

    runtime {

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk 200 SSD"
    }

}


task build_cluster_pseudobulk_matrices {
    input {
        Int preemptible_tries = 3
        String sample_id
        Array[File] quant_expr_files
        String docker
        Int memoryGB = 8
    }

    Int disksize = 20 + ceil(10 * length(quant_expr_files))

    String output_prefix = "~{sample_id}.LRAA.cluster_pseudobulk"
    # Machine shape (v0.47.1, generated text, identical in every N2D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). N2D rather than C3D because the task uses at most
    # 2 cores, and C3D has no 2-core shape, or stays on an HDD disk, and Compute Engine docs
    # appear not to offer pd-standard on C3D (unprobed).
    Int shape_cpu = 1
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 2 then 2
        else if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 32 then 32
        else if shape_cpu <= 48 then 48
        else if shape_cpu <= 64 then 64
        else if shape_cpu <= 80 then 80
        else 96
    Int shape_mem_tier = if shape_memGB <= 16 then 2
        else if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 256 then 32
        else if shape_memGB <= 384 then 48
        else if shape_memGB <= 512 then 64
        else if shape_memGB <= 640 then 80
        else 96
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    String machine_type = if shape_memGB <= shape_cores then "n2d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "n2d-standard-~{shape_cores}"
        else "n2d-highmem-~{shape_cores}"


    command <<<
        set -ex

        build_LRAA_expr_matrices.py \
            --output_prefix ~{output_prefix} \
            --quant_files ~{sep=' ' quant_expr_files}
    >>>

    output {
        File cluster_gene_counts_matrix = "~{output_prefix}.gene.counts.matrix"
        File cluster_gene_TPM_matrix = "~{output_prefix}.gene.TPM.matrix"
        File cluster_isoform_counts_matrix = "~{output_prefix}.isoform.counts.matrix"
        File cluster_isoform_TPM_matrix = "~{output_prefix}.isoform.TPM.matrix"
        File cluster_isoform_counts_forDiffIsoUsage = "~{output_prefix}.isoform.counts.matrix.forDiffIsoUsage"
    }

    runtime {

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} SSD"
    }
}


task sc_build_sparse_matrices {
    input {
        Int preemptible_tries = 3
        String sample_id
        File tracking_file
        String docker
        # PROVISIONAL -- see LRAA-build_sparse_matrices_from_tracking.wdl.
        Int memoryGB = 64
        String csv_engine = "direct"
        Int gzip_level = 1
        Boolean emit_dense_counts = false
    }

    Int disksize = 50 + ceil(2 * size(tracking_file, "GB"))

    String output_prefix = "~{sample_id}.LRAA.sc"
    # Machine shape (v0.47.1, generated text, identical in every C3D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). The command's thread count stays on the
    # requested cpu even when the shape has more cores. Batch rejects 8 GiB on c3d-highcpu-4
    # (7.5 GiB passes), hence 7 for that tier: such a request steps up to c3d-standard-4.
    Int shape_cpu = 3
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 30 then 30
        else if shape_cpu <= 60 then 60
        else if shape_cpu <= 90 then 90
        else 180
    Int shape_mem_tier = if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 240 then 30
        else if shape_memGB <= 480 then 60
        else if shape_memGB <= 720 then 90
        else 180
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    Int shape_highcpu_memGB = if shape_cores == 4 then 7
        else if shape_cores == 30 then 59
        else if shape_cores == 60 then 118
        else if shape_cores == 90 then 177
        else if shape_cores == 180 then 354
        else shape_cores * 2
    String machine_type = if shape_memGB <= shape_highcpu_memGB then "c3d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "c3d-standard-~{shape_cores}"
        else "c3d-highmem-~{shape_cores}"


    command <<<
        set -ex

        singlecell_tracking_to_sparse_matrix.py \
            --tracking ~{tracking_file} \
            --output_prefix ~{output_prefix} \
            --csv_engine ~{csv_engine} \
            --gzip_level ~{gzip_level} \
            ~{true="--emit_dense_counts" false="" emit_dense_counts}

        # Tar the generated sparse matrix directories for compact output.
        # NOTE: archive filenames use "." rather than the "^" that
        # singlecell_tracking_to_sparse_matrix.py uses for the source directory
        # names: these tar.gz files are passed as File inputs to downstream tasks,
        # and Apptainer's --bind spec parser mis-splits paths containing a literal
        # "^" (Docker's bind mounts are unaffected, but the WDL must work under both).
        if command -v pigz >/dev/null 2>&1; then
            tar --use-compress-program="pigz -~{gzip_level}" -cvf "~{output_prefix}.gene-sparseM.tar.gz" "~{output_prefix}^gene-sparseM" || true
            tar --use-compress-program="pigz -~{gzip_level}" -cvf "~{output_prefix}.isoform-sparseM.tar.gz" "~{output_prefix}^isoform-sparseM" || true
            tar --use-compress-program="pigz -~{gzip_level}" -cvf "~{output_prefix}.splice_pattern-sparseM.tar.gz" "~{output_prefix}^splice_pattern-sparseM" || true
        else
            tar -zcvf "~{output_prefix}.gene-sparseM.tar.gz" "~{output_prefix}^gene-sparseM" || true
            tar -zcvf "~{output_prefix}.isoform-sparseM.tar.gz" "~{output_prefix}^isoform-sparseM" || true
            tar -zcvf "~{output_prefix}.splice_pattern-sparseM.tar.gz" "~{output_prefix}^splice_pattern-sparseM" || true
        fi
    >>>

    output {
        File mapping_file = "~{output_prefix}.gene_transcript_splicehashcode.tsv"
        File? gene_counts = "~{output_prefix}.gene_cell_counts.tsv"
        File? isoform_counts = "~{output_prefix}.isoform_cell_counts.tsv"
        File? splice_pattern_counts = "~{output_prefix}.splice_pattern_cell_counts.tsv"

        # tar.gz of sparse matrix directories
        File gene_sparse_dir_tgz = "~{output_prefix}.gene-sparseM.tar.gz"
        File isoform_sparse_dir_tgz = "~{output_prefix}.isoform-sparseM.tar.gz"
        File splice_pattern_sparse_dir_tgz = "~{output_prefix}.splice_pattern-sparseM.tar.gz"
    }

    runtime {

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        cpu: 3
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} SSD"
    }
}


task require_annot_gtf {
    input {
        Int preemptible_tries = 3
        File annot_gtf
        # REQUIRED, and the omission of this input was a real Terra failure:
        # GCP Batch rejects any task with no image ("No container image found in
        # either 'container' or 'docker' runtime attributes"), while miniwdl
        # substitutes a default image and runs the task happily. So this task
        # passed every local test -- including three full by_chunk single-cell
        # smoke runs -- and failed only on Terra, in quant_only mode, which is
        # the one configuration that reaches it. Declared without a default so a
        # future call site cannot omit it silently.
        String docker
    }
    # Machine shape (v0.47.1, generated text, identical in every N2D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). N2D rather than C3D because the task uses at most
    # 2 cores, and C3D has no 2-core shape, or stays on an HDD disk, and Compute Engine docs
    # appear not to offer pd-standard on C3D (unprobed).
    Int shape_cpu = 1
    Int shape_memGB = 1
    Int shape_cpu_tier = if shape_cpu <= 2 then 2
        else if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 32 then 32
        else if shape_cpu <= 48 then 48
        else if shape_cpu <= 64 then 64
        else if shape_cpu <= 80 then 80
        else 96
    Int shape_mem_tier = if shape_memGB <= 16 then 2
        else if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 256 then 32
        else if shape_memGB <= 384 then 48
        else if shape_memGB <= 512 then 64
        else if shape_memGB <= 640 then 80
        else 96
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    String machine_type = if shape_memGB <= shape_cores then "n2d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "n2d-standard-~{shape_cores}"
        else "n2d-highmem-~{shape_cores}"


    command <<<
        set -euo pipefail
        echo "Annotation GTF provided: ~{annot_gtf}"
    >>>

    output {
        File ok = annot_gtf
    }

    runtime {

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 SSD"
    }
}


