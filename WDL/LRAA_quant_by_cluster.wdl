version 1.0

import "subwdls/Normalize_bam.wdl" as NormBam
import "subwdls/partition_bam_by_cell_cluster.wdl" as PartitionBam
import "LRAA.wdl" as LRAA

# LRAA_quant_by_cluster.wdl
# Workflow for cluster-specific quantification with normalized splice graph
#
# Strategy:
# 1. Normalize each input BAM in parallel
# 2. Merge all normalized BAMs
# 3. Normalize the merged BAM (final normalization)
# 4. Quantify each original BAM in parallel using the normalized merged BAM for splice graph construction

workflow LRAA_quant_by_cluster {
    input {
        String sample_id
        File referenceGenome
        File annot_gtf
        
        # Either provide bam_files directly OR provide inputBAM + cell_clusters_info to partition
        Array[File]? bam_files
        File? inputBAM
        File? cell_clusters_info
        Array[File]? pre_normalized_cluster_bams
        Array[File]? pre_normalized_cluster_bais

        # INTERNAL PLUMBING, not a user knob: the ONE cut plan for the WHOLE RUN,
        # emitted by the calling workflow BEFORE its phase-1 discovery so that the
        # caller's initial pass, its per-cluster discovery and this final quant all
        # cut at identical positions. When supplied, the emit_shared_chunk_plan call
        # below does not run; see it for why one plan per RUN rather than one per
        # phase, and for what a standalone call to this workflow still gets.
        File? internal_chunk_plan
        
        Boolean HiFi = false
        # Chunk geometry, forwarded to the per-cluster quantification runs. Unset
        # means LRAA's own defaults (10 Mb spacing, 1 Mb window).
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        String scattering = "by_chromosome"
        String? oversimplify
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int normalize_max_cov_level = 1000

        # Turn off 3' end weighting in the EM. LRAA_wf has declared this since the
        # negated-flag threading landed, but no single-cell entry point forwarded it,
        # so it was unreachable from Terra on every cluster path -- and 3'-biased
        # chemistry is exactly where you would want to measure without it.
        Boolean no_weight_reads_by_3prime_agreement = false
        
        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
        
        # Chromosome splitting parameters for LRAA quantification
        String main_chromosomes = "" # Set to split by chromosomes, leave empty to run without splitting
        
        # Cores per LRAA task: the task's cpu request AND the --cpu_budget it divides
        # across work units. There is no second knob to multiply it by.
        Int cpu = 2
        Int memoryGB_normalize = 8
        Int memoryGB_merge = 16
        # Optional override for direct quant-only LRAA.wdl calls per cluster.
        Int? memoryGB_quant
        # Optional override for chromosome-sharded quant-only LRAA.wdl workers only.
        # This applies only when main_chromosomes is non-empty for the per-cluster quantification calls below.
        Int? memoryGB_quant_scattered
        # Optional by_chunk overrides: make-chunks, chunk leaves, and merge need
        # independently sized resources. Unset preserves LRAA.wdl's defaults.
        Int? chunkMakeChunksCpu
        Int? chunkMakeChunksMemoryGB
        Int? chunkCpu
        Int? chunkMemoryGB
        Int? chunkMergeCpu
        Int? chunkMergeMemoryGB
        # Used only for chromosome-sharded per-cluster quantification runs.
        # OPTIONAL: forwarded to LRAA.wdl's cpuScattered, which resolves
        # select_first([cpuScattered, shard_cpu_computed]). A concrete value here
        # disables that per-shard estimate for every per-cluster quant.
        Int? cpu_scattered
        # The one-off shared-chunk-plan task below. Cut selection runs one concurrent
        # unit per contig under --cpu_budget, so cpu IS that budget.
        Int cpu_chunk_plan = 16
        # ChunkedRun charges a make-chunks unit 300 MiB in its own concurrency guard
        # (PREP_UNIT_PEAK_MIB, pylib/ChunkedRun.py:922), so 16 concurrent contig
        # selections are estimated at 4.7 GiB and 8 GiB is that plus headroom. The
        # guard reads MemAvailable and REDUCES concurrency rather than letting the box
        # OOM, so undersizing this costs wall time, not the run.
        Int memoryGB_chunk_plan = 8
        # WORKING disk for the per-cluster LRAA_wf calls below: they hand it to
        # subwdls/LRAA_runner.wdl, which spends it as `disks: "local-disk N HDD"` --
        # not the boot disk, which that task pins at 30. 256 is the same default
        # LRAA.wdl, LRAA-singlecell.wdl and LRAA-cell_cluster_guided.wdl already
        # carry, so leaving it unset changes nothing about how this workflow ran
        # before; a caller that already holds the value can now raise it.
        Int diskSizeGB = 256
        
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    }

    # If bam_files not provided, partition the input BAM by cell clusters
    if (!defined(bam_files)) {
        call PartitionBam.partition_bam_by_cell_cluster {
            input:
                sample_id = sample_id,
                cell_clusters_info = select_first([cell_clusters_info]),
                inputBAM = select_first([inputBAM]),
                cell_barcode_tag = cell_barcode_tag,
                # Same restriction as the cluster-guided caller applies to its own
                # partition call. Needed here too: this branch runs when this
                # workflow is called standalone without bam_files, and its output
                # feeds the same normalize/merge/normalize chain below.
                main_chromosomes = main_chromosomes,
                docker = docker,
                memoryGB = memoryGB_normalize
        }
    }

    # Select the BAM files to use: either provided or partitioned
    Array[File] cluster_bams = select_first([bam_files, partition_bam_by_cell_cluster.partitioned_bams])
    Boolean use_pre_normalized_bams = defined(pre_normalized_cluster_bams)

    if (use_pre_normalized_bams) {
        call validate_pre_normalized_inputs {
            input:
                cluster_bams = cluster_bams,
                normalized_bams = select_first([pre_normalized_cluster_bams]),
                normalized_bais = pre_normalized_cluster_bais,
                docker = docker
        }
    }

    # Step 1: Normalize each input BAM in parallel
    if (!use_pre_normalized_bams) {
        scatter (i in range(length(cluster_bams))) {
            String cluster_label = basename(cluster_bams[i], ".bam")
            
            call NormBam.normalize_bam_by_strand as normalize_cluster_bam {
                input:
                    input_bam = cluster_bams[i],
                    normalize_max_cov_level = normalize_max_cov_level,
                    label = cluster_label,
                    docker = docker,
                    memoryGB = memoryGB_normalize
            }
        }
    }

    Array[File] normalized_cluster_bams_use = select_first([pre_normalized_cluster_bams, normalize_cluster_bam.normalized_bam])
    Array[File] normalized_cluster_bais_use = select_first([pre_normalized_cluster_bais, normalize_cluster_bam.normalized_bai])

    # Step 2: Merge all normalized BAMs
    call merge_bams {
        input:
            normalized_bams = normalized_cluster_bams_use,
            output_basename = sample_id + ".clusters_merged",
            docker = docker,
            memoryGB = memoryGB_merge
    }

    # Step 3: Normalize the merged BAM (final normalization for splice graph)
    call NormBam.normalize_bam_by_strand as normalize_merged_bam {
        input:
            input_bam = merge_bams.merged_bam,
            normalize_max_cov_level = normalize_max_cov_level,
            label = sample_id + ".clusters_merged_final",
            docker = docker,
            memoryGB = memoryGB_normalize
    }

    # Step 3b: ONE chunk plan, shared by every cluster -- the FALLBACK for a
    # standalone call to this workflow.
    #
    # Two properties make the clusters comparable. All are quantified against one
    # consolidated GTF, and all build their splice graph from the one normalized BAM
    # produced above -- but cut POSITIONS are otherwise chosen per run. Left alone,
    # each cluster selects cuts on its OWN BAM, gets its own chunk boundaries, and
    # slices that shared normalized BAM at those boundaries, so the splice graph each
    # cluster actually saw differs and their boundary overhang drops differ with it.
    # The second property is then not achieved at all.
    #
    # Runs ONLY when the caller did not already emit one. Emitting here serves this
    # final quant and nothing earlier, so per-cluster DISCOVERY still selected its own
    # cuts -- measured at 76 distinct geometries over 270 extractions on the
    # ref-guided smoke run. LRAA-cell_cluster_guided.wdl now emits the plan before its
    # phase-1 discovery and passes it in as internal_chunk_plan, so a cluster-guided
    # run has ONE geometry end to end and this call is skipped. It stays for the
    # workflow's own entry point: a standalone LRAA_quant_by_cluster must not fall back
    # to per-cluster geometry.
    #
    # Selected on the WHOLE pre-partition BAM, which is an unthinned SUPERSET of every
    # cluster BAM: a position no read spans in the superset is spanned in no subset, so
    # one plan is safe for all of them. Selecting on normalize_merged_bam's output
    # instead would be worse than per-cluster selection -- normalization thins reads,
    # so a position that looks free there can still be spanned by a raw cluster read,
    # which extraction then DROPS with nothing reporting the loss, breaking
    # chunked-vs-unchunked quant parity.
    #
    # Produced for by_chunk AND by_chromosome, which is why the gate below tests only
    # for "off". by_chromosome is not an unchunked path: subwdls/LRAA_runner.wdl passes
    # --chunk unconditionally, so each chromosome shard chunks its contig inside itself
    # and, unplanned, picks those cuts from that cluster's own reads. Since
    # by_chromosome is this workflow's DEFAULT scattering, leaving it out would mean the
    # default single-cell final quant still had per-cluster geometry. off is the one
    # mode with nothing to share -- one whole-genome invocation per cluster, chunked
    # internally, and LRAA.wdl's validate_scattering refuses a plan there.
    #
    # And a plan is only possible when this workflow was handed the pre-partition BAM.
    # With bam_files supplied directly and no inputBAM there is no single superset
    # FILE -- the union of the cluster BAMs is not a file that exists, and merging them
    # here to manufacture one would be inventing an input. Nothing is passed in that
    # case, so ChunkedRun's own refusal (--bam_for_sg without shared cut geometry)
    # fails the run with the remedy named, rather than quietly giving every cluster its
    # own boundaries. The remedy for a caller is to pass inputBAM alongside bam_files,
    # which is what LRAA-cell_cluster_guided.wdl now does.
    if (scattering != "off" && defined(inputBAM) && !defined(internal_chunk_plan)) {
        call emit_shared_chunk_plan {
            input:
                inputBAM = select_first([inputBAM]),
                referenceGenome = referenceGenome,
                # THIS workflow's annot_gtf, which is the SAME gtf the per-cluster
                # calls below quantify against. That equality is what makes emitting
                # HERE safe; it is also why emitting here cannot serve an earlier
                # phase, since a plan selected on the consolidated GTF does not exist
                # until discovery has finished. The run-wide plan the caller emits
                # instead is selected on the REFERENCE annotation, before any novel
                # model exists, and that is safe for a different reason: discovery
                # then runs INSIDE those chunks, so no novel isoform can be built
                # spanning a cut -- no chunk ever held the reads that would support
                # one. Measured on the smoke runs, 0 of 994 consolidated ref-guided
                # models and 0 of 865 de novo models straddle any of the 9 cut
                # positions, and those were discovered under separate per-phase
                # geometries, which is the harsher test.
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

    # The geometry every cluster below cuts on: the caller's run-wide plan when it
    # emitted one, else the fallback above. Absent when scattering is off, or when no
    # pre-partition BAM was available to select from.
    File? shared_chunk_plan_use = if defined(internal_chunk_plan) then internal_chunk_plan else emit_shared_chunk_plan.chunk_plan

    # Step 4: Quantify each original BAM in parallel using normalized merged BAM for splice graph
    scatter (i in range(length(cluster_bams))) {
        String cluster_sample_id = basename(cluster_bams[i], ".bam")
        
        call LRAA.LRAA_wf as LRAA_quant_cluster {
            input:
                sample_id = cluster_sample_id,
                referenceGenome = referenceGenome,
                inputBAM = cluster_bams[i],
                # THREE bams, three roles, and they must be three different files:
                #
                #   inputBAM                 = this cluster's FULL reads. Pass 2
                #                              assigns these, so the counts reported
                #                              are cluster-local.
                #   internal_bam_for_sg      = the shared merged normalized bam. The
                #                              splice graph ONLY, taken as given and
                #                              never reconstructed, so every cluster is
                #                              quantified against one structure set.
                #   internal_bam_for_priors  = THIS cluster's own normalized reads.
                #                              Pass-1 theta.
                #
                # The third one is what keeps the first two honest. Unset, pass 1 falls
                # back to the splice-graph bam under stream_reads, which is ONE file
                # shared by every cluster: theta comes from POOLED evidence, so each
                # cluster's ambiguous reads are apportioned by every other cluster's
                # expression while its totals still look cluster-local (observed: 32
                # clusters all reporting reads_total 94,908). The normalized per-cluster
                # bams already exist -- step 1 above produced them for the merge -- so
                # this costs no new normalization. Supplying the sg bam also means LRAA
                # does not normalize again: LRAA.wdl derives no_norm from it being
                # defined, which is why there is no no_norm argument here any more.
                #
                # num_total_reads is deliberately NOT passed either. LRAA_wf counts the
                # bam it is handed, which here is this cluster's, so the TPM denominator
                # is THIS CLUSTER's read total and per-cluster TPM stays
                # cluster-relative. Forwarding a whole-library count would silently
                # rescale every cluster's TPM by that cluster's share of the library --
                # numbers that still look like TPM. The omission IS the semantics.
                internal_bam_for_sg = normalize_merged_bam.normalized_bam,
                internal_bam_for_sg_index = normalize_merged_bam.normalized_bai,
                internal_bam_for_priors = normalized_cluster_bams_use[i],
                internal_bam_for_priors_index = normalized_cluster_bais_use[i],
                # The ONE chunk geometry every cluster cuts on. Without it each cluster
                # would slice the shared splice-graph BAM at its own boundaries.
                internal_chunk_plan = shared_chunk_plan_use,
                annot_gtf = annot_gtf,
                quant_only = true,
                # single-cell: there is nothing to export, since the pre-normalized
                # splice-graph bam above means no normalization runs inside LRAA, and this
                # workflow surfaces its normalized BAMs from Normalize_bam.wdl instead
                retain_normalized_splice_graph_bam = false,
                no_EM = false,
                no_weight_reads_by_3prime_agreement = no_weight_reads_by_3prime_agreement,
                HiFi = HiFi,
                oversimplify = oversimplify,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                main_chromosomes = main_chromosomes,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                scattering = scattering,
                approx_MB_per_cut = approx_MB_per_cut,
                approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                cpu = cpu,
                cpuScattered = cpu_scattered,
                memoryGB = memoryGB_quant,
                memoryGBPerWorkerScattered = memoryGB_quant_scattered,
                chunkMakeChunksCpu = chunkMakeChunksCpu,
                chunkMakeChunksMemoryGB = chunkMakeChunksMemoryGB,
                chunkCpu = chunkCpu,
                chunkMemoryGB = chunkMemoryGB,
                chunkMergeCpu = chunkMergeCpu,
                chunkMergeMemoryGB = chunkMergeMemoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker
        }

        # Cluster identity for this cluster's read-assignment summary, carried
        # EXPLICITLY rather than by position.
        #
        # LRAA_wf.mergedReadAssignmentSummary is File? and genuinely can be absent:
        # LRAA copies the table only when stage 6 actually merged one (LRAA:3821-3831),
        # and LRAA_runner.wdl declares it File? for the same reason. select_all()
        # COMPACTS the gathered Array[File?], so one cluster without a summary shifts
        # every later element and the i-th surviving file is no longer cluster i.
        # Gating this String on the SAME defined() predicate select_all applies means
        # both arrays are filtered by one condition, in one order, and stay aligned
        # pair for pair. Do not replace this with a filename parse: the basename is
        # <cluster>.<LRAA.quant-only|ref-guided|ref-free>.read_assignment.summary.tsv,
        # so recovering the cluster means stripping a suffix that varies with run mode.
        if (defined(LRAA_quant_cluster.mergedReadAssignmentSummary)) {
            String cluster_id_for_summary = cluster_sample_id
        }
    }

    # Filtered by ONE predicate, so element k of each names the same cluster.
    Array[File] cluster_read_assignment_summaries = select_all(LRAA_quant_cluster.mergedReadAssignmentSummary)
    Array[String] cluster_read_assignment_ids = select_all(cluster_id_for_summary)

    # One TSV for the whole per-cluster phase: every cluster's per-chunk rows labelled
    # with their chunk, each cluster's own total, and one aggregate across clusters.
    # Gated on there being anything to collate -- the util refuses an empty input list
    # rather than writing an empty table, and every cluster's summary can legitimately
    # be absent on a build of LRAA that predates the per-unit summaries.
    if (length(cluster_read_assignment_summaries) > 0) {
        call collate_read_assignment_summaries as collate_cluster_read_assignment_summaries {
            input:
                summaryFiles = cluster_read_assignment_summaries,
                clusterIds = cluster_read_assignment_ids,
                outputFilePrefix = sample_id + ".clusters",
                # Each cluster's pass-1 reads were thinned by ITS OWN coverage
                # normalization, so the all_clusters total is a sum over independently
                # thinned populations and is not comparable to a pooled figure.
                population = "cluster_thinned",
                docker = docker
        }
    }

    output {
        # Per-cluster quantification outputs
        Array[File] quant_exprs = LRAA_quant_cluster.mergedQuantExpr
        Array[File] quant_trackings = LRAA_quant_cluster.mergedQuantTracking
        Array[Array[File]] read_assignment_shard_summaries_by_cluster = LRAA_quant_cluster.shardReadAssignmentSummaries
        Array[File] read_assignment_merged_summaries_by_cluster = cluster_read_assignment_summaries
        # The cluster ids of the array above, same order, same length. Published
        # because select_all() drops absent summaries and nothing else in this output
        # set records WHICH clusters survived that filter.
        Array[String] read_assignment_merged_summary_cluster_ids = cluster_read_assignment_ids
        # Every cluster's read-assignment accounting in ONE table: level=chunk rows
        # carrying cluster_id and chunk_id, one level=cluster row per cluster, one
        # level=all_clusters row. Rows of different level MUST NOT be summed; the file
        # says so on its first line.
        File? collated_read_assignment_summary = collate_cluster_read_assignment_summaries.collatedSummary
        
        # Intermediate outputs (for debugging/reuse)
        Array[File]? partitioned_bams = partition_bam_by_cell_cluster.partitioned_bams
        Array[File] normalized_cluster_bams = normalized_cluster_bams_use
        Array[File] normalized_cluster_bais = normalized_cluster_bais_use
        File merged_bam = merge_bams.merged_bam
        File merged_bai = merge_bams.merged_bai
        File normalized_merged_bam = normalize_merged_bam.normalized_bam
        File normalized_merged_bai = normalize_merged_bam.normalized_bai
        # The geometry every cluster above was cut on, whether this workflow emitted it
        # or the caller did. Absent when scattering is off, or when no pre-partition BAM
        # was available to select it from.
        File? shared_chunk_plan = shared_chunk_plan_use
    }
}

task validate_pre_normalized_inputs {
    input {
        Array[File] cluster_bams
        Array[File] normalized_bams
        Array[File]? normalized_bais
        String docker
    }

    command <<<
        set -euo pipefail

        cluster_count=~{length(cluster_bams)}
        normalized_count=~{length(normalized_bams)}
        if [[ "$cluster_count" -ne "$normalized_count" ]]; then
            echo "ERROR: cluster_bams count ($cluster_count) must equal pre_normalized_cluster_bams count ($normalized_count)." >&2
            exit 1
        fi

        ~{if defined(normalized_bais) then "bai_count=" + length(select_first([normalized_bais])) + "\nif [[ \"$normalized_count\" -ne \"$bai_count\" ]]; then\n    echo \"ERROR: pre_normalized_cluster_bams count ($normalized_count) must equal pre_normalized_cluster_bais count ($bai_count).\" >&2\n    exit 1\nfi" else "echo \"No pre_normalized_cluster_bais provided; proceeding without explicit BAI list check.\""}
    >>>

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 20 SSD"
    }
}


task merge_bams {
    input {
        Array[File] normalized_bams
        String output_basename
        String docker
        Int memoryGB = 16
        Int cpu = 4
    }

    Int disksize = 100 + ceil(3 * size(normalized_bams, "GB"))

    command <<<
        set -euo pipefail

        # --no-PG: samtools appends one @PG record PER EXISTING CHAIN TIP, so
        # merging N inputs concatenates all N @PG chains and then adds a record
        # for every resulting tip. This is the largest of the three merge
        # generations in the cluster-guided path -- measured on 13 cluster BAMs
        # each carrying 81,600 records, it produced 1,212,224 (1,060,800
        # inherited + 151,424 added). The chain originates upstream of LRAA:
        # XP132160.ucsc.bam arrives with 34,976 @PG records across 5,824
        # parallel chains, one per minimap2 alignment shard, never collapsed.
        #
        # Left unchecked the final merged BAM's header reached 2,727,296 records
        # = 1,188,154,439 bytes of UNCOMPRESSED SAM header text, in a BAM of
        # 3.68 GiB on disk (the header's own on-disk footprint is smaller, being
        # BGZF-compressed -- not a like-for-like ratio). Because resolving a region
        # name to a tid forces a full header parse, each per-chromosome region
        # query in Partition_data_by_chromosome cost ~5 min of pure header parsing
        # (~2h07m per cluster, invariant of alignment count: chr1 at 1.67M and
        # chrY at 27K alignments measured within 2.2% of each other).
        #
        # This suppresses only the records THIS merge would add; inherited chains
        # still concatenate, so --no-PG alone merely caps growth -- measured on
        # the dirty-input path it would still leave 1,818,752 records / ~792 MB.
        #
        # When inputs come from normalize_cluster_bam they are already collapsed
        # by util/misc/collapse_bam_pg_header.py, and --no-PG alone then holds the
        # result at zero (verified: 13 collapsed BAMs merge to 0 @PG records).
        # But callers may supply `pre_normalized_cluster_bams`, which SKIPS that
        # normalizer entirely -- those chains arrive uncollapsed and concatenate
        # here (84 records in the same 13-way merge with dirty inputs). So this
        # task must be able to collapse too.
        #
        # Done with samtools alone, deliberately: this task runs in a SEPARATELY
        # PINNED container image while this .wdl is read from the local checkout,
        # so invoking a newly added helper script from here would break whenever
        # the two are out of step. `reheader -c` + `-d` verified present in the
        # pinned samtools 1.13.
        samtools merge --no-PG \
            -@ ~{cpu} \
            -o ~{output_basename}.bam \
            ~{sep=' ' normalized_bams}

        # Is there anything to collapse? Counted exactly from the BGZF header,
        # streaming, in bounded memory. Three approaches were rejected, each
        # measured on an 85 MB / 200,000-record header:
        #   samtools view -H  -- appends one record per chain tip, so it cannot
        #                        report its own input (1 for an empty chain, 85
        #                        for an 84-record one, measured in this image),
        #                        and has been seen to hang at this size
        #   read+decode+split -- 191 MB peak RSS (2.2x the header)
        #   pysam to_dict()   -- 481 MB peak RSS (5.6x the header)
        #   streaming (below)  -- 17 MB peak, flat in header size
        # Scaled to the 1,188,154,439-byte header this step exists to remove,
        # the two materializing forms would need roughly 2.7 GB and 6.7 GB --
        # i.e. the check would OOM on exactly the input it is meant to detect.
        #
        # The gate is on the MEASURED header, not on which input route was
        # taken. Provenance would be the wrong discriminator: inputs coming
        # through normalize_cluster_bam are only pre-collapsed if that image
        # contains the collapse, so a newer .wdl running against an older image
        # would skip a chain that is genuinely present. Measuring cannot drift.
        n_pg=$(python3 -c '
import gzip, struct, sys
CHUNK = 1 << 20
with gzip.open(sys.argv[1], "rb") as fh:
    assert fh.read(4) == b"BAM\1", "not a BAM file"
    remaining = struct.unpack("<i", fh.read(4))[0]
    total = 0
    prev = b"\n"
    while remaining > 0:
        chunk = fh.read(min(CHUNK, remaining))
        if not chunk:
            break
        remaining -= len(chunk)
        window = prev + chunk
        total += window.count(b"\n@PG")
        prev = window[-3:]
print(total)
' ~{output_basename}.bam)

        if [ "$n_pg" -eq 0 ]; then
            # The normal route: inputs were already collapsed upstream, and
            # --no-PG held the merge at zero. Skip both the full-file tag scan
            # and the rewrite -- there is nothing to remove, and this BAM is the
            # largest artifact the workflow produces.
            echo "@PG chain already empty in ~{output_basename}.bam; no collapse needed"
        else
            # Dropping @PG definitions is only sound if no alignment references
            # one via a PG:Z: tag. Externally supplied BAMs are exactly where
            # such tags are plausible, so this is checked, never assumed.
            #
            # On a hit this FAILS rather than continuing uncollapsed. That is a
            # deliberate choice against precedent from this very defect: the
            # header cost was a log-visible warning nobody read for months,
            # burning ~27.5 h of a single whole-genome run. A warning is the
            # mechanism already known to fail here. Failing is also not a dead
            # end -- the remedy is named in the message, and
            # util/misc/collapse_bam_pg_header.py --force exists for a caller
            # who has decided the tags are expendable.
            n_pg_tagged=$(samtools view -@ ~{cpu} -c -d PG ~{output_basename}.bam)
            if [ "$n_pg_tagged" -ne 0 ]; then
                echo "ERROR: $n_pg_tagged alignment records in ~{output_basename}.bam" \
                     "carry a PG:Z: tag, so the accumulated @PG header chain ($n_pg" \
                     "records) cannot be collapsed without leaving those" \
                     "references dangling." >&2
                echo "Refusing to emit this BAM: left uncollapsed, per-chromosome" \
                     "region queries against it are dominated by header parsing" \
                     "(~5 min each once the header reached 2.7M records, ~2h07m" \
                     "per cluster), which is the defect this step exists to prevent." >&2
                echo "Remedies: (a) strip the PG:Z: tags upstream, or (b) collapse" \
                     "the inputs yourself with util/misc/collapse_bam_pg_header.py" \
                     "(--force accepts dangling references), then pass them as" \
                     "pre_normalized_cluster_bams." >&2
                exit 2
            fi

            # -P so reheader adds no @PG of its own; -c filters the existing
            # header in-stream, so no full header dump is needed.
            samtools reheader -P -c 'grep -v "^@PG"' ~{output_basename}.bam \
                > ~{output_basename}.pgcollapsed.bam
            mv ~{output_basename}.pgcollapsed.bam ~{output_basename}.bam
            echo "collapsed $n_pg @PG records in ~{output_basename}.bam"
        fi

        # Index AFTER any collapse: rewriting the header shifts every BGZF
        # virtual offset, so an index built before it is invalid against the
        # result and region queries fail with "Invalid BGZF header at offset N".
        samtools index -@ ~{cpu} ~{output_basename}.bam

        echo "Merged ~{length(normalized_bams)} BAMs into ~{output_basename}.bam"
    >>>

    output {
        File merged_bam = "~{output_basename}.bam"
        File merged_bai = "~{output_basename}.bam.bai"
    }

    runtime {
        docker: docker
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} SSD"
    }
}


task emit_shared_chunk_plan {
    # Cut selection ONLY -- no extraction, no chunk directories, no mini BAMs. The
    # single output is the geometry every consumer then applies with --chunk_plan, so
    # all of them cut at identical positions.
    #
    # Called from THREE places, at most once per run: here as the fallback for a
    # standalone final quant, and -- for a whole cluster-guided run -- once before
    # phase 1 from LRAA-cell_cluster_guided.wdl or LRAA-singlecell.wdl, which then
    # thread the file down so the initial pass, per-cluster discovery and the final
    # quant share ONE geometry. Defined here rather than in a subworkflow of its own
    # because those two already import this file, and the task belongs beside the
    # workflow whose comparability argument it exists to serve.
    input {
        File inputBAM
        File referenceGenome
        # The gtf whose models placement must not cut through. OPTIONAL: de novo
        # single-cell has no reference annotation to select against, and the selector
        # treats the annotation constraint as optional -- with none supplied every
        # grid-aligned position in the window is admissible on that axis. See the call
        # site for which annotation belongs here.
        File? annot_gtf
        # Restricts cut selection to the contigs the run will actually process.
        # Empty means every reference the fasta and the bam header agree on.
        #
        # This was omitted at first, on the reasoning that a plan covering extra
        # contigs is harmless because validation is per-PROCESSED-contig -- true
        # for correctness, wrong for cost. MEASURED on a PBMC run restricted to
        # chr21/chr22/chrM: the unrestricted plan selected 475 chunks across every
        # contig in GRCh38 including the _random scaffolds, for a run that
        # processed 11.
        String main_chromosomes = ""
        Boolean HiFi
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        # This workflow does not expose the mapping-quality floors, so the per-cluster
        # LRAA.wdl calls run at ITS defaults and make_chunks forwards those to
        # ChunkedRun explicitly. Cut selection filters on the same resolved floor
        # (ChunkedRun.resolve_min_mapping_quality), and ChunkedRun's own argparse
        # default comes from LRAA_Globals.config rather than from LRAA.wdl -- so
        # omitting these would select cuts at one floor while every consumer extracts
        # at another, and the consuming run refuses that mismatch. Mirrored here; if
        # LRAA.wdl's defaults change or this workflow starts exposing them, they have
        # to be threaded to this task too.
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        Int cpu
        Int memoryGB
        String docker
    }

    # 2x the localized inputs, against make_chunks' 3x: nothing here holds a partition
    # of the BAM. The outputs are the plan plus the per-contig cut artifacts, including
    # a severed-reads BAM carrying only reads that span a candidate cut.
    Float inputsGB = size(inputBAM, "GB") + size(referenceGenome, "GB") + (if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0)
    Float diskRawGB = 2.0 * inputsGB + 50.0
    Int diskGB = if diskRawGB > 100.0 then ceil(diskRawGB) else 100

    command <<<
    set -euo pipefail

    # Symlinked into a writable dir for the same reason make_chunks does it: the input
    # mounts can be read-only, and both `samtools index` and faidx write beside their
    # argument.
    mkdir -p inputs work
    ln -s ~{inputBAM} inputs/input.bam
    if [[ ! -e inputs/input.bam.bai && ! -e inputs/input.bam.csi ]]; then
        samtools index -@ ~{cpu} inputs/input.bam
    fi
    ln -s ~{referenceGenome} inputs/genome.fa
    samtools faidx inputs/genome.fa

    # Emitted through the LRAA DRIVER, not by calling ChunkedRun.py directly.
    # Cut placement filters on the RESOLVED min_per_id, and the --HiFi preset that
    # raises it to 97.0 lives in the driver (LRAA:355) -- the standalone module's
    # resolve_min_per_id just reads LRAA_Globals.config, which nothing has
    # presetted when it runs on its own. Calling the module directly therefore
    # selected at min_per_id 80 under --HiFi while every per-cluster consumer
    # resolved 97.0, and the consumers correctly refused the plan as a different
    # partition. Observed end to end on this fixture. Going through the driver
    # makes the emitting and consuming resolution the SAME CODE, so they cannot
    # drift again -- including for --config_update, which the module would also
    # have ignored.
    #
    # --quant_only WITH an annotation, discovery WITHOUT one. De novo single-cell has
    # no reference annotation to select against, and --quant_only with no --gtf has
    # nothing to quantify and is refused by LRAA. The mode reaches cut selection only
    # through resolve_min_mapping_quality (pylib/ChunkedRun.py:1399-1421) -- quant-only
    # resolves the final-quant floor, discovery the MIN of the two -- and this workflow
    # exposes neither floor, so both are 0 and the two modes resolve to the same
    # effective filter. If it ever does expose them the modes stop agreeing, and the
    # emitting run then has to be in the same mode as its consumers.
    #
    # No --num_total_reads: selection has no use for the TPM denominator, and it is
    # per-cluster regardless.
    #
    # CONTIG-RESTRICTED via --restrict_to_chromosomes when main_chromosomes is set.
    #
    # This block previously said "NO contig restriction, deliberately", on the
    # grounds that the driver takes a single --contig rather than a list and that
    # a plan covering extra contigs is harmless because validation is per
    # PROCESSED contig. The harmlessness is real -- a consumer ignores contigs it
    # does not touch -- but the conclusion was wrong: it paid whole-genome cut
    # selection for a run that needed three contigs. MEASURED on PBMC restricted
    # to chr21/chr22/chrM, the unrestricted plan held 475 chunks over every
    # contig in GRCh38 including _random scaffolds, against 11 processed.
    #
    # The driver DOES accept a list: --restrict_to_chromosomes (space or comma
    # separated) has always existed and is now forwarded to ChunkedRun's own
    # --contigs. So restriction and driver-resolved settings are no longer a
    # choice -- which they were, and why this ran unrestricted.
    /usr/local/src/LRAA/LRAA \
        --bam inputs/input.bam \
        --genome inputs/genome.fa \
        ~{if defined(annot_gtf) then "--gtf " + annot_gtf else ""} \
        ~{if defined(annot_gtf) then "--quant_only" else ""} \
        --output_prefix shared_plan \
        --chunk \
        --chunk_work_dir work \
        --cpu_budget ~{cpu} \
        --emit_cut_plan shared_cut_plan.json \
        ~{if main_chromosomes != "" then "--restrict_to_chromosomes '" + main_chromosomes + "'" else ""} \
        ~{true="--HiFi" false="" HiFi} \
        --min_mapping_quality ~{min_mapping_quality} \
        --min_mapping_quality_for_final_quant ~{min_mapping_quality_for_final_quant} \
        ~{if defined(approx_MB_per_cut) then "--approx_MB_per_cut " + approx_MB_per_cut else ""} \
        ~{if defined(approx_MB_per_cut_wiggle_window) then "--approx_MB_per_cut_wiggle_window " + approx_MB_per_cut_wiggle_window else ""}

    test -s shared_cut_plan.json
    >>>

    output {
        File chunk_plan = "shared_cut_plan.json"
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        # Non-preemptible unlike the per-chunk leaves: this is the one task every
        # cluster job waits on, and a preemption restarts cut selection over the whole
        # pre-partition BAM. Same reasoning as make_chunks in
        # subwdls/LRAA_chunk_scatter.wdl.
        preemptible: 0
        disks: "local-disk ~{diskGB} SSD"
    }
}



task collate_read_assignment_summaries {
    # Cross-invocation collation of LRAA's per-invocation read_assignment summaries.
    # Called once per PHASE, never across phases: an init/pseudobulk pass and a
    # per-cluster pass count the same underlying reads under different normalization,
    # so their totals legitimately disagree and one table holding both would invite a
    # false reconciliation.
    input {
        # One merged read_assignment.summary.tsv per LRAA invocation, chunked or not.
        Array[File] summaryFiles
        # Paired BY POSITION with summaryFiles; the util exits 2 naming both counts if
        # the lengths differ, so a caller cannot silently mislabel rows.
        Array[String] clusterIds
        String outputFilePrefix
        # WHICH read population these counts are over, written as the table's last
        # column. Deliberately REQUIRED here rather than defaulted: the util defaults
        # it to "unspecified" so an unwired caller is honest, and a WDL call site that
        # forgot it should fail the check instead. The label exists because an
        # all_clusters reads_total is a sum over per-cluster populations, each thinned
        # independently by coverage normalization, and is NOT comparable to a
        # pooled-normalized figure -- measured ~95,064 pooled against 198,415 summed
        # over 32 clusters, ~2x apart and both correct. Naming the population makes
        # that non-comparability machine-visible so nobody diffs the two files.
        String population
        # The aggregate row across every input. Turn it OFF for a single-invocation
        # phase, where it would merely restate that one invocation's row.
        Boolean emitAllClusters = true
        String docker
    }

    # Text collation over a few hundred small TSVs; the inputs are the only thing on
    # disk that scales, and the output is smaller than their sum.
    Int diskGB = ceil(size(summaryFiles, "GB") * 3.0) + 10

    command <<<
    set -euo pipefail

    # On PATH in the image via util/misc, so invoked bare.
    collate_read_assignment_summaries.py \
        ~{sep=" " prefix("--summary ", summaryFiles)} \
        ~{sep=" " prefix("--cluster_id ", clusterIds)} \
        --population ~{population} \
        ~{true="" false="--no_all_clusters" emitAllClusters} \
        --output "~{outputFilePrefix}.read_assignment.summary.collated.tsv"
    >>>

    output {
        File collatedSummary = "~{outputFilePrefix}.read_assignment.summary.collated.tsv"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk ~{diskGB} SSD"
    }
}