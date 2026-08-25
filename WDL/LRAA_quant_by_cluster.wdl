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
        
        Boolean HiFi = false
        # Chunk geometry, forwarded to the per-cluster quantification runs. Unset
        # means LRAA's own defaults (10 Mb spacing, 1 Mb window).
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        String scattering = "by_chromosome"
        String? oversimplify
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int normalize_max_cov_level = 1000
        
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

    # Step 3b: ONE chunk plan, shared by every cluster.
    #
    # Two properties make the clusters comparable. All are quantified against one
    # consolidated GTF, and all build their splice graph from the one normalized BAM
    # produced above -- but cut POSITIONS are otherwise chosen per run. Left alone,
    # each cluster selects cuts on its OWN BAM, gets its own chunk boundaries, and
    # slices that shared normalized BAM at those boundaries, so the splice graph each
    # cluster actually saw differs and their boundary overhang drops differ with it.
    # The second property is then not achieved at all.
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
    if (scattering != "off" && defined(inputBAM)) {
        call emit_shared_chunk_plan {
            input:
                inputBAM = select_first([inputBAM]),
                referenceGenome = referenceGenome,
                # THIS workflow's annot_gtf, which is the SAME gtf the per-cluster calls
                # below quantify against -- for the cluster-guided path that is
                # lraa_merge_gtf_task's consolidated output, not the reference
                # annotation. Load-bearing, not incidental: cut selection treats an
                # annotated transcript model as indivisible and will not cut inside one,
                # so a plan selected against the reference GTF can place a boundary
                # through a NOVEL model that phase-1 discovery found and the merge kept.
                # Every cluster would then quantify a severed model -- a wrong answer
                # rather than a smaller one, and precisely on the novel isoforms
                # discovery exists to find. Do not hoist this task earlier in the graph
                # to where only the reference annotation is available.
                annot_gtf = annot_gtf,
                HiFi = HiFi,
                approx_MB_per_cut = approx_MB_per_cut,
                approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                cpu = cpu_chunk_plan,
                memoryGB = memoryGB_chunk_plan,
                docker = docker
        }
    }

    # Step 4: Quantify each original BAM in parallel using normalized merged BAM for splice graph
    scatter (i in range(length(cluster_bams))) {
        String cluster_sample_id = basename(cluster_bams[i], ".bam")
        
        call LRAA.LRAA_wf as LRAA_quant_cluster {
            input:
                sample_id = cluster_sample_id,
                referenceGenome = referenceGenome,
                inputBAM = cluster_bams[i],
                # The ONE splice graph every cluster is quantified against: cluster bams
                # normalized, merged, then normalized again above. inputBAM stays this
                # cluster's own ORIGINAL bam, so the structure set is global and the reads
                # counted are cluster-local. Supplying it also means LRAA does not
                # normalize again -- LRAA.wdl derives no_norm from it being defined, which
                # is why there is no no_norm argument here any more.
                internal_bam_for_sg = normalize_merged_bam.normalized_bam,
                internal_bam_for_sg_index = normalize_merged_bam.normalized_bai,
                # The ONE chunk geometry every cluster cuts on, selected above from the
                # unthinned pre-partition BAM. Without it each cluster would slice the
                # shared splice-graph BAM at its own boundaries.
                internal_chunk_plan = emit_shared_chunk_plan.chunk_plan,
                annot_gtf = annot_gtf,
                quant_only = true,
                # single-cell: there is nothing to export, since the pre-normalized
                # splice-graph bam above means no normalization runs inside LRAA, and this
                # workflow surfaces its normalized BAMs from Normalize_bam.wdl instead
                retain_normalized_splice_graph_bam = false,
                no_EM = false,
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
                docker = docker
        }
    }

    output {
        # Per-cluster quantification outputs
        Array[File] quant_exprs = LRAA_quant_cluster.mergedQuantExpr
        Array[File] quant_trackings = LRAA_quant_cluster.mergedQuantTracking
        Array[Array[File]] read_assignment_shard_summaries_by_cluster = LRAA_quant_cluster.shardReadAssignmentSummaries
        Array[File] read_assignment_merged_summaries_by_cluster = select_all(LRAA_quant_cluster.mergedReadAssignmentSummary)
        
        # Intermediate outputs (for debugging/reuse)
        Array[File]? partitioned_bams = partition_bam_by_cell_cluster.partitioned_bams
        Array[File] normalized_cluster_bams = normalized_cluster_bams_use
        Array[File] normalized_cluster_bais = normalized_cluster_bais_use
        File merged_bam = merge_bams.merged_bam
        File merged_bai = merge_bams.merged_bai
        File normalized_merged_bam = normalize_merged_bam.normalized_bam
        File normalized_merged_bai = normalize_merged_bam.normalized_bai
        # The geometry every cluster above was cut on. Absent when scattering is off, or
        # when no pre-partition BAM was available to select it from.
        File? shared_chunk_plan = emit_shared_chunk_plan.chunk_plan
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
    # single output is the geometry every per-cluster quant run then applies with
    # --chunk_plan, so all of them cut at identical positions.
    input {
        File inputBAM
        File referenceGenome
        # The consolidated GTF the per-cluster quants target; see the call site for why
        # it must be that one and not the reference annotation.
        File annot_gtf
        # No main_chromosomes: the plan is emitted unrestricted on purpose. See the
        # command block for why that is the safe superset for every consumer.
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
    Float inputsGB = size(inputBAM, "GB") + size(referenceGenome, "GB") + size(annot_gtf, "GB")
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
    # --quant_only with the merged GTF, matching every per-cluster call in this
    # workflow: the mode changes the effective mapping-quality floor cut selection
    # filters on, so the emitting run has to be in the same mode as the consuming
    # ones. No --num_total_reads: selection has no use for the TPM denominator,
    # and it is per-cluster regardless.
    #
    # NO contig restriction, deliberately. The driver takes a single --contig, not
    # a list, and a plan is validated per contig the CONSUMING run processes --
    # extra contigs in the plan are ignored, a missing one is refused. So the
    # unrestricted plan is the safe superset for every consumer, whether it holds
    # one chromosome (by_chromosome) or main_chromosomes (by_chunk). The cost is
    # selection over scaffolds nothing quantifies, once per run.
    /usr/local/src/LRAA/LRAA \
        --bam inputs/input.bam \
        --genome inputs/genome.fa \
        --gtf ~{annot_gtf} \
        --quant_only \
        --output_prefix shared_plan \
        --chunk \
        --chunk_work_dir work \
        --cpu_budget ~{cpu} \
        --emit_cut_plan shared_cut_plan.json \
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
