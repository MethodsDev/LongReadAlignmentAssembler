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
        Int cpu_scattered = 2
        
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

    # Step 4: Quantify each original BAM in parallel using normalized merged BAM for splice graph
    scatter (i in range(length(cluster_bams))) {
        String cluster_sample_id = basename(cluster_bams[i], ".bam")
        
        call LRAA.LRAA_wf as LRAA_quant_cluster {
            input:
                sample_id = cluster_sample_id,
                referenceGenome = referenceGenome,
                inputBAM = cluster_bams[i],
                bam_for_sg = normalize_merged_bam.normalized_bam,
                annot_gtf = annot_gtf,
                quant_only = true,
                no_norm = true,
                # single-cell: no_norm above already means there is nothing to export, and this
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
