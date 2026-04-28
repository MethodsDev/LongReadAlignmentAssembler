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
        String? oversimplify
        Int normalize_max_cov_level = 1000
        
        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
        
        # Chromosome splitting parameters for LRAA quantification
        String main_chromosomes = "" # Set to split by chromosomes, leave empty to run without splitting
        
        Int num_threads_per_worker = 2
        Int num_parallel_contigs = 3
        Int memoryGB_normalize = 8
        Int memoryGB_merge = 16
        Int? memoryGB_quant
        Int? memoryGB_quant_scattered
        Int num_threads_per_worker_scattered = 2
        
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }

    # If bam_files not provided, partition the input BAM by cell clusters
    if (!defined(bam_files)) {
        call PartitionBam.partition_bam_by_cell_cluster {
            input:
                sample_id = sample_id,
                cell_clusters_info = select_first([cell_clusters_info]),
                inputBAM = select_first([inputBAM]),
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
                normalized_bais = pre_normalized_cluster_bais
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
                no_EM = false,
                HiFi = HiFi,
                oversimplify = oversimplify,
                main_chromosomes = main_chromosomes,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                numThreadsPerWorker = num_threads_per_worker,
                numThreadsPerWorkerScattered = num_threads_per_worker_scattered,
                num_parallel_contigs = num_parallel_contigs,
                memoryGB = memoryGB_quant,
                memoryGBPerWorkerScattered = memoryGB_quant_scattered,
                docker = docker
        }
    }

    output {
        # Per-cluster quantification outputs
        Array[File] quant_exprs = LRAA_quant_cluster.mergedQuantExpr
        Array[File] quant_trackings = LRAA_quant_cluster.mergedQuantTracking
        
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
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
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

        samtools merge \
            -@ ~{cpu} \
            -o ~{output_basename}.bam \
            ~{sep=' ' normalized_bams}

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
