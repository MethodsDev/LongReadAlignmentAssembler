version 1.0

import "Normalize_bam.wdl" as NormBam
import "../LRAA.wdl" as LRAA

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
        Array[File] bam_files
        
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
        Int memoryGB_quant = 32
        Int memoryGB_quant_scattered = 32
        Int num_threads_per_worker_scattered = 2
        
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }

    # Step 1: Normalize each input BAM in parallel
    scatter (i in range(length(bam_files))) {
        String cluster_label = basename(bam_files[i], ".bam")
        
        call NormBam.normalize_bam_by_strand as normalize_cluster_bam {
            input:
                input_bam = bam_files[i],
                normalize_max_cov_level = normalize_max_cov_level,
                label = cluster_label,
                docker = docker,
                memoryGB = memoryGB_normalize
        }
    }

    # Step 2: Merge all normalized BAMs
    call merge_bams {
        input:
            normalized_bams = normalize_cluster_bam.normalized_bam,
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
    scatter (i in range(length(bam_files))) {
        String cluster_sample_id = basename(bam_files[i], ".bam")
        
        call LRAA.LRAA_wf as LRAA_quant_cluster {
            input:
                sample_id = cluster_sample_id,
                referenceGenome = referenceGenome,
                inputBAM = bam_files[i],
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
        Array[File] normalized_cluster_bams = normalize_cluster_bam.normalized_bam
        Array[File] normalized_cluster_bais = normalize_cluster_bam.normalized_bai
        File merged_bam = merge_bams.merged_bam
        File merged_bai = merge_bams.merged_bai
        File normalized_merged_bam = normalize_merged_bam.normalized_bam
        File normalized_merged_bai = normalize_merged_bam.normalized_bai
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
