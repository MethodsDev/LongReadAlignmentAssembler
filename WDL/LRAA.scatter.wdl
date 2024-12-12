version 1.0

import "subwdls/LRAA_ID_ref_free.wdl" as IDRefFree
import "subwdls/LRAA_ID_ref_guided.wdl" as IDRefGuided
import "subwdls/LRAA_Quant.wdl" as Quant

workflow CombinedWorkflow {
     input {
        String sample_id
        File inputBAM
        File referenceGenome
        File? referenceGTF
        String mode
        String dockerImage = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Int numThreads = 8
        Int memoryGB = 32
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        Boolean LRAA_no_norm = false
        Int LRAA_min_mapping_quality = 20
        Float min_per_id = 98.0
        Boolean no_EM = false
        Float min_isoform_fraction = 0.01
        Float min_monoexonic_TPM = 1.0
        Boolean no_filter_internal_priming = false
        Float min_alt_splice_freq = 0.01
        Float min_alt_unspliced_freq = 0.01
    }

    Int diskSizeGB = 128

    String outputDir = "LRAA_out"

    
    call count_bam {
      input:
        bam = inputBAM
    }
        
    if (mode == "ID_ref_free_Quant_mode") {

        call IDRefFree.lraaWorkflow as IDRefFreeWorkflow {
          input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                num_total_reads = count_bam.count,
                min_per_id = min_per_id,
                no_EM = no_EM,
                min_isoform_fraction = min_isoform_fraction,
                min_monoexonic_TPM = min_monoexonic_TPM,
                no_filter_internal_priming = no_filter_internal_priming,
                min_alt_splice_freq = min_alt_splice_freq,
                min_alt_unspliced_freq = min_alt_unspliced_freq
          
        }

        call Quant.lraaWorkflow as QuantFree2Workflow {
             input:
                sample_id = sample_id,
                inputBAMArray = IDRefFreeWorkflow.splitBAMs,
                referenceGenomeArray = IDRefFreeWorkflow.splitFASTAs,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = IDRefFreeWorkflow.mergedReffreeGTF,
                main_chromosomes = main_chromosomes,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                num_total_reads = count_bam.count,
                min_per_id = min_per_id,
                no_EM = no_EM
        }
    }

    if (mode == "ID_ref_guided_Quant_mode") {

        File guaranteedRefGuided = select_first([referenceGTF])

        call IDRefGuided.lraaWorkflow as IDRefGuidedWorkflow {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                referenceAnnotation_reduced = guaranteedRefGuided,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                main_chromosomes = main_chromosomes,
                LRAA_no_norm = LRAA_no_norm,
                num_total_reads = count_bam.count,
                min_per_id =  min_per_id,
                no_EM = no_EM,
                min_isoform_fraction = min_isoform_fraction,
                min_monoexonic_TPM = min_monoexonic_TPM,
                no_filter_internal_priming = no_filter_internal_priming,
                min_alt_splice_freq = min_alt_splice_freq,
                min_alt_unspliced_freq = min_alt_unspliced_freq
        }
        
        call Quant.lraaWorkflow as QuantGuided2Workflow {
          input:
                sample_id = sample_id,
                inputBAMArray = IDRefGuidedWorkflow.splitBAMs,
                referenceGenomeArray = IDRefGuidedWorkflow.splitFASTAs,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = IDRefGuidedWorkflow.mergedRefguidedGTF,
                main_chromosomes = main_chromosomes,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                num_total_reads = count_bam.count,
                min_per_id = min_per_id,
                no_EM = no_EM
        }
    }

    if (mode == "Quant_only_mode") {

        File guaranteedRefQuantOnly = select_first([referenceGTF])

        call Quant.lraaWorkflow as QuantOnlyWorkflow {
             input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                referenceGenome = referenceGenome,
                numThreads = numThreads,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = dockerImage,
                referenceAnnotation_full = guaranteedRefQuantOnly,
                main_chromosomes = main_chromosomes,
                LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                num_total_reads = count_bam.count,
                min_per_id = min_per_id,
                no_EM = no_EM
        }
    }

    output {
        File? UpdatedGTF = if (mode == "ID_ref_free_Quant_mode") then IDRefFreeWorkflow.mergedReffreeGTF  else if (mode == "ID_ref_guided_Quant_mode") then IDRefGuidedWorkflow.mergedRefguidedGTF else "null"
        File? Quant = if (mode == "ID_ref_free_Quant_mode") then QuantFree2Workflow.mergedQuantExpr else if (mode == "ID_ref_guided_Quant_mode") then QuantGuided2Workflow.mergedQuantExpr else QuantOnlyWorkflow.mergedQuantExpr
        File? Tracking = if (mode == "ID_ref_free_Quant_mode") then QuantFree2Workflow.mergedQuantTracking else if (mode == "ID_ref_guided_Quant_mode") then QuantGuided2Workflow.mergedQuantTracking else QuantOnlyWorkflow.mergedQuantTracking
    }
}




task count_bam {
  input {
    File bam
  }

  command <<<
    set -ex
    samtools view -c ~{bam}

  >>>
  runtime {
    docker: "quay.io/ucsc_cgl/samtools"
    disks: "local-disk " + ceil(2 * size(bam, "GB") ) + " HDD"
    memory: "4G"
  }
  output {
    Int count = read_int(stdout())
  }
}
