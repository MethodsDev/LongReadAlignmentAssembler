version 1.0

import "Partition_data_by_chromosome.wdl" as PartByChr
import "LRAA_runner.wdl" as LRAA_runner


task mergeResults {
    input {
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        String outputFilePrefix
        String docker
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -eo pipefail

        quant_expr_output="~{outputFilePrefix}.quant.expr"
        for file in ~{sep=' ' quantExprFiles}; do
            if [[ ! -f "$quant_expr_output" ]]; then
                cp "$file" "$quant_expr_output"
            else
                tail -n +2 "$file" >> "$quant_expr_output"
            fi
        done

        quant_tracking_output="~{outputFilePrefix}.quant.tracking"
        for file in ~{sep=' ' quantTrackingFiles}; do
            if [[ ! -f "$quant_tracking_output" ]]; then
                cp "$file" "$quant_tracking_output"
            else
                tail -n +2 "$file" >> "$quant_tracking_output"
            fi
        done
    >>>

    output {
        File mergedQuantExprFile = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}.quant.tracking"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

workflow LRAA_Quant_wf {
     input {
        String sample_id
                 
        File? referenceGenome 
        File? inputBAM
        File? annot_gtf
        
        Array[File]? referenceGenomeArray
        Array[File]? inputBAMArray
        Array[File]? annotGTFArray
         
        String? main_chromosomes # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        
        Int num_total_reads
        Float? min_per_id
        Boolean no_EM
        Int? min_mapping_quality
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        
        
    }

    Boolean run_without_splitting = (!defined(inputBAMArray) && !defined(main_chromosomes) && defined(inputBAM))

    
    if (!run_without_splitting && !defined(inputBAMArray) && defined(inputBAM)) {

        ## Split inputs by main chromosomes

        String nonOptional_main_chromosomes = select_first([main_chromosomes, ""])
        
        call PartByChr.partition_by_chromosome as splitByChr {
            input:
                inputBAM = inputBAM,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                chromosomes_want_partitioned = nonOptional_main_chromosomes,
            
                docker = docker,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
            }
     }

     if (!run_without_splitting) {
        # Scatter computes
         
        Array[File] nonOptionalInputBAMArray = select_first([inputBAMArray, splitByChr.chromosomeBAMs, [] ] )
        Array[File] nonOptionalReferenceGenomeArray = select_first([referenceGenomeArray, splitByChr.chromosomeFASTAs,  [] ] )
        Array[File] annotGTFs = select_first([annotGTFArray, splitByChr.chromosomeGTFs, [] ] )
        
         
        scatter (i in range(length(nonOptionalInputBAMArray))) {
            call LRAA_runner.LRAA_runner as LRAA_quant_scatter {
                input:
                    sample_id = sample_id,
                    shardno = i,
                    inputBAM = nonOptionalInputBAMArray[i],
                    genome_fasta = nonOptionalReferenceGenomeArray[i],
                    annot_gtf = annotGTFs[i],
                    num_total_reads = num_total_reads,
                    min_per_id = min_per_id,
                    quant_only = true,
                    no_norm = true,
                    no_EM = no_EM,
                    numThreads = numThreads,
                    min_mapping_quality = min_mapping_quality,
                    docker = docker,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }

        call mergeResults {
            input:
                quantExprFiles = LRAA_quant_scatter.LRAA_quant_expr,
                quantTrackingFiles = LRAA_quant_scatter.LRAA_quant_tracking,
                outputFilePrefix = sample_id,
                docker = docker,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
    }

    if (run_without_splitting) {

        File nonOptionalInputBAM = select_first([inputBAM, ""])
        File nonOptionalReferenceGenome = select_first([referenceGenome, ""])
            
        call LRAA_runner.LRAA_runner as LRAA_quant {
            input:
                sample_id = sample_id,
                inputBAM = nonOptionalInputBAM,
                genome_fasta = nonOptionalReferenceGenome,
                annot_gtf = annot_gtf,
                num_total_reads = num_total_reads,
                min_per_id = min_per_id,
                quant_only = true,
                no_norm = true,
                no_EM = no_EM,
                numThreads = numThreads,
                min_mapping_quality = min_mapping_quality,
                docker = docker,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
        
    }
    
    output {
        File mergedQuantExpr = select_first([mergeResults.mergedQuantExprFile, LRAA_quant.LRAA_quant_expr, ""]) 
        File mergedQuantTracking = select_first([mergeResults.mergedQuantTrackingFile, LRAA_quant.LRAA_quant_tracking, ""])
    }
}
