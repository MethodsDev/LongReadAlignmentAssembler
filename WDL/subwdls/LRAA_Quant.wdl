version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        String docker
        Int numThreads
        File referenceGenome
        File referenceAnnotation_full
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        set -eo pipefail
        mkdir -p split_bams
        
        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index -@ ~{numThreads} ~{inputBAM}
        fi
        
        for chr in ~{main_chromosomes}; do
            samtools view -@ ~{numThreads} -b ~{inputBAM} $chr > split_bams/$chr.bam
            samtools faidx ~{referenceGenome} $chr > split_bams/$chr.genome.fasta
            
            if [ -f "~{referenceAnnotation_full}" ]; then
                cat ~{referenceAnnotation_full} | perl -lane 'if ($F[0] eq "'$chr'") { print; }' > split_bams/$chr.full.annot.gtf
            fi
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_bams/*.genome.fasta")
        Array[File] fullAnnotations = glob("split_bams/*.full.annot.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task splitGTFByChromosome {
    input {
        File referenceAnnotation_full
        String main_chromosomes
        String docker
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -eo pipefail
        mkdir -p split_gtfs

        for chr in ~{main_chromosomes}; do
            cat ~{referenceAnnotation_full} | awk -v chr=$chr '$1 == chr' > split_gtfs/$chr.gtf
        done
    >>>

    output {
        Array[File] chromosomeGTFs = glob("split_gtfs/*.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task lraaPerChromosome {
    input {
        File inputBAM
        File referenceGenome
        String OutDir
        String docker
        Int numThreads
        Boolean? LRAA_no_norm
        Int? LRAA_min_mapping_quality
        File referenceAnnotation_full
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    String min_mapping_quality_flag = "--min_mapping_quality=" + select_first([LRAA_min_mapping_quality, 0])
    
    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir -p ~{OutDir}/Quant_noEM_minMapQ
    
        /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.quant \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 ~{min_mapping_quality_flag} --CPU 1
    >>>
    
    output {
        File lraaQuantExpr = "~{OutDir}/Quant_noEM_minMapQ/LRAA.quant.quant.expr"
        File lraaQuantTracking = "~{OutDir}/Quant_noEM_minMapQ/LRAA.quant.quant.tracking"
    }
    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: "~{numThreads}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task mergeResults {
    input {
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        String outputFilePrefix
        String docker
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        set -eo pipefail

        quant_expr_output="~{outputFilePrefix}_merged_quant.expr"
        for file in ~{sep=' ' quantExprFiles}; do
            if [[ ! -f "$quant_expr_output" ]]; then
                cp "$file" "$quant_expr_output"
            else
                tail -n +2 "$file" >> "$quant_expr_output"
            fi
        done

        quant_tracking_output="~{outputFilePrefix}_merged_quant.tracking"
        for file in ~{sep=' ' quantTrackingFiles}; do
            if [[ ! -f "$quant_tracking_output" ]]; then
                cp "$file" "$quant_tracking_output"
            else
                tail -n +2 "$file" >> "$quant_tracking_output"
            fi
        done
    >>>

    output {
        File mergedQuantExprFile = "~{outputFilePrefix}_merged_quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}_merged_quant.tracking"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

workflow lraaWorkflow {
    input {
        File? inputBAM
        Array[File]? inputBAMArray
        File? referenceGenome
        Array[File]? referenceGenomeArray
        Int? LRAA_min_mapping_quality
        Boolean? LRAA_no_norm
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File referenceAnnotation_full
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    }

    String OutDir = "LRAA_out"

    if (defined(inputBAM)) {
        File nonOptionalInputBAM = select_first([inputBAM, ""])
        File nonOptionalReferenceGenome = select_first([referenceGenome, ""])

        call splitBAMByChromosome {
            input:
                inputBAM = nonOptionalInputBAM,
                main_chromosomes = main_chromosomes,
                docker = docker,
                numThreads = numThreads,
                referenceGenome = nonOptionalReferenceGenome,
                referenceAnnotation_full = referenceAnnotation_full,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }

        scatter (i in range(length(splitBAMByChromosome.chromosomeBAMs))) {
            call lraaPerChromosome {
                input:
                    inputBAM = splitBAMByChromosome.chromosomeBAMs[i],
                    referenceGenome = splitBAMByChromosome.chromosomeFASTAs[i],
                    OutDir = OutDir,
                    docker = docker,
                    numThreads = numThreads,
                    LRAA_no_norm = LRAA_no_norm,
                    LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                    referenceAnnotation_full = splitBAMByChromosome.fullAnnotations[i],
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }
    }

    if (defined(inputBAMArray) && defined(referenceGenomeArray)) {
        Array[File] nonOptionalInputBAMArray = select_first([inputBAMArray, []])
        Array[File] nonOptionalReferenceGenomeArray = select_first([referenceGenomeArray, []])

        call splitGTFByChromosome {
            input:
                referenceAnnotation_full = referenceAnnotation_full,
                main_chromosomes = main_chromosomes,
                docker = docker,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }

        scatter (j in range(length(nonOptionalInputBAMArray))) {
            call lraaPerChromosome as lraaPerChromosomeArray {
                input:
                    inputBAM = nonOptionalInputBAMArray[j],
                    referenceGenome = nonOptionalReferenceGenomeArray[j],
                    OutDir = OutDir,
                    docker = docker,
                    numThreads = numThreads,
                    LRAA_no_norm = LRAA_no_norm,
                    LRAA_min_mapping_quality = LRAA_min_mapping_quality,
                    referenceAnnotation_full = splitGTFByChromosome.chromosomeGTFs[j],
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }
    }

    Array[File] quantExprFiles = if defined(inputBAM) then select_first([lraaPerChromosome.lraaQuantExpr, []]) else select_first([lraaPerChromosomeArray.lraaQuantExpr, []])
    Array[File] quantTrackingFiles = if defined(inputBAM) then select_first([lraaPerChromosome.lraaQuantTracking, []]) else select_first([lraaPerChromosomeArray.lraaQuantTracking, []])
    

    call mergeResults {
        input:
            quantExprFiles = quantExprFiles,
            quantTrackingFiles = quantTrackingFiles,
            outputFilePrefix = "merged",
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    output {
        File mergedQuantExpr = mergeResults.mergedQuantExprFile
        File mergedQuantTracking = mergeResults.mergedQuantTrackingFile
    }
}
