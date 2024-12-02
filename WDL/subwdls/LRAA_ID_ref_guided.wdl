version 1.0

task splitBAMAndGTFByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        File referenceGenome
        String docker
        Int numThreads
        File referenceAnnotation_reduced
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        set -eo pipefail
        mkdir -p split_bams
        mkdir -p split_gtfs
        
        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index -@ ~{numThreads} ~{inputBAM}
        fi
        
        for chr in ~{main_chromosomes}; do
            samtools view -@ ~{numThreads} -b ~{inputBAM} $chr > split_bams/$chr.bam
            samtools faidx ~{referenceGenome} $chr > split_bams/$chr.genome.fasta
            cat ~{referenceAnnotation_reduced} | awk -v chr=$chr '$1 == chr' > split_gtfs/$chr.gtf
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_bams/*.genome.fasta")
        Array[File] chromosomeGTFs = glob("split_gtfs/*.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: "~{numThreads}"
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
        File referenceAnnotation_reduced
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String no_norm_flag = if defined(LRAA_no_norm) && LRAA_no_norm then "--no_norm" else ""
    
    command <<<
        bash ~{monitoringScript} > monitoring.log &

        mkdir -p ~{OutDir}/ID_refguided
    
        /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_refguided/LRAA_refguided \
                                 ~{no_norm_flag} --CPU 1 --gtf ~{referenceAnnotation_reduced} 
    >>>
    
    output {
        File lraaID_refguided_GTF = "~{OutDir}/ID_refguided/LRAA_refguided.gtf"
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
        Array[File] gtfFiles
        String outputFilePrefix
        String docker
        Int memoryGB
        Int diskSizeGB
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        set -eo pipefail

        gtf_output="~{outputFilePrefix}_merged.gtf"
        touch "$gtf_output"

        gtf_files_str="~{sep=' ' gtfFiles}"

        for file in $gtf_files_str; do
            cat "$file" >> "$gtf_output"
        done
    >>>

    output {
        File mergedGtfFile = "~{outputFilePrefix}_merged.gtf"
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
        Array[File]? referenceGenomeArray
        File referenceAnnotation_reduced
        File referenceGenome
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        Boolean? LRAA_no_norm
    }

    String OutDir = "LRAA_out"

    if (defined(inputBAM)) {
        File nonOptionalInputBAM = select_first([inputBAM, ""])

        call splitBAMAndGTFByChromosome {
            input:
                inputBAM = nonOptionalInputBAM,
                main_chromosomes = main_chromosomes,
                docker = docker,
                numThreads = numThreads,
                referenceGenome = referenceGenome,
                referenceAnnotation_reduced = referenceAnnotation_reduced,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }

        scatter (i in range(length(splitBAMAndGTFByChromosome.chromosomeBAMs))) {
            call lraaPerChromosome {
                input:
                    inputBAM = splitBAMAndGTFByChromosome.chromosomeBAMs[i],
                    referenceGenome = splitBAMAndGTFByChromosome.chromosomeFASTAs[i],
                    referenceAnnotation_reduced = splitBAMAndGTFByChromosome.chromosomeGTFs[i],
                    OutDir = OutDir,
                    docker = docker,
                    numThreads = numThreads,
                    LRAA_no_norm = LRAA_no_norm,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }
    }

    if (defined(inputBAMArray) && defined(referenceGenomeArray)) {
        Array[File] nonOptionalInputBAMArray = select_first([inputBAMArray, []])
        Array[File] nonOptionalReferenceGenomeArray = select_first([referenceGenomeArray, []])
        Array[File] nonOptionalChromosomeGTFs = select_first([splitBAMAndGTFByChromosome.chromosomeGTFs, []])

        scatter (j in range(length(nonOptionalInputBAMArray))) {
            call lraaPerChromosome as lraaPerChromosomeArray {
                input:
                    inputBAM = nonOptionalInputBAMArray[j],
                    referenceGenome = nonOptionalReferenceGenomeArray[j],
                    referenceAnnotation_reduced = nonOptionalChromosomeGTFs[j],
                    OutDir = OutDir,
                    docker = docker,
                    numThreads = numThreads,
                    LRAA_no_norm = LRAA_no_norm,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }
    }

    Array[File] nonOptionalGtfFiles = if defined(inputBAM) then select_first([lraaPerChromosome.lraaID_refguided_GTF, []]) else select_first([lraaPerChromosomeArray.lraaID_refguided_GTF, []])

    call mergeResults {
        input:
            gtfFiles = nonOptionalGtfFiles,
            outputFilePrefix = "merged",
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    output {
        File mergedRefguidedGTF = mergeResults.mergedGtfFile
        Array[File]? splitBAMs = if defined(inputBAM) then splitBAMAndGTFByChromosome.chromosomeBAMs else []
        Array[File]? splitFASTAs = if defined(inputBAM) then splitBAMAndGTFByChromosome.chromosomeFASTAs else []
    }
}
