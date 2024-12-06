version 1.0

task splitBAMByChromosome {
    input {
        File inputBAM
        String main_chromosomes
        String docker
        Int numThreads
        File referenceGenome
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -eo pipefail
        mkdir -p split_bams
        
        if [ ! -f "~{inputBAM}.bai" ]; then
            samtools index -@ ~{numThreads} ~{inputBAM}
        fi
        
        for chr in ~{main_chromosomes}; do
            samtools view -@ ~{numThreads} -b ~{inputBAM} $chr > split_bams/$chr.bam
            samtools faidx ~{referenceGenome} $chr > split_bams/$chr.genome.fasta
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_bams/*.genome.fasta")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        cpu: "~{numThreads}"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task lraaPerChromosome {
  input {
        String sample_id
        Int shardno
        File inputBAM
        Int num_total_reads
        File referenceGenome
        String OutDir
        String docker
        Int numThreads
        Boolean LRAA_no_norm
        Int memoryGB
        Int diskSizeGB
    }

    String no_norm_flag = if (LRAA_no_norm) then "--no_norm" else ""
    
    command <<<

        mkdir -p ~{OutDir}/ID_reffree
    
        /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reffree/~{sample_id}.shardno-~{shardno}.LRAA_ref-free \
                                 ~{no_norm_flag} --CPU 1 \
                                 --num_total_reads ~{num_total_reads}
    >>>
    
    output {
        File lraaID_reffree_GTF = "~{OutDir}/ID_reffree/~{sample_id}.shardno-~{shardno}.LRAA_ref-free.gtf"
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
    }

    command <<<

        set -eo pipefail

        gtf_output="~{outputFilePrefix}.LRAA_ref-free.gtf"
        touch "$gtf_output"

        gtf_files_str="~{sep=' ' gtfFiles}"

        for file in $gtf_files_str; do
            cat "$file" >> "$gtf_output"
        done
    >>>

    output {
        File mergedGtfFile = "~{outputFilePrefix}.LRAA_ref-free.gtf"
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
        String sample_id
        File? inputBAM
        Array[File]? inputBAMArray
        Array[File]? referenceGenomeArray
        File referenceGenome
        Int num_total_reads
        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        String main_chromosomes = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        Boolean LRAA_no_norm
    }

    String OutDir = "LRAA_out"

    if (defined(inputBAM)) {
        File nonOptionalInputBAM = select_first([inputBAM, ""])

        call splitBAMByChromosome {
            input:
                inputBAM = nonOptionalInputBAM,
                main_chromosomes = main_chromosomes,
                docker = docker,
                numThreads = numThreads,
                referenceGenome = referenceGenome,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }

        scatter (i in range(length(splitBAMByChromosome.chromosomeBAMs))) {
            call lraaPerChromosome {
                 input:
                    sample_id = sample_id,
                    shardno = i,
                    inputBAM = splitBAMByChromosome.chromosomeBAMs[i],
                    referenceGenome = splitBAMByChromosome.chromosomeFASTAs[i],
                    num_total_reads = num_total_reads,
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

        scatter (j in range(length(nonOptionalInputBAMArray))) {
            call lraaPerChromosome as lraaPerChromosomeArray {
                input:
                    sample_id = sample_id,
                    shardno = j,
                    inputBAM = nonOptionalInputBAMArray[j],
                    referenceGenome = nonOptionalReferenceGenomeArray[j],
                    num_total_reads = num_total_reads,
                    OutDir = OutDir,
                    docker = docker,
                    numThreads = numThreads,
                    LRAA_no_norm = LRAA_no_norm,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }
    }

    Array[File] nonOptionalGtfFiles = if defined(inputBAM) then select_first([lraaPerChromosome.lraaID_reffree_GTF, []]) else select_first([lraaPerChromosomeArray.lraaID_reffree_GTF, []])

    call mergeResults {
        input:
            gtfFiles = nonOptionalGtfFiles,
            outputFilePrefix = sample_id,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }
    
    output {
        File mergedReffreeGTF = mergeResults.mergedGtfFile
        Array[File]? splitBAMs = if defined(inputBAM) then splitBAMByChromosome.chromosomeBAMs else []
        Array[File]? splitFASTAs = if defined(inputBAM) then splitBAMByChromosome.chromosomeFASTAs else []
    }
}
