version 1.0

task partition_by_chromosome_task {
    input {
        File? inputBAM
        File? genome_fasta
    File? annot_gtf
    String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."

    String docker
    Int samtools_threads = 16
    }

    command <<<
        set -euo pipefail

        export PARTITION_SAMTOOLS_THREADS=~{samtools_threads}

        ulimit -n 8192

        partition_data_by_chromosome.py \
            ~{if defined(inputBAM) then "--input-bam " + inputBAM else ""} \
            ~{if defined(genome_fasta) then "--genome-fasta " + genome_fasta else ""} \
            ~{if defined(annot_gtf) then "--annot-gtf " + annot_gtf else ""} \
            --chromosomes ~{chromosomes_want_partitioned} \
            --samtools-threads ~{samtools_threads} \
            --bam-out-dir split_bams \
            --fasta-out-dir split_fastas \
            --gtf-out-dir split_gtfs
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File] chromosomeFASTAs = glob("split_fastas/*.genome.fasta")
        Array[File] chromosomeGTFs = glob("split_gtfs/*.annot.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 50
        cpu: samtools_threads
        memory: "24 GiB"
        preemptible: 0
        disks: "local-disk " + (if ceil((size(inputBAM, "GB") + size(genome_fasta, "GB") + size(annot_gtf, "GB")) * 2.2 + 20) > 150 then ceil((size(inputBAM, "GB") + size(genome_fasta, "GB") + size(annot_gtf, "GB")) * 2.2 + 20) else 150) + " SSD"
    }
}


workflow partition_by_chromosome {
  input {
        File? inputBAM
        File? genome_fasta
        File? annot_gtf
        String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
            Int samtools_threads = 16
    }

    call partition_by_chromosome_task {
        input:
          inputBAM=inputBAM,
          genome_fasta=genome_fasta,
          annot_gtf=annot_gtf,
          chromosomes_want_partitioned=chromosomes_want_partitioned,
                    docker = docker,
                    samtools_threads = samtools_threads

      
    }

    output {
        Array[File] chromosomeBAMs = partition_by_chromosome_task.chromosomeBAMs
        Array[File] chromosomeFASTAs = partition_by_chromosome_task.chromosomeFASTAs
        Array[File] chromosomeGTFs = partition_by_chromosome_task.chromosomeGTFs
    }    
}


