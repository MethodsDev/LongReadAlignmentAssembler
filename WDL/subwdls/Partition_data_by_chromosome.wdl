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

    Float bam_size_gb = if defined(inputBAM) then size(inputBAM, "GB") else 0.0
    Float fasta_size_gb = if defined(genome_fasta) then size(genome_fasta, "GB") else 0.0
    Float gtf_size_gb = if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0
    Float estimated_disk = ceil((bam_size_gb + fasta_size_gb + gtf_size_gb) * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

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
        cpu: samtools_threads
        disks: "local-disk " + disk_gb_int + " SSD"
        memory: "24G"
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
            inputBAM = inputBAM,
            genome_fasta = genome_fasta,
            annot_gtf = annot_gtf,
            chromosomes_want_partitioned = chromosomes_want_partitioned,
            docker = docker,
            samtools_threads = samtools_threads
    }

    output {
        Array[File] chromosomeBAMs = partition_by_chromosome_task.chromosomeBAMs
        Array[File] chromosomeFASTAs = partition_by_chromosome_task.chromosomeFASTAs
        Array[File] chromosomeGTFs = partition_by_chromosome_task.chromosomeGTFs
    }    
}


