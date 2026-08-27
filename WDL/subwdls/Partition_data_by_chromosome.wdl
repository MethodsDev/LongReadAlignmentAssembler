version 1.0

task partition_by_chromosome_task {
    input {
        File? inputBAM
        File? bam_for_sg
        File? genome_fasta
        File? annot_gtf
        String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."

        String docker
        # CORE BUDGET for this task, not a thread count -- see samtools_extra_threads
        # below, which is derived from it. The name is kept for callers that already
        # bind it.
        #
        # Was 16. On a 16-core single-worker swarm a cpu:16 reservation cannot be
        # placed until the box is empty, and this task runs once per per-cluster
        # LRAA_wf -- 14 to 32 times in a single-cell run -- for a few seconds of
        # `samtools view` each, so it repeatedly drained the whole box. 5 is where
        # the scaling measurement below stops paying.
        Int samtools_threads = 5
        Int? memoryGB
    }

    # samtools' -@ is ADDITIONAL threads, so N there means N+1 running. Spend the
    # budget minus the main thread, and never more than 4 additional however large
    # a budget a caller passes, because `samtools view` stops scaling there.
    #
    # Measured on a 1.30 GiB slice of a real library (chr1 of a 21.2 GiB mouse
    # isoseqsim BAM), best of two, page cache warm, wall seconds by -@ N:
    #   -@ 2 -> 3.81   -@ 3 -> 2.58   -@ 4 -> 1.97   -@ 5 -> 1.63   -@ 8 -> 1.67
    # Per RESERVED core (N+1) that is 0.66 / 0.77 / 0.81 / 0.81 / 0.53 of linear,
    # so budget 5 (-@ 4) is the knee and 9+ is strictly worse. A budget of 4 would
    # mean -@ 3, which is 31% slower than -@ 4 -- hence 5 and not 4.
    #
    # WDL 1.0 has no min(); it is a 1.1 builtin and both miniwdl and womtool reject
    # it here, so this is the conditional form.
    Int samtools_extra_threads = if samtools_threads - 1 < 4 then samtools_threads - 1 else 4

    Float bam_size_gb = if defined(inputBAM) then size(inputBAM, "GB") else 0.0
    Float bam_for_sg_size_gb = if defined(bam_for_sg) then size(bam_for_sg, "GB") else 0.0
    Float fasta_size_gb = if defined(genome_fasta) then size(genome_fasta, "GB") else 0.0
    Float gtf_size_gb = if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0
    Float estimated_disk = ceil((bam_size_gb + bam_for_sg_size_gb + fasta_size_gb + gtf_size_gb) * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

    # Dynamic memory: 0.5× (BAM + splice-graph BAM + FASTA), floor 24 GiB
    Float mem_raw_partition = 0.5 * (bam_size_gb + bam_for_sg_size_gb + fasta_size_gb)
    Int computed_memoryGB = if mem_raw_partition > 24.0 then ceil(mem_raw_partition) else 24
    Int effective_memoryGB = select_first([memoryGB, computed_memoryGB])

    command <<<
        set -euo pipefail

        export PARTITION_SAMTOOLS_THREADS=~{samtools_extra_threads}

        ulimit -n 8192

        partition_data_by_chromosome.py \
            ~{if defined(inputBAM) then "--input-bam " + inputBAM else ""} \
            ~{if defined(bam_for_sg) then "--bam-for-sg " + bam_for_sg else ""} \
            ~{if defined(genome_fasta) then "--genome-fasta " + genome_fasta else ""} \
            ~{if defined(annot_gtf) then "--annot-gtf " + annot_gtf else ""} \
            --chromosomes ~{chromosomes_want_partitioned} \
            --samtools-threads ~{samtools_extra_threads} \
            --bam-out-dir split_bams \
            --bam-for-sg-out-dir split_bams_for_sg \
            --fasta-out-dir split_fastas \
            --gtf-out-dir split_gtfs
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File]? chromosomeBAMsForSG = if defined(bam_for_sg) then glob("split_bams_for_sg/*.bam") else []
        Array[File] chromosomeFASTAs = glob("split_fastas/*.genome.fasta")
        Array[File] chromosomeGTFs = glob("split_gtfs/*.annot.gtf")
    }

    runtime {
        docker: docker
        bootDiskSizeGb: 50
        cpu: samtools_threads
        memory: effective_memoryGB + " GiB"
        preemptible: 0
        disks: "local-disk " + disk_gb_int + " SSD"
    }
}


workflow partition_by_chromosome {
    input {
        File? inputBAM
        File? bam_for_sg
        File? genome_fasta
        File? annot_gtf
        String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
        # Core budget forwarded to the task; kept in step with the task's own
        # default, whose comment carries the measurement.
        Int samtools_threads = 5
    }

    call partition_by_chromosome_task {
        input:
            inputBAM = inputBAM,
            bam_for_sg = bam_for_sg,
            genome_fasta = genome_fasta,
            annot_gtf = annot_gtf,
            chromosomes_want_partitioned = chromosomes_want_partitioned,
            docker = docker,
            samtools_threads = samtools_threads
    }

    output {
        Array[File] chromosomeBAMs = partition_by_chromosome_task.chromosomeBAMs
        Array[File]? chromosomeBAMsForSG = partition_by_chromosome_task.chromosomeBAMsForSG
        Array[File] chromosomeFASTAs = partition_by_chromosome_task.chromosomeFASTAs
        Array[File] chromosomeGTFs = partition_by_chromosome_task.chromosomeGTFs
    }    
}


