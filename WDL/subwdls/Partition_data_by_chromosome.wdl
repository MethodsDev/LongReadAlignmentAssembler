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

        # Contigs extracted CONCURRENTLY. Default 1 is the historical serial pass,
        # and it is the default for a reason: this task also runs once per cluster --
        # 14 to 32 times in a single-cell run, seconds of work each -- and that is
        # exactly the case a wide reservation ruins, which is why cpu was cut to 5
        # above. Only the caller that owns the ONE big top-level partition should
        # raise this.
        #
        # MEASURED without it, on a 188 GB library (1.51 B mapped reads, 25 contigs)
        # on a 28-core host: partition_by_chromosome_task ran 27+ minutes with ONE
        # core busy, because the script's four-way pool fans out over JOB TYPES --
        # BAM, BAM_FOR_SG, FASTA, GTF -- and three of the four are trivial or absent,
        # so the whole task is a serial pass over the largest bam in the pipeline.
        # It precedes all shard work, on an otherwise idle box.
        Int partition_workers = 1
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

    # What the script can actually run at once: partition_workers concurrent
    # `samtools view`, each needing its -@ count PLUS its own main thread, plus the
    # two single-threaded FASTA/GTF jobs that run alongside them. Reserving less
    # oversubscribes -- and under unprivileged Apptainer there is no cgroup, so a
    # reservation does not CAP anything and an oversubscribed task competes with its
    # neighbours for real. Reserving more repeats the mistake that cut this task to 5:
    # cpu nobody uses cannot be placed until the box drains.
    #
    # At the default partition_workers = 1 this is samtools_threads + 2, one core more
    # than the historical 5, which is honest rather than new: the FASTA and GTF jobs
    # always ran alongside the bam pass and were never counted.
    Int partition_cpu = partition_workers * (samtools_extra_threads + 1) + 2

    Float bam_size_gb = if defined(inputBAM) then size(inputBAM, "GB") else 0.0
    Float bam_for_sg_size_gb = if defined(bam_for_sg) then size(bam_for_sg, "GB") else 0.0
    Float fasta_size_gb = if defined(genome_fasta) then size(genome_fasta, "GB") else 0.0
    Float gtf_size_gb = if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0
    Float estimated_disk = ceil((bam_size_gb + bam_for_sg_size_gb + fasta_size_gb + gtf_size_gb) * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

    # Dynamic memory: 0.5x (BAM + splice-graph BAM + FASTA), floor 24 GiB.
    #
    # Concurrency adds to this, because each `samtools view` holds its own read and
    # BGZF buffers: partition_workers of them run at once. 1 GiB per additional worker
    # is a deliberately loose allowance -- samtools' per-process footprint at -@ 4 is
    # far below that -- chosen because under unprivileged Apptainer there is NO cgroup,
    # so this reservation caps nothing and an underestimate becomes a host OOM rather
    # than a task failure. UNMEASURED at concurrency above 1; the peak wants recording
    # from a real wide run before the allowance is tightened or a default is raised.
    Float mem_raw_partition = 0.5 * (bam_size_gb + bam_for_sg_size_gb + fasta_size_gb)
    Int computed_memoryGB = if mem_raw_partition > 24.0 then ceil(mem_raw_partition) else 24
    Int concurrency_memoryGB = computed_memoryGB + (partition_workers - 1)
    Int effective_memoryGB = select_first([memoryGB, concurrency_memoryGB])

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
            --num-workers ~{partition_workers} \
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
        cpu: partition_cpu
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
        # Contigs extracted concurrently. Default 1 is the historical serial pass; see
        # the task's own comment for why the default is not higher (this task also runs
        # once per cluster, 14 to 32 times per single-cell run). The task's cpu
        # reservation is derived from this, so raising it here raises the reservation
        # rather than oversubscribing it.
        Int partition_workers = 1
    }

    call partition_by_chromosome_task {
        input:
            inputBAM = inputBAM,
            bam_for_sg = bam_for_sg,
            genome_fasta = genome_fasta,
            annot_gtf = annot_gtf,
            chromosomes_want_partitioned = chromosomes_want_partitioned,
            docker = docker,
            samtools_threads = samtools_threads,
            partition_workers = partition_workers
    }

    output {
        Array[File] chromosomeBAMs = partition_by_chromosome_task.chromosomeBAMs
        Array[File]? chromosomeBAMsForSG = partition_by_chromosome_task.chromosomeBAMsForSG
        Array[File] chromosomeFASTAs = partition_by_chromosome_task.chromosomeFASTAs
        Array[File] chromosomeGTFs = partition_by_chromosome_task.chromosomeGTFs
    }    
}


