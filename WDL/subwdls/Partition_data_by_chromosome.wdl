version 1.0

task partition_by_chromosome_task {
    input {
        File? inputBAM
        File? bam_for_sg
        File? genome_fasta
        File? annot_gtf
        String chromosomes_want_partitioned # ex. "chr1 chr2 chr3 ..."

        String docker

        # Cores the VM is sized for. The task runs alone on its own VM on Terra, so
        # every core it is given is a core it can use: the script spreads its
        # extractions over all of them, bounded by what it can usefully run (one
        # extraction per contig, -@ threads each). The machine is the next C3D size
        # at or above this, and the script is told the cores of THAT machine, not
        # of this request.
        #
        # Size it to the work: the per-cluster call splits one cluster bam (0.5 to
        # 1.8 GB measured) and 4 is plenty; the call that splits the shared inputs
        # (the whole library, or the shared splice-graph bam, plus the genome fasta
        # and annotation) is the wide one.
        Int cpu = 4

        # samtools' -@ for ONE extraction, as a core budget (-@ is this minus one). 1 means
        # no helper threads: each extraction is one process on one core, and the cores of
        # the machine go to MORE extractions at once instead. That is the better use of
        # them once there are many contigs. MEASURED over all 25 contigs of a 6.6 GB bam
        # on 6 cores, level 4 output:
        #     1 worker  x -@ 3   164 s
        #     2 workers x -@ 1   325 s
        #     3 workers x -@ 0   171 s
        #     4 workers x -@ 0   128 s
        # and at the shipped defaults: the shared split (a 6.6 GB splice-graph bam, the 3.1 GB
        # genome fasta and the 1.4 GB annotation) on 16 cores with 14 workers took 54 s at a
        # 974 MB peak, and a per-cluster call (one 518 MB cluster bam) on 4 cores with 2
        # workers took 20 s at 114 MB -- against 230 to 340 s per cluster when every cluster
        # re-split the shared inputs. So one core per extraction scales almost linearly with workers, where helper
        # threads inside one extraction do not (a single extraction's best case on a 1.3 GiB
        # slice was -@ 4 at 1.97 s against 3.81 s at -@ 2, the knee of what one extraction
        # can use). Raise it only for a call with few, huge contigs.
        Int samtools_threads = 1

        # Contigs extracted at once. Unset fills the machine: as many as fit in its
        # cores after the two single-threaded FASTA/GTF jobs that run alongside.
        Int? partition_workers

        # Per-contig BAMs are intermediates read once by a shard. Measured on a
        # 518 MB cluster bam over all 25 contigs at -@ 4: default (6) 518 MB in 22 s,
        # level 1 620 MB in 11 s, level 4 543 MB in 14 s -- 95% of the compression
        # for 64% of the time. See util/partition_data_by_chromosome.py.
        Int bam_compression_level = 4

        # MEASURED peak RSS is small and flat: 129 to 209 MB for the bam extractions at 1
        # to 4 workers on a 6.6 GB bam, and 1.2 to 1.4 GiB for a call that also splits the
        # genome fasta and the annotation (one cluster bam, 7 cores). Region extraction does
        # not hold the bam. This is headroom, not a model of the input.
        Int memoryGB = 4

        # Spot attempts before falling back to a standard VM. Every call is short -- 20 s for
        # a per-cluster split on 4 cores, 54 s for the shared split on 16 -- so a preemption
        # costs little. A whole-library split of a very large bam is the one that could run
        # long; lower this for it if it is measured past ~10 minutes.
        Int preemptible_tries = 3
    }

    # samtools' -@ is ADDITIONAL threads, so N there means N+1 running. Never more
    # than 4 additional for one extraction, because `samtools view` stops scaling
    # there (the measurement above).
    #
    # WDL 1.0 has no min() or max(); both are 1.1 builtins and miniwdl and womtool
    # reject them here, so these are the conditional forms.
    Int samtools_extra_threads_requested = if samtools_threads - 1 < 4 then samtools_threads - 1 else 4

    Float bam_size_gb = if defined(inputBAM) then size(inputBAM, "GB") else 0.0
    Float bam_for_sg_size_gb = if defined(bam_for_sg) then size(bam_for_sg, "GB") else 0.0
    Float fasta_size_gb = if defined(genome_fasta) then size(genome_fasta, "GB") else 0.0
    Float gtf_size_gb = if defined(annot_gtf) then size(annot_gtf, "GB") else 0.0
    Float estimated_disk = ceil((bam_size_gb + bam_for_sg_size_gb + fasta_size_gb + gtf_size_gb) * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu = cpu
    Int c3d_mem = memoryGB
    Int c3d_cpu_tier = if c3d_cpu <= 4 then 4
        else if c3d_cpu <= 8 then 8
        else if c3d_cpu <= 16 then 16
        else if c3d_cpu <= 30 then 30
        else if c3d_cpu <= 60 then 60
        else if c3d_cpu <= 90 then 90
        else if c3d_cpu <= 180 then 180
        else 360
    Int c3d_mem_tier = if c3d_mem <= 32 then 4
        else if c3d_mem <= 64 then 8
        else if c3d_mem <= 128 then 16
        else if c3d_mem <= 240 then 30
        else if c3d_mem <= 480 then 60
        else if c3d_mem <= 720 then 90
        else if c3d_mem <= 1440 then 180
        else 360
    Int c3d_effective_cpu = if c3d_cpu_tier >= c3d_mem_tier then c3d_cpu_tier else c3d_mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int c3d_highcpu_ram = if c3d_effective_cpu == 4 then 8
        else if c3d_effective_cpu == 8 then 16
        else if c3d_effective_cpu == 16 then 32
        else if c3d_effective_cpu == 30 then 59
        else if c3d_effective_cpu == 60 then 118
        else if c3d_effective_cpu == 90 then 177
        else if c3d_effective_cpu == 180 then 354
        else 708
    String c3d_machine_type = if c3d_mem <= c3d_highcpu_ram
        then "c3d-highcpu-${c3d_effective_cpu}"
        else if c3d_mem <= c3d_effective_cpu * 4
        then "c3d-standard-${c3d_effective_cpu}"
        else "c3d-highmem-${c3d_effective_cpu}"

    # One extraction plus the two FASTA/GTF jobs has to fit the machine, or the script
    # refuses its reservation: so the helper threads give way first when the machine is
    # small. WDL 1.0 has no min(), hence the conditional form.
    Int extra_threads_that_fit = if c3d_effective_cpu - 3 < 0 then 0 else c3d_effective_cpu - 3
    Int samtools_extra_threads = if samtools_extra_threads_requested < extra_threads_that_fit then samtools_extra_threads_requested else extra_threads_that_fit

    # Workers that fit the machine: its cores less the 2 for the FASTA/GTF jobs, over
    # the runnable threads of one extraction. Derived from the MACHINE'S cores, so it
    # is computed after the tier above -- the tier takes cpu and memory, never the
    # worker count, which is what keeps this from being circular.
    Int workers_that_fit_raw = (c3d_effective_cpu - 2) / (samtools_extra_threads + 1)
    Int workers_that_fit = if workers_that_fit_raw < 1 then 1 else workers_that_fit_raw
    Int effective_partition_workers = if defined(partition_workers) then (if select_first([partition_workers]) < 1 then 1 else select_first([partition_workers])) else workers_that_fit

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
            --num-workers ~{effective_partition_workers} \
            --reserved-cpu ~{c3d_effective_cpu} \
            --bam-compression-level ~{bam_compression_level} \
            --bam-out-dir split_bams \
            --bam-for-sg-out-dir split_bams_for_sg \
            --fasta-out-dir split_fastas \
            --gtf-out-dir split_gtfs

        # Chromosome name -> file, one tab-separated line per output. The maps below
        # are read from these, so a consumer can look a contig up by NAME instead of
        # trusting that two globbed arrays line up by position.
        name_map() {
            local dir="$1" suffix="$2"
            local f
            for f in "$dir"/*"$suffix"; do
                [ -e "$f" ] || continue
                printf '%s\t%s\n' "$(basename "$f" "$suffix")" "$f"
            done
        }
        name_map split_fastas .genome.fasta > fastas_by_name.tsv
        name_map split_gtfs .annot.gtf > gtfs_by_name.tsv
        ~{if defined(bam_for_sg) then "name_map split_bams_for_sg .bam > sg_bams_by_name.tsv" else ": > sg_bams_by_name.tsv"}
        # A map with nothing to say (no annotation, or no splice-graph bam) must still be a
        # readable map: Cromwell's read_map rejects an EMPTY file ("TSV must be 2 columns to
        # convert to a Map") where miniwdl returns {}, and that failed a whole run after the
        # partition itself had succeeded. The placeholder row points at the file itself so it
        # is a real File to delocalize; nothing looks it up, because every consumer reads
        # these maps only when the matching input (annot_gtf, bam_for_sg) was given.
        for f in fastas_by_name.tsv gtfs_by_name.tsv sg_bams_by_name.tsv; do
            [ -s "$f" ] || printf 'none\t%s\n' "$f" > "$f"
        done
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_bams/*.bam")
        Array[File]? chromosomeBAMsForSG = if defined(bam_for_sg) then glob("split_bams_for_sg/*.bam") else []
        Array[File] chromosomeFASTAs = glob("split_fastas/*.genome.fasta")
        Array[File] chromosomeGTFs = glob("split_gtfs/*.annot.gtf")

        # The same files keyed by chromosome name. Empty for an input that was not
        # supplied.
        Map[String, File] chromosomeFASTAsByName = read_map("fastas_by_name.tsv")
        Map[String, File] chromosomeGTFsByName = read_map("gtfs_by_name.tsv")
        Map[String, File] chromosomeBAMsForSGByName = read_map("sg_bams_by_name.tsv")
    }

    runtime {
        docker: docker
        predefinedMachineType: c3d_machine_type
        bootDiskSizeGb: 50
        preemptible: preemptible_tries
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
        # Forwarded to the task, whose comments carry the measurements.
        Int cpu = 4
        Int samtools_threads = 1
        Int? partition_workers
        Int bam_compression_level = 4
        Int preemptible_tries = 3
    }

    call partition_by_chromosome_task {
        input:
            inputBAM = inputBAM,
            bam_for_sg = bam_for_sg,
            genome_fasta = genome_fasta,
            annot_gtf = annot_gtf,
            chromosomes_want_partitioned = chromosomes_want_partitioned,
            docker = docker,
            cpu = cpu,
            samtools_threads = samtools_threads,
            partition_workers = partition_workers,
            bam_compression_level = bam_compression_level,
            preemptible_tries = preemptible_tries
    }

    output {
        Array[File] chromosomeBAMs = partition_by_chromosome_task.chromosomeBAMs
        Array[File]? chromosomeBAMsForSG = partition_by_chromosome_task.chromosomeBAMsForSG
        Array[File] chromosomeFASTAs = partition_by_chromosome_task.chromosomeFASTAs
        Array[File] chromosomeGTFs = partition_by_chromosome_task.chromosomeGTFs
        Map[String, File] chromosomeFASTAsByName = partition_by_chromosome_task.chromosomeFASTAsByName
        Map[String, File] chromosomeGTFsByName = partition_by_chromosome_task.chromosomeGTFsByName
        Map[String, File] chromosomeBAMsForSGByName = partition_by_chromosome_task.chromosomeBAMsForSGByName
    }    
}


