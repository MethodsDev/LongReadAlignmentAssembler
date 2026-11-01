version 1.0

# Normalize_bam.wdl
# Task: run a provided Python normalization script to cap per-base coverage by strand,
# then index the resulting BAM with samtools. 

task normalize_bam_by_strand {
  input {
    Int preemptible_tries = 3
    File input_bam
    Int normalize_max_cov_level
    String label = ""
    # runtime knobs
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    # Also the worker count handed to the script, NOT just the runtime request.
    # A container's cpu quota does not narrow what the kernel reports as
    # available, so a script left to size itself would start one worker per HOST
    # core inside a 2-cpu task. 8 rather than 2 because the work this task does is
    # divisible and was not being divided: the strand split ran as one thread over
    # the whole bam (measured: 15.3 min of a 72.6 min step over 48.1 M records,
    # with 27 of 28 cores idle), and the depth normalization that follows it ran
    # as TWO units, one per strand bam, however many workers were asked for.
    # Both are now one unit per populated reference, so the floor of each is its
    # largest contig -- chr1 at ~10% of a whole-genome bam's records -- and past
    # ~10 workers there is nothing left to give either of them.
    Int cpu = 8
    # Up to `cpu` normalization units are live at once, but a unit now holds ONE
    # contig's int64 depth array (contig_length/depth_window entries) plus that
    # contig's junction tally, where a per-strand unit held one array per contig
    # in the file. MEASURED over the 55 M-record whole-genome strand bam pair of a
    # real cluster-guided run, as summed VmHWM across live workers: 1.05 GiB for
    # the two-way pass (0.52 GiB in its largest worker) against 0.56 GiB for 130
    # per-contig units on 8 workers (0.08 GiB in its largest). Summed VmHWM is a
    # FLOOR -- it adds high-water marks that need not have been simultaneous and
    # is roughly 2.6x under a true cgroup peak -- so 8 GiB stays, unchanged: the
    # measurement says this stage got cheaper, not that the request can shrink.
    Int memoryGB = 8
    # BGZF level of the strand bams and per-contig parts, which are read once and deleted.
    # The final bam is written by samtools merge at its own level. Measured with the threaded
    # PG scan on a 6.6 GB merged bam: 547 s -> 353 s, records, header and size identical.
    Int intermediate_compression_level = 1
  }

  # derive a safe base name in WDL (avoid putting conditional logic inside the command string)
  String base = if label == "" then basename(input_bam) else label
    # Machine shape (v0.47.1, generated text, identical in every SSD task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). The command's thread count stays on the
    # requested cpu even when the shape has more cores. Batch rejects 8 GiB on c3d-highcpu-4
    # (7.5 GiB passes), hence 7 for that tier: such a request steps up to c3d-standard-4.
    Int shape_cpu = cpu
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 30 then 30
        else if shape_cpu <= 60 then 60
        else if shape_cpu <= 90 then 90
        else 180
    Int shape_mem_tier = if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 240 then 30
        else if shape_memGB <= 480 then 60
        else if shape_memGB <= 720 then 90
        else 180
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    Int shape_highcpu_memGB = if shape_cores == 4 then 7
        else if shape_cores == 30 then 59
        else if shape_cores == 60 then 118
        else if shape_cores == 90 then 177
        else if shape_cores == 180 then 354
        else shape_cores * 2
    String machine_type = if shape_memGB <= shape_highcpu_memGB then "c3d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "c3d-standard-~{shape_cores}"
        else "c3d-highmem-~{shape_cores}"


  command <<<
set -euo pipefail

# Run normalization script (script is expected in PATH inside the docker image) and index output
normalize_bam_by_strand.py --input_bam "~{input_bam}" --normalize_max_cov_level ~{normalize_max_cov_level} --output_bam "~{base}.norm_~{normalize_max_cov_level}.bam" --num_workers ~{cpu} --intermediate_compression_level ~{intermediate_compression_level}
samtools index -@ ~{cpu} "~{base}.norm_~{normalize_max_cov_level}.bam"

echo "WDL: produced ~{base}.norm_~{normalize_max_cov_level}.bam and ~{base}.norm_~{normalize_max_cov_level}.bam.bai"
>>>

  # Two stages now write per-contig parts, and neither overlaps the other: the
  # split's parts are gone before the normalization's exist. The high-water mark
  # is the normalization's concatenation -- localized input + two strand bams +
  # every normalized part + the two normalized strand bams -- which is the same
  # shape the merge that follows it already had (input + two strand bams + two
  # normalized bams + merged output), and normalized bams are smaller than the
  # strand bams they were thinned from. MEASURED on a real 5.6 GB cluster-guided
  # input: split peak 3.05x, normalization peak 3.85x, merge peak 3.85x. So ~5x
  # is unchanged and still covers it; every part set is removed as soon as it is
  # concatenated.
  Int disksize = 50 + ceil(5 * size(input_bam, "GB"))

  output {
    File normalized_bam = "~{base}.norm_~{normalize_max_cov_level}.bam"
    File normalized_bai = "~{base}.norm_~{normalize_max_cov_level}.bam.bai"
  }

  runtime {

      predefinedMachineType: machine_type

      preemptible: preemptible_tries
    docker: docker
    bootDiskSizeGb: 30
    cpu: "~{cpu}"
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} SSD"
  }
}

workflow NormalizeBam {
  input {
    File input_bam
    Int normalize_max_cov_level
    String label = ""
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    Int cpu = 8
  }

  call normalize_bam_by_strand {
    input: input_bam=input_bam, normalize_max_cov_level=normalize_max_cov_level, label=label, docker=docker, cpu=cpu
  }

  output {
    File normalized_bam = normalize_bam_by_strand.normalized_bam
    File normalized_bai = normalize_bam_by_strand.normalized_bai
  }
}
