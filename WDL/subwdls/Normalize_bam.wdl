version 1.0

# Normalize_bam.wdl
# Task: run a provided Python normalization script to cap per-base coverage by strand,
# then index the resulting BAM with samtools. 

task normalize_bam_by_strand {
  input {
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
    # BGZF level of the strand bams and per-contig parts, which are read once and deleted;
    # the final bam is never written at it. MEASURED on a 6.6 GB merged bam over 34.5 M
    # records at 8 workers: 545 s at the default against 377 s with this at 1 plus the
    # collapse's scan threaded, records and header identical. See
    # util/normalize_bam_by_strand.py.
    Int intermediate_compression_level = 1
    # Spot attempts before falling back to a standard VM. MEASURED 353 s at 8 workers on a
    # 6.6 GB merged bam (the real merged-bam normalize took 544 s before the speedups), so a
    # preemption costs minutes.
    Int preemptible_tries = 3
  }

  # derive a safe base name in WDL (avoid putting conditional logic inside the command string)
  String base = if label == "" then basename(input_bam) else label

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

  command <<<
set -euo pipefail

# Run normalization script (script is expected in PATH inside the docker image) and index output
normalize_bam_by_strand.py --input_bam "~{input_bam}" --normalize_max_cov_level ~{normalize_max_cov_level} --output_bam "~{base}.norm_~{normalize_max_cov_level}.bam" --num_workers ~{c3d_effective_cpu} --intermediate_compression_level ~{intermediate_compression_level}
samtools index -@ ~{c3d_effective_cpu} "~{base}.norm_~{normalize_max_cov_level}.bam"

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
    docker: docker
    predefinedMachineType: c3d_machine_type
    bootDiskSizeGb: 30
    preemptible: preemptible_tries
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
