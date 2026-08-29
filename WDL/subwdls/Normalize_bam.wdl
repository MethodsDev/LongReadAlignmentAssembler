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
    # with 27 of 28 cores idle), and the two depth normalizations are independent
    # files. Past ~10 there is nothing left to give it -- the split's floor is its
    # largest contig, chr1 at ~10% of a whole-genome bam's records.
    Int cpu = 8
    # Both depth normalizations now run at once, so both peaks are live together:
    # each holds one int64 array per contig of contig_length/depth_window entries
    # plus its junction tally, ~0.5 GiB per job for a whole human genome at the
    # default 100 bp window.
    Int memoryGB = 8
  }

  # derive a safe base name in WDL (avoid putting conditional logic inside the command string)
  String base = if label == "" then basename(input_bam) else label

  command <<<
set -euo pipefail

# Run normalization script (script is expected in PATH inside the docker image) and index output
normalize_bam_by_strand.py --input_bam "~{input_bam}" --normalize_max_cov_level ~{normalize_max_cov_level} --output_bam "~{base}.norm_~{normalize_max_cov_level}.bam" --num_workers ~{cpu}
samtools index -@ ~{cpu} "~{base}.norm_~{normalize_max_cov_level}.bam"

echo "WDL: produced ~{base}.norm_~{normalize_max_cov_level}.bam and ~{base}.norm_~{normalize_max_cov_level}.bam.bai"
>>>

  # Live at once on this disk: the localized input, the two strand bams, the
  # per-contig parts the strand split concatenates them from, the two normalized
  # strand bams, and the merged output -- ~5x the input, where it used to be ~4x.
  # The parts are removed as soon as they are concatenated, but they are on disk
  # while the split runs.
  Int disksize = 50 + ceil(5 * size(input_bam, "GB"))

  output {
    File normalized_bam = "~{base}.norm_~{normalize_max_cov_level}.bam"
    File normalized_bai = "~{base}.norm_~{normalize_max_cov_level}.bam.bai"
  }

  runtime {
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
