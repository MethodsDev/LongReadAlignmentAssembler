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
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    Int cpu = 2
    Int memoryGB = 8
  }

  # derive a safe base name in WDL (avoid putting conditional logic inside the command string)
  String base = if label == "" then basename(input_bam) else label

  command <<<
set -euo pipefail

# Run normalization script (script is expected in PATH inside the docker image) and index output
normalize_bam_by_strand.py --input_bam "~{input_bam}" --normalize_max_cov_level ~{normalize_max_cov_level} --output_bam "~{base}.norm_~{normalize_max_cov_level}.bam"
samtools index "~{base}.norm_~{normalize_max_cov_level}.bam"

echo "WDL: produced ~{base}.norm_~{normalize_max_cov_level}.bam and ~{base}.norm_~{normalize_max_cov_level}.bam.bai"
>>>

  # Compute disk needs: 50 + 4 * size(input_bam in GB)
  Int disksize = 50 + ceil(4 * size(input_bam, "GB"))

  output {
    File normalized_bam = "~{base}.norm_~{normalize_max_cov_level}.bam"
    File normalized_bai = "~{base}.norm_~{normalize_max_cov_level}.bam.bai"
  }

  runtime {
    docker: docker
    bootDiskSizeGb: disksize
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
  }

  call normalize_bam_by_strand {
    input: input_bam=input_bam, normalize_max_cov_level=normalize_max_cov_level, label=label
  }

  output {
    File normalized_bam = normalize_bam_by_strand.normalized_bam
    File normalized_bai = normalize_bam_by_strand.normalized_bai
  }
}
