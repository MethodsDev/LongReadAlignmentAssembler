version 1.0

workflow FilterGoodCells {
  input {
    String sample_id
    File gene_sparse_tar_gz
    File? isoform_sparse_tar_gz
    File? splice_pattern_sparse_tar_gz
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"
    Int memoryGB = 32
    
    # Filter parameters.
    #
    # `seed` seeds the emptyDrops Monte Carlo p-value simulation.  It has a
    # default because an unset seed is what made two runs of one library disagree
    # on 8 of ~15k cells, which Seurat then amplified to 55; a seed that only
    # helps when someone remembers to pass it leaves the pipeline
    # non-reproducible for everyone who does not.  Spelling matches the Seurat
    # subworkflow's `Int seed = 1` -> `--seed`.
    Float fdr_threshold = 0.01
    Int? lower_threshold
    Int seed = 1
  }

  String output_prefix = sample_id + ".genes.filtered"

  call run_filter_good_cells {
    input:
      gene_sparse_tar_gz = gene_sparse_tar_gz,
      isoform_sparse_tar_gz = isoform_sparse_tar_gz,
      splice_pattern_sparse_tar_gz = splice_pattern_sparse_tar_gz,
      output_prefix = output_prefix,
      docker = docker,
      memoryGB = memoryGB,
      fdr_threshold = fdr_threshold,
      seed = seed,
      lower_threshold = lower_threshold
  }

  output {
    File filtered_gene_sparse_tar_gz = run_filter_good_cells.filtered_gene_sparse_tar_gz
    File? filtered_isoform_sparse_tar_gz = run_filter_good_cells.filtered_isoform_sparse_tar_gz
    File? filtered_splice_pattern_sparse_tar_gz = run_filter_good_cells.filtered_splice_pattern_sparse_tar_gz
    File good_cell_barcodes = run_filter_good_cells.good_cell_barcodes
    File filtering_summary = run_filter_good_cells.filtering_summary
  }
}


task run_filter_good_cells {
  input {
      Int preemptible_tries = 3
    File gene_sparse_tar_gz
    File? isoform_sparse_tar_gz
    File? splice_pattern_sparse_tar_gz
    String output_prefix
    String docker
    Int memoryGB = 32
    Float fdr_threshold
    Int seed
    Int? lower_threshold
  }

  Int disksize = 50 + ceil(2 * (size(gene_sparse_tar_gz, "GB") + size(isoform_sparse_tar_gz, "GB") + size(splice_pattern_sparse_tar_gz, "GB")))

    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu = 1
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
    set -ex

    # Extract the input matrices into fixed directory names, all at the same time.
    # Optional inputs are declared here (empty string when absent) so the checks below
    # stay plain shell.
    isoform_tar="~{default="" isoform_sparse_tar_gz}"
    splice_tar="~{default="" splice_pattern_sparse_tar_gz}"
    ISOFORM_ARGS=""
    SPLICE_PATTERN_ARGS=""

    extract() {  # tarball, name
      mkdir "$2-sparseM-input" "$2-sparseM-filtered"
      tar -xzf "$1" --no-same-owner --strip-components=1 -C "$2-sparseM-input"
    }
    extract ~{gene_sparse_tar_gz} gene &
    pids="$!"
    if [ -n "${isoform_tar}" ]; then
      extract "${isoform_tar}" isoform &
      pids="${pids} $!"
      ISOFORM_ARGS="--isoform_matrix_dir isoform-sparseM-input --isoform_output_dir isoform-sparseM-filtered"
    fi
    if [ -n "${splice_tar}" ]; then
      extract "${splice_tar}" splice_pattern &
      pids="${pids} $!"
      SPLICE_PATTERN_ARGS="--splice_pattern_matrix_dir splice_pattern-sparseM-input --splice_pattern_output_dir splice_pattern-sparseM-filtered"
    fi
    for pid in ${pids}; do wait "${pid}"; done

    # Run the filter_good_cells.R script
    filter_good_cells.R \
      --matrix_dir gene-sparseM-input \
      --output_dir gene-sparseM-filtered \
      --fdr_threshold ~{fdr_threshold} \
      --seed ~{seed} \
      ~{if defined(lower_threshold) then "--lower " + lower_threshold else ""} \
      ${ISOFORM_ARGS} \
      ${SPLICE_PATTERN_ARGS} \
      > filter_good_cells.log 2>&1 || {
        echo "filter_good_cells.R failed; tailing log" >&2
        tail -n 200 filter_good_cells.log >&2
        exit 1
      }

    # The matrix files inside are already gzipped, so the outer gzip only has to wrap them
    # (level 1 gives the same size as level 6, +0.02%); one pipeline per matrix, in
    # parallel. The barcodes are extracted at the same time.
    pack() {  # directory, tarball
      tar -cf - "$1" | pigz -1 -p 2 > "$2"
    }
    pids=""
    pack gene-sparseM-filtered ~{output_prefix}.gene-sparseM.tar.gz &
    pids="$!"
    if [ -d isoform-sparseM-filtered ]; then
      pack isoform-sparseM-filtered ~{output_prefix}.isoform-sparseM.tar.gz &
      pids="${pids} $!"
    fi
    if [ -d splice_pattern-sparseM-filtered ]; then
      pack splice_pattern-sparseM-filtered ~{output_prefix}.splice_pattern-sparseM.tar.gz &
      pids="${pids} $!"
    fi
    # Extract the good cell barcodes to a standalone file
    zcat gene-sparseM-filtered/barcodes.tsv.gz > ~{output_prefix}.good_cell_barcodes.txt
    for pid in ${pids}; do wait "${pid}"; done

    # The summary also stays inside the gene tarball (packed above); move it out now.
    mv gene-sparseM-filtered/filtering_summary.tsv ~{output_prefix}.filtering_summary.tsv
  >>>

  output {
    File filtered_gene_sparse_tar_gz = "~{output_prefix}.gene-sparseM.tar.gz"
    File? filtered_isoform_sparse_tar_gz = "~{output_prefix}.isoform-sparseM.tar.gz"
    File? filtered_splice_pattern_sparse_tar_gz = "~{output_prefix}.splice_pattern-sparseM.tar.gz"
    File good_cell_barcodes = "~{output_prefix}.good_cell_barcodes.txt"
    File filtering_summary = "~{output_prefix}.filtering_summary.tsv"
  }

  runtime {

      preemptible: preemptible_tries
    docker: docker
    predefinedMachineType: c3d_machine_type
    disks: "local-disk " + disksize + " SSD"
  }
}
