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
    # Machine shape (v0.47.1, generated text, identical in every small or HDD task). cpu: and
    # memory: stay as they are, so miniwdl's reservations do not change; this only names the
    # Terra VM (miniwdl ignores predefinedMachineType). N2D rather than C3D here because the
    # task uses at most 2 cores or an HDD disk: C3D has no 2-core shape, and Compute Engine
    # docs appear not to offer pd-standard on C3D (unprobed). Disks stay as they were until a
    # probe settles that.
    Int shape_cpu = 1
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 2 then 2
        else if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 32 then 32
        else if shape_cpu <= 48 then 48
        else if shape_cpu <= 64 then 64
        else if shape_cpu <= 80 then 80
        else 96
    Int shape_mem_tier = if shape_memGB <= 16 then 2
        else if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 256 then 32
        else if shape_memGB <= 384 then 48
        else if shape_memGB <= 512 then 64
        else if shape_memGB <= 640 then 80
        else 96
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    String machine_type = if shape_memGB <= shape_cores then "n2d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "n2d-standard-~{shape_cores}"
        else "n2d-highmem-~{shape_cores}"


  command <<<
    set -ex

    # Extract gene sparse matrix into a fixed directory name
    mkdir gene-sparseM-input
    tar -xzf ~{gene_sparse_tar_gz} --no-same-owner --strip-components=1 -C gene-sparseM-input

    # Create output directory for filtered gene matrix
    mkdir gene-sparseM-filtered

    # Prepare optional arguments for isoform and splice pattern matrices
    ISOFORM_ARGS=""
    SPLICE_PATTERN_ARGS=""

    # Extract and prepare isoform matrix if provided
    if [ -n "~{isoform_sparse_tar_gz}" ] && [ "~{isoform_sparse_tar_gz}" != "" ]; then
      mkdir isoform-sparseM-input
      tar -xzf ~{isoform_sparse_tar_gz} --no-same-owner --strip-components=1 -C isoform-sparseM-input
      mkdir isoform-sparseM-filtered
      ISOFORM_ARGS="--isoform_matrix_dir isoform-sparseM-input --isoform_output_dir isoform-sparseM-filtered"
    fi

    # Extract and prepare splice pattern matrix if provided
    if [ -n "~{splice_pattern_sparse_tar_gz}" ] && [ "~{splice_pattern_sparse_tar_gz}" != "" ]; then
      mkdir splice_pattern-sparseM-input
      tar -xzf ~{splice_pattern_sparse_tar_gz} --no-same-owner --strip-components=1 -C splice_pattern-sparseM-input
      mkdir splice_pattern-sparseM-filtered
      SPLICE_PATTERN_ARGS="--splice_pattern_matrix_dir splice_pattern-sparseM-input --splice_pattern_output_dir splice_pattern-sparseM-filtered"
    fi

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

    # Extract the good cell barcodes to a standalone file
    zcat gene-sparseM-filtered/barcodes.tsv.gz > ~{output_prefix}.good_cell_barcodes.txt

    # Tar and compress the filtered sparse matrix directory
    tar -czf ~{output_prefix}.gene-sparseM.tar.gz gene-sparseM-filtered

    # Copy the summary file to the expected output location
    cp gene-sparseM-filtered/filtering_summary.tsv ~{output_prefix}.filtering_summary.tsv

    # Tar and compress the filtered isoform matrix if it was processed
    if [ -d isoform-sparseM-filtered ]; then
      tar -czf ~{output_prefix}.isoform-sparseM.tar.gz isoform-sparseM-filtered
    fi

    # Tar and compress the filtered splice pattern matrix if it was processed
    if [ -d splice_pattern-sparseM-filtered ]; then
      tar -czf ~{output_prefix}.splice_pattern-sparseM.tar.gz splice_pattern-sparseM-filtered
    fi
  >>>

  output {
    File filtered_gene_sparse_tar_gz = "~{output_prefix}.gene-sparseM.tar.gz"
    File? filtered_isoform_sparse_tar_gz = "~{output_prefix}.isoform-sparseM.tar.gz"
    File? filtered_splice_pattern_sparse_tar_gz = "~{output_prefix}.splice_pattern-sparseM.tar.gz"
    File good_cell_barcodes = "~{output_prefix}.good_cell_barcodes.txt"
    File filtering_summary = "~{output_prefix}.filtering_summary.tsv"
  }

  runtime {

      predefinedMachineType: machine_type

      preemptible: preemptible_tries
    docker: docker
    memory: memoryGB + " GB"
    disks: "local-disk " + disksize + " HDD"
  }
}
