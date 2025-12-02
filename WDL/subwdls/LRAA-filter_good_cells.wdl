version 1.0

workflow FilterGoodCells {
  input {
    String sample_id
    File gene_sparse_tar_gz
    File? isoform_sparse_tar_gz
    File? splice_pattern_sparse_tar_gz
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    Int memoryGB = 32
    
    # Filter parameters
    Float fdr_threshold = 0.01
    Int? lower_threshold
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
    File gene_sparse_tar_gz
    File? isoform_sparse_tar_gz
    File? splice_pattern_sparse_tar_gz
    String output_prefix
    String docker
    Int memoryGB = 32
    Float fdr_threshold
    Int? lower_threshold
  }

  Int disksize = 50 + ceil(2 * (size(gene_sparse_tar_gz, "GB") + size(isoform_sparse_tar_gz, "GB") + size(splice_pattern_sparse_tar_gz, "GB")))

  command <<<
    set -ex

    # Extract gene sparse matrix into a fixed directory name
    mkdir gene-sparseM-input
    tar -xzf ~{gene_sparse_tar_gz} --strip-components=1 -C gene-sparseM-input

    # Create output directory for filtered gene matrix
    mkdir gene-sparseM-filtered

    # Prepare optional arguments for isoform and splice pattern matrices
    ISOFORM_ARGS=""
    SPLICE_PATTERN_ARGS=""

    # Extract and prepare isoform matrix if provided
    if [ -n "~{isoform_sparse_tar_gz}" ] && [ "~{isoform_sparse_tar_gz}" != "" ]; then
      mkdir isoform-sparseM-input
      tar -xzf ~{isoform_sparse_tar_gz} --strip-components=1 -C isoform-sparseM-input
      mkdir isoform-sparseM-filtered
      ISOFORM_ARGS="--isoform_matrix_dir isoform-sparseM-input --isoform_output_dir isoform-sparseM-filtered"
    fi

    # Extract and prepare splice pattern matrix if provided
    if [ -n "~{splice_pattern_sparse_tar_gz}" ] && [ "~{splice_pattern_sparse_tar_gz}" != "" ]; then
      mkdir splice_pattern-sparseM-input
      tar -xzf ~{splice_pattern_sparse_tar_gz} --strip-components=1 -C splice_pattern-sparseM-input
      mkdir splice_pattern-sparseM-filtered
      SPLICE_PATTERN_ARGS="--splice_pattern_matrix_dir splice_pattern-sparseM-input --splice_pattern_output_dir splice_pattern-sparseM-filtered"
    fi

    # Run the filter_good_cells.R script
    filter_good_cells.R \
      --matrix_dir gene-sparseM-input \
      --output_dir gene-sparseM-filtered \
      --fdr_threshold ~{fdr_threshold} \
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
    docker: docker
    memory: memoryGB + " GB"
    disks: "local-disk " + disksize + " HDD"
  }
}
