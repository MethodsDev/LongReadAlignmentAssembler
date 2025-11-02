version 1.0

workflow BuildSparseMatricesFromTracking {
  input {
    String sample_id
    File tracking_file
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    Int memoryGB = 16
  }

  call sc_build_sparse_matrices_from_tracking as build_sc_sparse_matrices {
    input:
      sample_id = sample_id,
      tracking_file = tracking_file,
      docker = docker,
      memoryGB = memoryGB
  }

  output {
    File mapping_file = build_sc_sparse_matrices.mapping_file
    File gene_counts = build_sc_sparse_matrices.gene_counts
    File isoform_counts = build_sc_sparse_matrices.isoform_counts
    File splice_pattern_counts = build_sc_sparse_matrices.splice_pattern_counts

    File gene_sparse_dir_tgz = build_sc_sparse_matrices.gene_sparse_dir_tgz
    File isoform_sparse_dir_tgz = build_sc_sparse_matrices.isoform_sparse_dir_tgz
    File splice_pattern_sparse_dir_tgz = build_sc_sparse_matrices.splice_pattern_sparse_dir_tgz
  }
}


task sc_build_sparse_matrices_from_tracking {
  input {
    String sample_id
    File tracking_file
    String docker
    Int memoryGB = 16
  }

  Int disksize = 50 + ceil(2 * size(tracking_file, "GB"))

  String output_prefix = "~{sample_id}.LRAA.sc"

  command <<<
    set -ex

    singlecell_tracking_to_sparse_matrix.py \
      --tracking ~{tracking_file} \
      --output_prefix ~{output_prefix}

    # Tar the generated sparse matrix directories for compact output
    tar -zcvf "~{output_prefix}^gene-sparseM.tar.gz" "~{output_prefix}^gene-sparseM" || true
    tar -zcvf "~{output_prefix}^isoform-sparseM.tar.gz" "~{output_prefix}^isoform-sparseM" || true
    tar -zcvf "~{output_prefix}^splice_pattern-sparseM.tar.gz" "~{output_prefix}^splice_pattern-sparseM" || true
  >>>

  output {
    File mapping_file = "~{output_prefix}.gene_transcript_splicehashcode.tsv"
    File gene_counts = "~{output_prefix}.gene_cell_counts.tsv"
    File isoform_counts = "~{output_prefix}.isoform_cell_counts.tsv"
    File splice_pattern_counts = "~{output_prefix}.splice_pattern_cell_counts.tsv"

    File gene_sparse_dir_tgz = "~{output_prefix}^gene-sparseM.tar.gz"
    File isoform_sparse_dir_tgz = "~{output_prefix}^isoform-sparseM.tar.gz"
    File splice_pattern_sparse_dir_tgz = "~{output_prefix}^splice_pattern-sparseM.tar.gz"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}
