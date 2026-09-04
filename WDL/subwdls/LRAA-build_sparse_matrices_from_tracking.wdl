version 1.0

workflow BuildSparseMatricesFromTracking {
  input {
    String sample_id
    File tracking_file
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"
    Int memoryGB = 16
    String csv_engine = "direct"
    Int gzip_level = 1
    Boolean emit_dense_counts = false
  }

  call sc_build_sparse_matrices_from_tracking as build_sc_sparse_matrices {
    input:
      sample_id = sample_id,
      tracking_file = tracking_file,
      docker = docker,
      memoryGB = memoryGB,
      csv_engine = csv_engine,
      gzip_level = gzip_level,
      emit_dense_counts = emit_dense_counts
  }

  output {
    File mapping_file = build_sc_sparse_matrices.mapping_file
    # Optional: only produced under emit_dense_counts. They are 111 GiB on a
    # 1.5 B-read library, and no workflow output block references them -- they
    # were delocalized and then dropped. See PLAN.sc_sparse_from_shards.md.
    File? gene_counts = build_sc_sparse_matrices.gene_counts
    File? isoform_counts = build_sc_sparse_matrices.isoform_counts
    File? splice_pattern_counts = build_sc_sparse_matrices.splice_pattern_counts

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
    # PROVISIONAL, not measured on this code path. The old build peaked at
    # 114.6 GiB on a 1.5 B-read library, so the previous default of 32 was 3.6x
    # low and a guaranteed OOM on any cgroup-backed backend. 64 is an estimate:
    # a faithful replica of this code measured 54.61 GiB, and removing the
    # explicit sum_duplicates() should take roughly 23 GiB off that -- but the
    # edited script itself has not been run at full scale. Re-derive from a real
    # run before treating this as a sizing rule; a size(tracking_file)-based
    # formula is likely the right end state.
    Int memoryGB = 64
    String csv_engine = "direct"
    Int gzip_level = 1
    Boolean emit_dense_counts = false
  }

  # 148 GB for a 49 GB tracking file. That was under-budgeted while the dense
  # *_cell_counts.tsv files were written (111 GiB of them); with those off the
  # peak occupancy is the localized input plus the uncompressed .mtx files.
  Int disksize = 50 + ceil(2 * size(tracking_file, "GB"))

  String output_prefix = "~{sample_id}.LRAA.sc"

  command <<<
    set -ex

    singlecell_tracking_to_sparse_matrix.py \
      --tracking ~{tracking_file} \
      --output_prefix ~{output_prefix} \
      --csv_engine ~{csv_engine} \
      --gzip_level ~{gzip_level} \
      ~{true="--emit_dense_counts" false="" emit_dense_counts}

    # Tar the generated sparse matrix directories for compact output
    # NOTE: the archive filenames below intentionally use "." rather than the "^"
    # that singlecell_tracking_to_sparse_matrix.py uses for the source directory
    # names: these tar.gz files are passed as File inputs to downstream tasks, and
    # Apptainer's --bind spec parser mis-splits paths containing a literal "^"
    # (Docker's bind mounts are unaffected, but the WDL must work under both).
    if command -v pigz >/dev/null 2>&1; then
      tar --use-compress-program="pigz -~{gzip_level}" -cvf "~{output_prefix}.gene-sparseM.tar.gz" "~{output_prefix}^gene-sparseM" || true
      tar --use-compress-program="pigz -~{gzip_level}" -cvf "~{output_prefix}.isoform-sparseM.tar.gz" "~{output_prefix}^isoform-sparseM" || true
      tar --use-compress-program="pigz -~{gzip_level}" -cvf "~{output_prefix}.splice_pattern-sparseM.tar.gz" "~{output_prefix}^splice_pattern-sparseM" || true
    else
      tar -zcvf "~{output_prefix}.gene-sparseM.tar.gz" "~{output_prefix}^gene-sparseM" || true
      tar -zcvf "~{output_prefix}.isoform-sparseM.tar.gz" "~{output_prefix}^isoform-sparseM" || true
      tar -zcvf "~{output_prefix}.splice_pattern-sparseM.tar.gz" "~{output_prefix}^splice_pattern-sparseM" || true
    fi
  >>>

  output {
    File mapping_file = "~{output_prefix}.gene_transcript_splicehashcode.tsv"
    File? gene_counts = "~{output_prefix}.gene_cell_counts.tsv"
    File? isoform_counts = "~{output_prefix}.isoform_cell_counts.tsv"
    File? splice_pattern_counts = "~{output_prefix}.splice_pattern_cell_counts.tsv"

    File gene_sparse_dir_tgz = "~{output_prefix}.gene-sparseM.tar.gz"
    File isoform_sparse_dir_tgz = "~{output_prefix}.isoform-sparseM.tar.gz"
    File splice_pattern_sparse_dir_tgz = "~{output_prefix}.splice_pattern-sparseM.tar.gz"
  }

  runtime {
    docker: docker
    # The script is single-threaded; only pigz uses more, and only for the
    # output compression at the end.
    cpu: 2
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}
