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


# ---------------------------------------------------------------------------
# Per-shard build + streaming merge.
#
# The task above builds the whole library's matrices in one single-threaded
# pass over the merged tracking file: 12.48 h and 114.6 GiB on a 1.5 B-read
# library. These two split that work by shard.
#
# A shard is a (cluster, contig) partition. Every (feature, cell) pair belongs
# to exactly one shard, because feature -> contig and cell -> cluster are both
# functions, so the merge never sums across shards. See
# PLAN.sc_sparse_from_shards.md.
#
# sc_build_shard_sparse runs on the sc image rather than inside
# LRAA_runner_task, because lraa-core carries numpy but not pandas or scipy.
# Called from inside LRAA.wdl's per-contig scatter, so each shard's build starts
# as soon as that shard's runner finishes and overlaps the other shards' work --
# there is no barrier. The alternative, adding pandas and scipy to lraa-base so
# the build can happen in the runner task itself, would avoid re-localizing the
# shard tracking file; it was not taken because it puts those dependencies into
# every published image.
# ---------------------------------------------------------------------------

task sc_build_shard_sparse {
  input {
    String shard_name
    File tracking_file
    String docker
    # PROVISIONAL. Scales with this shard's non-zeros, not the library's: the
    # largest VILLAGE shard (chr1, ~9% of 1.78 B rows) would hold roughly 0.8
    # GiB. Not measured on a real shard.
    Int memoryGB = 8
    Int cpu = 1
  }

  Int disksize = 20 + ceil(3 * size(tracking_file, "GB"))

  command <<<
    set -ex
    mkdir -p shard_out
    tracking_to_shard_sparse.py \
      --tracking ~{tracking_file} \
      --outdir shard_out \
      --shard "~{shard_name}"
    tar -cf "~{shard_name}.sc_shard_sparse.tar" -C shard_out .
  >>>

  output {
    File shard_sparse_tar = "~{shard_name}.sc_shard_sparse.tar"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}


task merge_sc_shard_sparse {
  input {
    String sample_id
    Array[File] shard_sparse_tars
    String docker
    # PROVISIONAL. The merge holds no matrix -- file handles, the barcode
    # remaps and one output row -- so this should be generous by a wide margin.
    # Not measured at library scale.
    Int memoryGB = 16
    Int gzip_level = 1
    # Set for cluster-guided runs, where one feature legitimately comes from
    # several shards (the clusters) carrying disjoint cells. Left false in basic
    # mode so that a feature appearing twice is reported as the error it is.
    Boolean shared_features = false
  }

  Int disksize = 100 + ceil(8 * size(shard_sparse_tars, "GB"))

  String output_prefix = "~{sample_id}.LRAA.sc"

  command <<<
    set -ex

    mkdir -p shards
    i=0
    for t in ~{sep=' ' shard_sparse_tars}; do
      i=$((i+1))
      d="shards/$(printf 's%04d' $i)"
      mkdir -p "$d"
      tar -xf "$t" -C "$d"
    done

    merge_shard_sparse_matrices.py \
      --shard_dirs shards/s* \
      --output_prefix "~{output_prefix}" \
      --gzip_level ~{gzip_level} \
      ~{true="--shared-features" false="" shared_features}

    # Same naming as sc_build_sparse_matrices_from_tracking: the archives use
    # "." where the directories use "^", because Apptainer's --bind spec parser
    # mis-splits paths containing a literal "^".
    for level in gene isoform splice_pattern; do
      if command -v pigz >/dev/null 2>&1; then
        tar --use-compress-program="pigz -~{gzip_level}" \
            -cvf "~{output_prefix}.${level}-sparseM.tar.gz" "~{output_prefix}^${level}-sparseM"
      else
        tar -zcvf "~{output_prefix}.${level}-sparseM.tar.gz" "~{output_prefix}^${level}-sparseM"
      fi
    done
  >>>

  output {
    File mapping_file = "~{output_prefix}.gene_transcript_splicehashcode.tsv"
    File gene_sparse_dir_tgz = "~{output_prefix}.gene-sparseM.tar.gz"
    File isoform_sparse_dir_tgz = "~{output_prefix}.isoform-sparseM.tar.gz"
    File splice_pattern_sparse_dir_tgz = "~{output_prefix}.splice_pattern-sparseM.tar.gz"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}
