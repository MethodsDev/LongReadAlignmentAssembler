version 1.0

workflow Incorporate_gene_symbols {
  input {
    String sample_id
    File reference_gtf
    File? final_gtf
    File final_sc_gene_sparse_tar_gz
    File final_sc_isoform_sparse_tar_gz
    File final_sc_splice_pattern_sparse_tar_gz
    File final_sc_gene_transcript_splicehash_mapping
    # Single-nucleus mode's gene-body matrices: gene_full and its spliced / unspliced /
    # ambiguous parts. Same gene ids as the gene matrix, so named the same way.
    File? final_sc_gene_full_sparse_tar_gz
    File? final_sc_gene_spliced_sparse_tar_gz
    File? final_sc_gene_unspliced_sparse_tar_gz
    File? final_sc_gene_ambiguous_sparse_tar_gz
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    Int gffcompare_memoryGB = 8
    Int integrate_memoryGB = 16
  }

  if (defined(final_gtf)) {
    call run_gffcompare {
      input:
        sample_id = sample_id,
        reference_gtf = reference_gtf,
        query_gtf = select_first([final_gtf]),
        docker = docker,
        memoryGB = gffcompare_memoryGB
    }
  }

  call incorporate_gene_symbols_sc as integrate_symbols {
    input:
      sample_id = sample_id,
      reference_gtf = reference_gtf,
      final_gtf = final_gtf,
      gene_sparse_tar_gz = final_sc_gene_sparse_tar_gz,
      isoform_sparse_tar_gz = final_sc_isoform_sparse_tar_gz,
      splice_pattern_sparse_tar_gz = final_sc_splice_pattern_sparse_tar_gz,
      gene_full_sparse_tar_gz = final_sc_gene_full_sparse_tar_gz,
      gene_spliced_sparse_tar_gz = final_sc_gene_spliced_sparse_tar_gz,
      gene_unspliced_sparse_tar_gz = final_sc_gene_unspliced_sparse_tar_gz,
      gene_ambiguous_sparse_tar_gz = final_sc_gene_ambiguous_sparse_tar_gz,
      id_mappings_tsv = final_sc_gene_transcript_splicehash_mapping,
      gffcompare_tracking = run_gffcompare.tracking,
      docker = docker,
      memoryGB = integrate_memoryGB
  }

  output {
    File? gffcompare_tracking = run_gffcompare.tracking
    File? gffcompare_stats = run_gffcompare.stats
    
    File? updated_gtf_with_gene_symbols = integrate_symbols.updated_gtf
    File updated_id_mappings = integrate_symbols.annotated_id_mappings
    File updated_gene_sparse_tar_gz = integrate_symbols.annotated_gene_sparse_tar_gz
    File updated_isoform_sparse_tar_gz = integrate_symbols.annotated_isoform_sparse_tar_gz
    File updated_splice_pattern_sparse_tar_gz = integrate_symbols.annotated_splice_pattern_sparse_tar_gz
    File? updated_gene_full_sparse_tar_gz = integrate_symbols.annotated_gene_full_sparse_tar_gz
    File? updated_gene_spliced_sparse_tar_gz = integrate_symbols.annotated_gene_spliced_sparse_tar_gz
    File? updated_gene_unspliced_sparse_tar_gz = integrate_symbols.annotated_gene_unspliced_sparse_tar_gz
    File? updated_gene_ambiguous_sparse_tar_gz = integrate_symbols.annotated_gene_ambiguous_sparse_tar_gz
  }
}


task run_gffcompare {
  input {
    Int preemptible_tries = 3
    String sample_id
    File reference_gtf
    File query_gtf
    String docker
    Int memoryGB = 8
  }

  Int disksize = 20 + ceil(2 * (size(reference_gtf, "GB") + size(query_gtf, "GB")))
  String output_prefix = "~{sample_id}.gffcmp"
    # Machine shape (v0.47.1, generated text, identical in every N2D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). N2D rather than C3D because the task uses at most
    # 2 cores, and C3D has no 2-core shape, or stays on an HDD disk, and Compute Engine docs
    # appear not to offer pd-standard on C3D (unprobed).
    Int shape_cpu = 2
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
    set -euo pipefail

    # Prepare reference annotation (ensure uncompressed)
    if [[ "~{reference_gtf}" == *.gz ]]; then
      gunzip -c ~{reference_gtf} > reference.gtf
    else
      cp ~{reference_gtf} reference.gtf
    fi

    # Prepare query GTF (ensure uncompressed)
    if [[ "~{query_gtf}" == *.gz ]]; then
      gunzip -c ~{query_gtf} > query.gtf
    else
      cp ~{query_gtf} query.gtf
    fi

    gffcompare -r reference.gtf -o ~{output_prefix} query.gtf > gffcompare.log 2>&1 || {
      echo "gffcompare failed; tailing log" >&2
      tail -n 200 gffcompare.log >&2
      exit 1
    }

    mv ~{output_prefix} ~{output_prefix}.stats
    mv gffcompare.log ~{output_prefix}.gffcompare.log
    
  >>>

  output {
    File tracking = "~{output_prefix}.tracking"
    File stats = "~{output_prefix}.stats"
    File log = "~{output_prefix}.gffcompare.log"
  }

  runtime {

      predefinedMachineType: machine_type

      preemptible: preemptible_tries
    docker: docker
    cpu: 2
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} SSD"
  }
}


task incorporate_gene_symbols_sc {
  input {
    Int preemptible_tries = 3
    String sample_id
    File reference_gtf
    File? final_gtf
    File gene_sparse_tar_gz
    File isoform_sparse_tar_gz
    File splice_pattern_sparse_tar_gz
    File? gene_full_sparse_tar_gz
    File? gene_spliced_sparse_tar_gz
    File? gene_unspliced_sparse_tar_gz
    File? gene_ambiguous_sparse_tar_gz
    File id_mappings_tsv
    File? gffcompare_tracking
    String docker
    Int memoryGB = 16
  }

  Int disksize = 50 + ceil(2 * (size(gene_sparse_tar_gz, "GB") + size(isoform_sparse_tar_gz, "GB") + size(splice_pattern_sparse_tar_gz, "GB")
      + size(gene_full_sparse_tar_gz, "GB") + size(gene_spliced_sparse_tar_gz, "GB") + size(gene_unspliced_sparse_tar_gz, "GB") + size(gene_ambiguous_sparse_tar_gz, "GB")))

  # "." rather than "^": this File crosses task/workflow boundaries and Apptainer's
  # --bind spec parser mis-splits paths containing a literal "^" (Docker is unaffected).
  String gene_sparse_tar_out = "~{sample_id}.withGeneSymbols.gene-sparseM.tar.gz"
  String isoform_sparse_tar_out = "~{sample_id}.withGeneSymbols.isoform-sparseM.tar.gz"
  String splice_sparse_tar_out = "~{sample_id}.withGeneSymbols.splice_pattern-sparseM.tar.gz"
  String updated_gtf_out = "~{sample_id}.withGeneSymbols.gtf"
  String updated_mapping_out = "~{sample_id}.gene_transcript_splicehashcode.withGeneSymbols.tsv"
    # Machine shape (v0.47.1, generated text, identical in every N2D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). N2D rather than C3D because the task uses at most
    # 2 cores, and C3D has no 2-core shape, or stays on an HDD disk, and Compute Engine docs
    # appear not to offer pd-standard on C3D (unprobed).
    Int shape_cpu = 2
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
  set -euo pipefail
  set -x

    # Ensure reference GTF is plain text
    if [[ "~{reference_gtf}" == *.gz ]]; then
      gunzip -c ~{reference_gtf} > reference.gtf
    else
      cp ~{reference_gtf} reference.gtf
    fi

    ~{if defined(final_gtf) then 
      "# Ensure final GTF is plain text and work on a local copy\n" +
      "    if [[ \"~{final_gtf}\" == *.gz ]]; then\n" +
      "      gunzip -c ~{final_gtf} > final.gtf\n" +
      "    else\n" +
      "      cp ~{final_gtf} final.gtf\n" +
      "    fi"
      else ""}

    cp ~{id_mappings_tsv} id_mappings.tsv
    cp ~{gene_sparse_tar_gz} gene_sparse.tar.gz
    cp ~{isoform_sparse_tar_gz} isoform_sparse.tar.gz
    cp ~{splice_pattern_sparse_tar_gz} splice_sparse.tar.gz

    set +o pipefail  # allow tar|head probing without SIGPIPE failures under pipefail
    gene_dir=$(tar -tzf gene_sparse.tar.gz | head -1 | cut -d/ -f1 | sed 's@^\./@@')
    isoform_dir=$(tar -tzf isoform_sparse.tar.gz | head -1 | cut -d/ -f1 | sed 's@^\./@@')
    splice_dir=$(tar -tzf splice_sparse.tar.gz | head -1 | cut -d/ -f1 | sed 's@^\./@@')
    set -o pipefail

    tar -xzf gene_sparse.tar.gz --no-same-owner
    tar -xzf isoform_sparse.tar.gz --no-same-owner
    tar -xzf splice_sparse.tar.gz --no-same-owner

    if [[ -z "${gene_dir}" || -z "${isoform_dir}" || -z "${splice_dir}" ]]; then
      echo "Failed to determine sparse matrix directory names from tarballs" >&2
      exit 1
    fi

    # Gene-body matrices (single-nucleus mode), each named like the gene matrix and
    # re-tarred under its own level name.
    gene_body_dirs=()
    for spec in "gene_full:~{gene_full_sparse_tar_gz}" "gene_spliced:~{gene_spliced_sparse_tar_gz}" \
                "gene_unspliced:~{gene_unspliced_sparse_tar_gz}" "gene_ambiguous:~{gene_ambiguous_sparse_tar_gz}"; do
      level="${spec%%:*}"
      tgz="${spec#*:}"
      [[ -n "${tgz}" ]] || continue
      mkdir "gb_${level}"
      tar -xzf "${tgz}" --no-same-owner -C "gb_${level}"
      d=$(ls -d gb_${level}/*/ | head -1)
      gene_body_dirs+=("${d%/}")
    done

    incorporate_gene_symbols_in_sc_features.py \
      --ref_gtf reference.gtf \
      --id_mappings id_mappings.tsv \
      --sparseM_dirs "${gene_dir}" "${isoform_dir}" "${splice_dir}" ${gene_body_dirs[@]+"${gene_body_dirs[@]}"} \
      ~{if defined(final_gtf) then "--LRAA_gtf final.gtf" else ""} \
      ~{if defined(gffcompare_tracking) then "--gffcompare_tracking " + gffcompare_tracking else ""}

    mv id_mappings.tsv.wAnnotIDs "~{updated_mapping_out}"
    
    # Only move updated GTF if final_gtf was provided
    if ~{defined(final_gtf)}; then
      mv final.gtf.updated.gtf "~{updated_gtf_out}"
    fi

    tar -zcf "~{gene_sparse_tar_out}" "${gene_dir}"
    tar -zcf "~{isoform_sparse_tar_out}" "${isoform_dir}"
    tar -zcf "~{splice_sparse_tar_out}" "${splice_dir}"
    for d in ${gene_body_dirs[@]+"${gene_body_dirs[@]}"}; do
      level="${d%%/*}"
      level="${level#gb_}"
      tar -zcf "~{sample_id}.withGeneSymbols.${level}-sparseM.tar.gz" -C "$(dirname "$d")" "$(basename "$d")"
    done

  >>>

  output {
    File? updated_gtf = "~{updated_gtf_out}"
    File annotated_id_mappings = "~{updated_mapping_out}"
    File annotated_gene_sparse_tar_gz = "~{gene_sparse_tar_out}"
    File annotated_isoform_sparse_tar_gz = "~{isoform_sparse_tar_out}"
    File annotated_splice_pattern_sparse_tar_gz = "~{splice_sparse_tar_out}"
    File? annotated_gene_full_sparse_tar_gz = "~{sample_id}.withGeneSymbols.gene_full-sparseM.tar.gz"
    File? annotated_gene_spliced_sparse_tar_gz = "~{sample_id}.withGeneSymbols.gene_spliced-sparseM.tar.gz"
    File? annotated_gene_unspliced_sparse_tar_gz = "~{sample_id}.withGeneSymbols.gene_unspliced-sparseM.tar.gz"
    File? annotated_gene_ambiguous_sparse_tar_gz = "~{sample_id}.withGeneSymbols.gene_ambiguous-sparseM.tar.gz"
  }

  runtime {

      predefinedMachineType: machine_type

      preemptible: preemptible_tries
    docker: docker
    cpu: 2
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} SSD"
  }
}
