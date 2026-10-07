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

    # Neither task can use more than a couple of cores and both are short, so they stay
    # on 2-vCPU N2D machines; memoryGB picks the flavour (highcpu 2, standard 8, highmem 16 GB).
    String n2d_machine_type = if memoryGB <= 2 then "n2d-highcpu-2"
        else if memoryGB <= 8 then "n2d-standard-2"
        else if memoryGB <= 16 then "n2d-highmem-2"
        else if memoryGB <= 32 then "n2d-highmem-4"
        else "n2d-highmem-8"

  command <<<
    set -euo pipefail

    # gffcompare (v0.12.6) cannot read gzipped GTFs, so decompress; a plain-text input
    # is just linked under the working name (no copy). The two decompressions run at
    # the same time.
    prep_gtf() {
      if [[ "$1" == *.gz ]]; then
        gunzip -c "$1" > "$2"
      else
        ln -s "$1" "$2"
      fi
    }
    prep_gtf ~{reference_gtf} reference.gtf &
    ref_pid=$!
    prep_gtf ~{query_gtf} query.gtf &
    query_pid=$!
    wait $ref_pid
    wait $query_pid

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

      preemptible: preemptible_tries
    docker: docker
    predefinedMachineType: n2d_machine_type
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
    File id_mappings_tsv
    File? gffcompare_tracking
    String docker
    Int memoryGB = 16
  }

  Int disksize = 50 + ceil(2 * (size(gene_sparse_tar_gz, "GB") + size(isoform_sparse_tar_gz, "GB") + size(splice_pattern_sparse_tar_gz, "GB")))

  # "." rather than "^": this File crosses task/workflow boundaries and Apptainer's
  # --bind spec parser mis-splits paths containing a literal "^" (Docker is unaffected).
  String gene_sparse_tar_out = "~{sample_id}.withGeneSymbols.gene-sparseM.tar.gz"
  String isoform_sparse_tar_out = "~{sample_id}.withGeneSymbols.isoform-sparseM.tar.gz"
  String splice_sparse_tar_out = "~{sample_id}.withGeneSymbols.splice_pattern-sparseM.tar.gz"
  String updated_gtf_out = "~{sample_id}.withGeneSymbols.gtf"
  String updated_mapping_out = "~{sample_id}.gene_transcript_splicehashcode.withGeneSymbols.tsv"

    # Neither task can use more than a couple of cores and both are short, so they stay
    # on 2-vCPU N2D machines; memoryGB picks the flavour (highcpu 2, standard 8, highmem 16 GB).
    String n2d_machine_type = if memoryGB <= 2 then "n2d-highcpu-2"
        else if memoryGB <= 8 then "n2d-standard-2"
        else if memoryGB <= 16 then "n2d-highmem-2"
        else if memoryGB <= 32 then "n2d-highmem-4"
        else "n2d-highmem-8"

  command <<<
  set -euo pipefail
  set -x

    # The scripts below read plain-text GTFs: decompress a .gz input, or just link a
    # plain one under the working name (no copy). Both decompressions run at the same time.
    prep_gtf() {
      if [[ "$1" == *.gz ]]; then
        gunzip -c "$1" > "$2"
      else
        ln -s "$1" "$2"
      fi
    }
    prep_gtf ~{reference_gtf} reference.gtf &
    ref_pid=$!
    final_gtf_path="~{default="" final_gtf}"
    if [[ -n "${final_gtf_path}" ]]; then
      prep_gtf "${final_gtf_path}" final.gtf &
      final_pid=$!
    fi
    wait $ref_pid
    if [[ -n "${final_gtf_path}" ]]; then wait $final_pid; fi

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

    incorporate_gene_symbols_in_sc_features.py \
      --ref_gtf reference.gtf \
      --id_mappings id_mappings.tsv \
      --sparseM_dirs "${gene_dir}" "${isoform_dir}" "${splice_dir}" \
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

  >>>

  output {
    File? updated_gtf = "~{updated_gtf_out}"
    File annotated_id_mappings = "~{updated_mapping_out}"
    File annotated_gene_sparse_tar_gz = "~{gene_sparse_tar_out}"
    File annotated_isoform_sparse_tar_gz = "~{isoform_sparse_tar_out}"
    File annotated_splice_pattern_sparse_tar_gz = "~{splice_sparse_tar_out}"
  }

  runtime {

      preemptible: preemptible_tries
    docker: docker
    predefinedMachineType: n2d_machine_type
    disks: "local-disk ~{disksize} SSD"
  }
}
