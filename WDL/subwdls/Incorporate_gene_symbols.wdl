version 1.0

workflow Incorporate_gene_symbols {
  input {
    String sample_id
    File reference_gtf
    File final_gtf
    File final_sc_gene_sparse_tar_gz
    File final_sc_isoform_sparse_tar_gz
    File final_sc_splice_pattern_sparse_tar_gz
    File final_sc_gene_transcript_splicehash_mapping
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    Int gffcompare_memoryGB = 8
    Int integrate_memoryGB = 16
  }

  call run_gffcompare {
    input:
      sample_id = sample_id,
      reference_gtf = reference_gtf,
      query_gtf = final_gtf,
      docker = docker,
      memoryGB = gffcompare_memoryGB
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
    File gffcompare_tracking = run_gffcompare.tracking
    File gffcompare_stats = run_gffcompare.stats
    
    File updated_gtf_with_gene_symbols = integrate_symbols.updated_gtf
    File updated_id_mappings = integrate_symbols.annotated_id_mappings
    File updated_gene_sparse_tar_gz = integrate_symbols.annotated_gene_sparse_tar_gz
    File updated_isoform_sparse_tar_gz = integrate_symbols.annotated_isoform_sparse_tar_gz
    File updated_splice_pattern_sparse_tar_gz = integrate_symbols.annotated_splice_pattern_sparse_tar_gz
  }
}


task run_gffcompare {
  input {
    String sample_id
    File reference_gtf
    File query_gtf
    String docker
    Int memoryGB = 8
  }

  Int disksize = 20 + ceil(2 * (size(reference_gtf, "GB") + size(query_gtf, "GB")))
  String output_prefix = "~{sample_id}.gffcmp"

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
  >>>

  output {
    File tracking = "~{output_prefix}.tracking"
    File stats = "~{output_prefix}.stats"
    File log = "gffcompare.log"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}


task incorporate_gene_symbols_sc {
  input {
    String sample_id
    File reference_gtf
    File final_gtf
    File gene_sparse_tar_gz
    File isoform_sparse_tar_gz
    File splice_pattern_sparse_tar_gz
    File id_mappings_tsv
    File gffcompare_tracking
    String docker
    Int memoryGB = 16
  }

  Int disksize = 50 + ceil(2 * (size(gene_sparse_tar_gz, "GB") + size(isoform_sparse_tar_gz, "GB") + size(splice_pattern_sparse_tar_gz, "GB")))

  String gene_sparse_tar_out = "~{sample_id}.withGeneSymbols^gene-sparseM.tar.gz"
  String isoform_sparse_tar_out = "~{sample_id}.withGeneSymbols^isoform-sparseM.tar.gz"
  String splice_sparse_tar_out = "~{sample_id}.withGeneSymbols^splice_pattern-sparseM.tar.gz"
  String updated_gtf_out = "~{sample_id}.withGeneSymbols.gtf"
  String updated_mapping_out = "~{sample_id}.gene_transcript_splicehashcode.withGeneSymbols.tsv"

  command <<<
    set -euo pipefail

    # Ensure reference GTF is plain text
    if [[ "~{reference_gtf}" == *.gz ]]; then
      gunzip -c ~{reference_gtf} > reference.gtf
    else
      cp ~{reference_gtf} reference.gtf
    fi

    # Ensure final GTF is plain text and work on a local copy
    if [[ "~{final_gtf}" == *.gz ]]; then
      gunzip -c ~{final_gtf} > final.gtf
    else
      cp ~{final_gtf} final.gtf
    fi

    cp ~{id_mappings_tsv} id_mappings.tsv
    cp ~{gene_sparse_tar_gz} gene_sparse.tar.gz
    cp ~{isoform_sparse_tar_gz} isoform_sparse.tar.gz
    cp ~{splice_pattern_sparse_tar_gz} splice_sparse.tar.gz

    gene_dir=$(tar -tzf gene_sparse.tar.gz | head -1 | cut -d/ -f1 | sed 's@^\./@@')
    isoform_dir=$(tar -tzf isoform_sparse.tar.gz | head -1 | cut -d/ -f1 | sed 's@^\./@@')
    splice_dir=$(tar -tzf splice_sparse.tar.gz | head -1 | cut -d/ -f1 | sed 's@^\./@@')

    tar -xzf gene_sparse.tar.gz
    tar -xzf isoform_sparse.tar.gz
    tar -xzf splice_sparse.tar.gz

    if [[ -z "${gene_dir}" || -z "${isoform_dir}" || -z "${splice_dir}" ]]; then
      echo "Failed to determine sparse matrix directory names from tarballs" >&2
      exit 1
    fi

    incorporate_gene_symbols_in_sc_features.py  \
      --ref_gtf reference.gtf \
      --id_mappings id_mappings.tsv \
      --sparseM_dirs "${gene_dir}" "${isoform_dir}" "${splice_dir}" \
      --LRAA_gtf final.gtf \
      --gffcompare_tracking ~{gffcompare_tracking}

    mv id_mappings.tsv.wAnnotIDs "${updated_mapping_out}"
    mv final.gtf.updated.gtf "${updated_gtf_out}"

    tar -zcf "${gene_sparse_tar_out}" "${gene_dir}"
    tar -zcf "${isoform_sparse_tar_out}" "${isoform_dir}"
    tar -zcf "${splice_sparse_tar_out}" "${splice_dir}"
  >>>

  output {
    File updated_gtf = "~{updated_gtf_out}"
    File annotated_id_mappings = "~{updated_mapping_out}"
    File annotated_gene_sparse_tar_gz = "~{gene_sparse_tar_out}"
    File annotated_isoform_sparse_tar_gz = "~{isoform_sparse_tar_out}"
    File annotated_splice_pattern_sparse_tar_gz = "~{splice_sparse_tar_out}"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}
