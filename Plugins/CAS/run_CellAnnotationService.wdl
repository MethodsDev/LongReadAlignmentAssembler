version 1.1

workflow RunCellAnnotationService {
  input {
    File matrix_dir_tarball
    String? matrix_dir_subpath
    String api_token
    String output_prefix_basename
    Int chunk_size = 500
    Float min_acceptable_score = 0.2
    Int top_k = 3
    String obs_prefix = "cas_cell_type"
    String? cas_model_name
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas:latest"
  }

  call RunCAS {
    input:
      matrix_dir_tarball = matrix_dir_tarball,
      matrix_dir_subpath = matrix_dir_subpath,
      api_token = api_token,
      output_prefix_basename = output_prefix_basename,
      chunk_size = chunk_size,
      min_acceptable_score = min_acceptable_score,
      top_k = top_k,
      obs_prefix = obs_prefix,
      cas_model_name = cas_model_name,
      docker_image = docker_image
  }

  output {
    File cas_h5ad = RunCAS.cas_h5ad
    File cas_tsv = RunCAS.cas_tsv
  }

  meta {
    description: "Run Cellarium CAS annotations on a 10x-style sparse matrix produced by LRAA."
  }

  parameter_meta {
    matrix_dir_tarball: {
      description: "Tar.gz archive of the TenX-style sparse matrix directory (create via: tar -czf matrix.tgz matrix_dir)."
    }
    matrix_dir_subpath: {
      description: "Optional relative path inside the tarball that contains matrix.mtx/barcodes/feature files; leave unset to auto-detect."
    }
    api_token: {
      description: "Cellarium CAS API token; store this securely in Terra (e.g., in a workspace secret)."
    }
    output_prefix_basename: {
      description: "Filename prefix for the generated outputs (without extension)."
    }
    chunk_size: {
      description: "Number of cells per CAS chunk." 
    }
    min_acceptable_score: {
      description: "Minimum acceptable evidence score when computing most granular calls."
    }
    top_k: {
      description: "Number of CAS cell type calls to retain per cell."
    }
    obs_prefix: {
      description: "Prefix for the CAS cell type columns in AnnData.obs."
    }
    cas_model_name: {
      description: "Optional CAS model override; leave unset to use Cellarium's default."
    }
    docker_image: {
      description: "Docker image containing run_CellAnnotationService.py and dependencies."
    }
  }
}

task RunCAS {
  input {
    File matrix_dir_tarball
    String? matrix_dir_subpath
    String matrix_dir_subpath_value = select_first([matrix_dir_subpath, ""])
    String api_token
    String output_prefix_basename
    Int chunk_size
    Float min_acceptable_score
    Int top_k
    String obs_prefix
    String? cas_model_name
    String cas_model_name_value = select_first([cas_model_name, ""])
    String cas_model_arg = if cas_model_name_value != "" then "--cas-model-name \"" + cas_model_name_value + "\"" else ""
    String docker_image
  }

  command <<<
    set -euo pipefail

    MATRIX_ROOT=$(mktemp -d cas_matrix.XXXX)
    tar -xzf "~{matrix_dir_tarball}" -C "$MATRIX_ROOT"

    if [[ -n "~{matrix_dir_subpath_value}" ]]; then
      CAS_MATRIX_DIR="$MATRIX_ROOT/~{matrix_dir_subpath_value}"
    else
      subdir_count=$(find "$MATRIX_ROOT" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')
      file_count=$(find "$MATRIX_ROOT" -mindepth 1 -maxdepth 1 -type f | wc -l | tr -d ' ')
      if [[ "$file_count" -eq 0 && "$subdir_count" -eq 1 ]]; then
        CAS_MATRIX_DIR=$(find "$MATRIX_ROOT" -mindepth 1 -maxdepth 1 -type d | head -n 1)
      else
        CAS_MATRIX_DIR="$MATRIX_ROOT"
      fi
    fi

    run_CellAnnotationService.py \
      --matrix-dir "$CAS_MATRIX_DIR" \
      --api-token "~{api_token}" \
      --output-prefix "~{output_prefix_basename}" \
      --chunk-size ~{chunk_size} \
      --min-acceptable-score ~{min_acceptable_score} \
      --top-k ~{top_k} \
      --obs-prefix "~{obs_prefix}" \
      ~{cas_model_arg}
  >>>

  runtime {
    docker: docker_image
    cpu: 4
    memory: "16 GiB"
    disks: "local-disk 200 HDD"
  }

  output {
    File cas_h5ad = output_prefix_basename + ".h5ad"
    File cas_tsv = output_prefix_basename + ".tsv"
  }

  meta {
    description: "Wrapper task that runs run_CellAnnotationService.py inside the CAS Docker image."
  }

  parameter_meta {
    matrix_dir_tarball: {
      description: "Tar.gz archive that contains the sparse matrix directory."
    }
    matrix_dir_subpath: {
      description: "Optional relative path inside the tarball to use as the matrix directory."
    }
    api_token: {
      description: "Cellarium CAS API token (sensitive)."
    }
    output_prefix_basename: {
      description: "Prefix used for the generated .h5ad and .tsv outputs."
    }
    chunk_size: {
      description: "Number of cells sent to CAS per chunk."
    }
    min_acceptable_score: {
      description: "Minimum acceptable evidence score."
    }
    top_k: {
      description: "Top-k CAS cell type calls retained."
    }
    obs_prefix: {
      description: "Prefix for CAS columns written to AnnData.obs."
    }
    cas_model_name: {
      description: "Optional CAS model name override."
    }
    docker_image: {
      description: "Docker image reference for the CAS runner."
    }
  }
}
