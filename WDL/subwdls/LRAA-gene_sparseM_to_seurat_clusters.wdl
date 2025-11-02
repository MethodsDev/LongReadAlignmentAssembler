version 1.0

workflow GeneSparseM_To_SeuratClusters {
  input {
    String sample_id
    File gene_sparse_tar_gz
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    Int memoryGB = 32

    # Seurat parameters (defaults aligned with the Rmd)
    Int min_cells = 10
    Int min_features = 1000
    Float percent_mt_max = 20.0
    String mt_pattern = "^MT-"
    Int npcs = 12
    Float resolution = 0.6
    Int n_variable_features = 2000
    Int seed = 1
  }

  String output_prefix = sample_id + ".genes"

  call run_seurat_from_gene_sparseM as run_seurat {
    input:
      gene_sparse_tar_gz = gene_sparse_tar_gz,
      output_prefix = output_prefix,
      docker = docker,
      memoryGB = memoryGB,
      min_cells = min_cells,
      min_features = min_features,
      percent_mt_max = percent_mt_max,
      mt_pattern = mt_pattern,
      npcs = npcs,
      resolution = resolution,
      n_variable_features = n_variable_features,
      seed = seed
  }

  output {
    File seurat_rds = run_seurat.seurat_rds
    File umap_pdf = run_seurat.umap_pdf
    File umap_with_clusters_tsv = run_seurat.umap_with_clusters_tsv
    File cluster_assignments_tsv = run_seurat.cluster_assignments_tsv
  }
}


task run_seurat_from_gene_sparseM {
  input {
    File gene_sparse_tar_gz
    String output_prefix
    String docker
    Int memoryGB = 32

    Int min_cells
    Int min_features
    Float percent_mt_max
    String mt_pattern
    Int npcs
    Float resolution
    Int n_variable_features
    Int seed
  }

  Int disksize = 50 + ceil(2 * size(gene_sparse_tar_gz, "GB"))

  command <<<'
    set -ex

    # Extract gene sparse matrix into a fixed directory name
    mkdir gene-sparseM
    tar -xzf ~{gene_sparse_tar_gz} --strip-components=1 -C gene-sparseM

    # Run the Seurat pipeline script
    Rscript gene_sparseM_to_seurat_clusters_and_umap.R \
      --sparseM_dir gene-sparseM \
      --output_prefix ~{output_prefix} \
      --min_cells ~{min_cells} \
      --min_features ~{min_features} \
      --percent_mt_max ~{percent_mt_max} \
      --mt_pattern '~{mt_pattern}' \
      --npcs ~{npcs} \
      --resolution ~{resolution} \
      --n_variable_features ~{n_variable_features} \
      --seed ~{seed} \
      > seurat_run.log 2>&1 || {
        echo "Seurat script failed; tailing log" >&2
        tail -n 200 seurat_run.log >&2
        exit 1
      }
  '''

  output {
    File seurat_rds = "~{output_prefix}-seurat_obj.rds"
    File umap_pdf = "~{output_prefix}-umap.pdf"
    File umap_with_clusters_tsv = "~{output_prefix}-cell_cluster_assignments.wUMAP.tsv"
    File cluster_assignments_tsv = "~{output_prefix}-cell_cluster_assignments.tsv"
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} HDD"
  }
}
