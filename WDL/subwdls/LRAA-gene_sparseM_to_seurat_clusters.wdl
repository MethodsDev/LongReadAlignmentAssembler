version 1.0

workflow GeneSparseM_To_SeuratClusters {
  input {
    String sample_id
    File gene_sparse_tar_gz
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"
    Int memoryGB = 32

    # Seurat parameters (defaults aligned with the Rmd)
    Int min_cells = 10
    Int min_features = 1000
    Float percent_mt_max = 20.0
    String mt_pattern = "^(MT-|mt-|g:(chrM|MT|M):)"
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
    File seurat_rds_initial = run_seurat.seurat_rds_initial
    File seurat_rds = run_seurat.seurat_rds
    File umap_pdf = run_seurat.umap_pdf
    File umap_with_clusters_tsv = run_seurat.umap_with_clusters_tsv
    File cluster_assignments_tsv = run_seurat.cluster_assignments_tsv
  }
}


task run_seurat_from_gene_sparseM {
  input {
      Int preemptible_tries = 3
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

    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int c3d_cpu = 4
    Int c3d_mem = memoryGB
    Int c3d_cpu_tier = if c3d_cpu <= 4 then 4
        else if c3d_cpu <= 8 then 8
        else if c3d_cpu <= 16 then 16
        else if c3d_cpu <= 30 then 30
        else if c3d_cpu <= 60 then 60
        else if c3d_cpu <= 90 then 90
        else if c3d_cpu <= 180 then 180
        else 360
    Int c3d_mem_tier = if c3d_mem <= 32 then 4
        else if c3d_mem <= 64 then 8
        else if c3d_mem <= 128 then 16
        else if c3d_mem <= 240 then 30
        else if c3d_mem <= 480 then 60
        else if c3d_mem <= 720 then 90
        else if c3d_mem <= 1440 then 180
        else 360
    Int c3d_effective_cpu = if c3d_cpu_tier >= c3d_mem_tier then c3d_cpu_tier else c3d_mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int c3d_highcpu_ram = if c3d_effective_cpu == 4 then 8
        else if c3d_effective_cpu == 8 then 16
        else if c3d_effective_cpu == 16 then 32
        else if c3d_effective_cpu == 30 then 59
        else if c3d_effective_cpu == 60 then 118
        else if c3d_effective_cpu == 90 then 177
        else if c3d_effective_cpu == 180 then 354
        else 708
    String c3d_machine_type = if c3d_mem <= c3d_highcpu_ram
        then "c3d-highcpu-${c3d_effective_cpu}"
        else if c3d_mem <= c3d_effective_cpu * 4
        then "c3d-standard-${c3d_effective_cpu}"
        else "c3d-highmem-${c3d_effective_cpu}"

  command <<<
    set -ex

    # Extract gene sparse matrix into a fixed directory name
    mkdir gene-sparseM
    tar -xzf ~{gene_sparse_tar_gz} --no-same-owner --strip-components=1 -C gene-sparseM

    # Run the Seurat pipeline script (executable on PATH)
    gene_sparseM_to_seurat_clusters_and_umap.R \
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
  >>>

  output {
    File seurat_rds_initial = "~{output_prefix}-seurat_obj.initial.rds"
    File seurat_rds = "~{output_prefix}-seurat_obj.rds"
    File umap_pdf = "~{output_prefix}-umap.pdf"
    File umap_with_clusters_tsv = "~{output_prefix}-cell_cluster_assignments.wUMAP.tsv"
    File cluster_assignments_tsv = "~{output_prefix}-cell_cluster_assignments.tsv"
  }

  runtime {

      preemptible: preemptible_tries
    docker: docker
    predefinedMachineType: c3d_machine_type
    disks: "local-disk ~{disksize} SSD"
  }
}
