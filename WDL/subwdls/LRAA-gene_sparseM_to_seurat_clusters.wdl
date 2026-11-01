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
    # Machine shape (v0.47.1, generated text, identical in every C3D task). cpu: and memory:
    # stay as they are, so miniwdl's reservations do not change; this only names the Terra VM
    # (miniwdl ignores predefinedMachineType). The command's thread count stays on the
    # requested cpu even when the shape has more cores. Batch rejects 8 GiB on c3d-highcpu-4
    # (7.5 GiB passes), hence 7 for that tier: such a request steps up to c3d-standard-4.
    Int shape_cpu = 4
    Int shape_memGB = memoryGB
    Int shape_cpu_tier = if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 30 then 30
        else if shape_cpu <= 60 then 60
        else if shape_cpu <= 90 then 90
        else 180
    Int shape_mem_tier = if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 240 then 30
        else if shape_memGB <= 480 then 60
        else if shape_memGB <= 720 then 90
        else 180
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    Int shape_highcpu_memGB = if shape_cores == 4 then 7
        else if shape_cores == 30 then 59
        else if shape_cores == 60 then 118
        else if shape_cores == 90 then 177
        else if shape_cores == 180 then 354
        else shape_cores * 2
    String machine_type = if shape_memGB <= shape_highcpu_memGB then "c3d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "c3d-standard-~{shape_cores}"
        else "c3d-highmem-~{shape_cores}"


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

      predefinedMachineType: machine_type

      preemptible: preemptible_tries
    docker: docker
    cpu: 4
    memory: "~{memoryGB} GiB"
    disks: "local-disk ~{disksize} SSD"
  }
}
