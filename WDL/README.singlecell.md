# LRAA single-cell end-to-end pipeline

This workflow orchestrates a full single-cell analysis using LRAA:

1) Initial transcript discovery and quantification from a single BAM (LRAA.wdl)
2) Build single-cell gene/isoform/splice-pattern sparse matrices from the initial tracking
3) Cluster cells using Seurat from the gene-level sparse matrix
4) Cluster-guided LRAA to refine transcripts and produce final single-cell outputs

Workflow file: `WDL/LRAA-singlecell.wdl`

## Inputs

Required:
- `sample_id` – label for outputs
- `referenceGenome` – FASTA file (indexed) for the reference genome
- `inputBAM` – single-cell long-read BAM with CB/UMI tags (default CB/XM)

Optional:
- `initial_annot_gtf` – GTF to guide the initial LRAA discovery (not required)
- `HiFi` – set true for HiFi reads
- `oversimplify` – comma-separated contigs (e.g., "chrM" or "chrM,MT") to oversimplify
- `main_chromosomes` – space-separated contigs to shard; empty runs direct mode
- `region` – restricts run to a specific region (forces direct mode)

Resources:
- `numThreadsPerWorker` and `numThreadsPerWorkerScattered` control LRAA worker threads for direct and chromosome-sharded runs, respectively
- `memoryGB` is an optional override for direct LRAA runs; when unset, `LRAA.wdl` computes memory as `max(64 GiB, ceil(1.5 x input BAM size in GiB))`
- `memoryGBPerWorkerScattered` is an optional override for chromosome-sharded LRAA workers only; when unset, each shard self-sizes in `subwdls/LRAA_runner.wdl` with a `32 GiB` floor and larger allocations for mid-size and large shard BAMs
- `main_chromosomes` determines whether LRAA runs direct or chromosome-sharded; if empty, `memoryGBPerWorkerScattered` has no effect because the run stays in direct mode
- In cluster-guided single-cell mode, the outer per-cluster scatter does not by itself activate `memoryGBPerWorkerScattered`; that setting only applies if each per-cluster `LRAA.wdl` call is also chromosome-sharded via `main_chromosomes`
- `diskSizeGB` defaults to `256`
- Sparse matrix build tuning: `memoryGBbuildSparseMatrices`, `memoryGBscSparseMatrices`, `sparseMatrixCsvEngine` (`python` by default), `sparseMatrixGzipLevel` (`1` by default for faster compression)

Seurat clustering parameters (defaults mirror included R pipeline):
- `min_cells` (10), `min_features` (1000), `percent_mt_max` (20.0), `mt_pattern` ("^MT-")
- `npcs` (12), `resolution` (0.6), `n_variable_features` (2000), `seed` (1)

See `WDL/example_inputs/LRAA-singlecell.inputs.json` for a starting template.

## Outputs

From initial discovery:
- `init_quant_expr` – initial pseudobulk quantification
- `init_quant_tracking` – initial read-to-isoform tracking (.gz)
- `init_gtf` – initial reconstructed transcripts GTF
- `init_sc_*_sparse_tar_gz` – tarballs of gene/isoform/splice-pattern sparse matrices
- `seurat_umap_pdf`, `seurat_cluster_assignments` – clustering diagnostics

Final cluster-guided results (primary deliverables):
- `final_gtf` – final transcript GTF (cluster-guided discovery + merge)
- `final_tracking` – merged per-cluster tracking (.gz)
- `final_sc_*_sparse_tar_gz` – final single-cell sparse matrices (gene/iso/splice)
- Cluster-level pseudobulk matrices and tarballs of per-cluster quant outputs

## How to run

With Cromwell:

```bash
java -jar cromwell.jar run WDL/LRAA-singlecell.wdl \
  --inputs WDL/example_inputs/LRAA-singlecell.inputs.json
```

With miniwdl (local backend):

```bash
miniwdl run WDL/LRAA-singlecell.wdl \
  sample_id=PBMC10k \
  referenceGenome=/path/genome.fa \
  inputBAM=/path/reads.bam
```

Notes:
- Ensure the reference FASTA is indexed (.fai) and the BAM is indexed (.bai). LRAA will attempt to index BAMs as needed.
- The container image includes R/Seurat and LRAA utilities invoked by the subworkflows.
- Sparse matrix construction now uses the single-pass Python implementation; the old `--parallel` mode is no longer used by the WDL tasks.
- If you already have a cell-cluster assignment TSV, you can run `WDL/LRAA-cell_cluster_guided.wdl` directly by providing `cell_clusters_info` and an annotation GTF.
