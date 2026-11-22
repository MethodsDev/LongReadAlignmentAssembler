# LRAA single-cell end-to-end pipeline

This workflow orchestrates a full single-cell analysis using LRAA:

1) Initial transcript discovery and quantification from a single BAM (LRAA.wdl)
2) Build single-cell gene/isoform/splice-pattern sparse matrices from the initial tracking
3) Cluster cells using Seurat from the gene-level sparse matrix
4) Cluster-guided LRAA to refine transcripts (discovery per cluster, then a cross-cluster merge)
5) One pooled quantification of all cells against the final transcripts, from which the final
   single-cell matrices and the per-cluster pseudobulk outputs are built

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
- `cpu` and `cpuScattered` are the core counts for direct and chromosome-sharded LRAA tasks. Each is the task's `cpu` runtime request and the `--cpu_budget` LRAA divides across its work units, so the cores requested and the cores LRAA spends are the same number. They replace `numThreadsPerWorker`, `numThreadsPerWorkerScattered` and `num_parallel_contigs`, which were multiplied rather than divided: `numThreadsPerWorker 5` alongside `num_parallel_contigs 3` let LRAA start three workers of five threads on a five-core task
- `memoryGB` overrides `memoryGB_whole_genome` (`32` GiB), the box a whole-genome (non-scattered) LRAA run gets. Fixed rather than derived from the input BAM: every run chunks, so peak RSS follows chunk concurrency and chunk width, not library size
- `memoryGBPerWorkerScattered` overrides `memoryGB_per_chromosome_shard` (`16` GiB), the box every chromosome shard gets. Uniform across shards on purpose -- a chunked shard's peak follows its chunk concurrency, which `max_cpu_per_chunked_shard` caps, not the length of its contig
- There is no `chunk` input: every LRAA run in these workflows chunks. That is what makes the two numbers above fixed rather than functions of input size, and it is why a workspace still binding `chunk` will fail on an unknown input
- `main_chromosomes` determines whether LRAA runs direct or chromosome-sharded; if empty, `memoryGBPerWorkerScattered` has no effect because the run stays in direct mode
- In cluster-guided single-cell mode, the outer per-cluster discovery scatter does not by itself activate `memoryGBPerWorkerScattered`; that setting only applies if each per-cluster `LRAA.wdl` call is also chromosome-sharded via `main_chromosomes`
- The cluster-guided final quantification is one whole-library LRAA run, sized like the initial pass: `cpuInit` sets its cpu, `memoryGBquantFinal` its memory when `scattering_final_quant` is `off` (unset: computed from the BAM), and `initial_partition_workers` its chromosome partition; `num_total_reads` is its RPM denominator as for the initial pass
- `diskSizeGB` defaults to `256`
- Sparse matrix build tuning: `memoryGBbuildSparseMatrices`, `memoryGBscSparseMatrices`, `sparseMatrixCsvEngine` (`direct` by default: no CSV parser, 4.1x faster than the python engine; `python` and `c` select the pandas readers), `emit_dense_counts` (off by default; the `*_cell_counts.tsv` files are 111 GiB on a 1.5 B-read library and no workflow output references them), `sparseMatrixGzipLevel` (`1` by default for faster compression)

Seurat clustering parameters (defaults mirror included R pipeline):
- `min_cells` (10), `min_features` (1000), `percent_mt_max` (20.0), `mt_pattern` (`^(MT-|mt-|g:(chrM|MT|M):)`; matches conventional mitochondrial gene symbols and annotation-free LRAA gene IDs on `chrM`, `MT`, or `M`)
- `npcs` (12), `resolution` (0.6), `n_variable_features` (2000), `seed` (1)

Single-nucleus libraries:
- `sc_nucseq` (false) – opt-in single-nucleus mode for snRNA-seq, where most reads are unspliced pre-mRNA. Never switched on from the data. When true:
  - the initial pass and the cluster-guided final quant also run LRAA's gene-body counting (config key `gene_body_counts`, see `docs/single_nucleus_intronic_gene_counts.md`): every read is credited to one gene, intronic and partly spliced reads included, and classed spliced (S), unspliced (U) or ambiguous (A). This holds in both discovery modes, de novo and ref-guided, where reads are credited to the genes of the models the run reports, and in `quant_only` mode, where they are credited to `initial_annot_gtf`'s genes. Isoform outputs are unchanged
  - `gene_full` (S + U + A) matrices are built next to the standard gene matrix, which stays as it is
  - empty-droplet filtering calls cells on `gene_full`, and Seurat clusters on it, with `min_features_nucseq` (1000) in place of `min_features` and variable genes ranked within bins of log10 gene span (gene_full counts grow with gene length). Everything downstream of the cluster assignments is unchanged
  - the key is merged into the config inside the LRAA task, so `init_config_update_json` stays unset and the `quant_only` reuse of the initial pass as the final quant still applies
  - inputs should be strand-specific (oriented 10x / Kinnex reads): an unspliced read can only take its aligned strand. LRAA logs a warning per contig-strand when more than 20% of spliced reads align against their transcribed strand, and the run total is written to `init_gene_body_strand_check` / `final_gene_body_strand_check`
  - the initial pass must run when clustering runs: a `precomputed_init_quant_tracking` carries no gene-body counts
- `min_features_nucseq` (1000) – Seurat's `min.features` under `sc_nucseq`; gene_full detects many more genes per nucleus than the standard matrix
- `min_features_nucseq_min_passing_cells` (200), `min_features_nucseq_min_passing_frac` (0.5), `min_features_nucseq_adaptive_quantile` (0.10), `min_features_nucseq_lower_bound` (100) – with `sc_nucseq` and `enable_filter_good_cells`, `min_features_nucseq` is lowered for a shallow library: when fewer nuclei than `..._min_passing_cells`, or less than `..._min_passing_frac` of them, reach it, Seurat uses their `..._adaptive_quantile` of detected genes instead, never below `..._lower_bound`. The value used and why are in `seurat_min_features`. An unfiltered matrix (all barcodes) keeps the fixed value

See `WDL/example_inputs/LRAA-singlecell.inputs.json` for a starting template.

## Outputs

From initial discovery:
- `init_quant_expr` – initial pseudobulk quantification
- `init_quant_tracking` – initial read-to-isoform tracking (.gz)
- `init_gtf` – initial reconstructed transcripts GTF
- `init_sc_*_sparse_tar_gz` – tarballs of gene/isoform/splice-pattern sparse matrices
- `seurat_umap_pdf`, `seurat_cluster_assignments` – clustering diagnostics

With `sc_nucseq = true` (all null otherwise):
- `init_gene_body_tracking` – the initial pass's gene-body sidecar (`<prefix>.quant.gene_body.tracking.gz`), one row per read and gene with its splice class
- `init_sc_gene_full_sparse_tar_gz` – gene x barcode counts including intronic reads; the matrix cell filtering and clustering used
- `init_sc_gene_spliced_sparse_tar_gz`, `init_sc_gene_unspliced_sparse_tar_gz`, `init_sc_gene_ambiguous_sparse_tar_gz` – its S, U and A parts; all four share one feature and one barcode list, so S + U + A = gene_full entry for entry
- `init_gene_body_cell_stats` – per barcode: `gene_full`, `spliced`, `unspliced`, `ambiguous`, `n_genes_full`, `spliced_frac`, `unspliced_frac` (the intronic fraction; high in nuclei)
- `filtered_gene_full_sparse_tar_gz` – gene_full over the cells emptyDrops called on it (`filtered_gene_sparse_tar_gz` is then the standard gene matrix over the same cells)
- `seurat_gene_spans_tsv`, `seurat_variable_features_tsv` – the gene spans used, and each gene's span bin, within-bin score and whether it was selected
- `seurat_min_features` – (with filtered nuclei) `min_features` requested and used, nuclei passing each, the settings, and the reason
- `init_gene_body_strand_check` – one row: units (contig-strands), `spliced_reads_checked`, `spliced_reads_antisense`, `antisense_frac`, the thresholds, and `verdict` (`stranded`, `unstranded_or_reverse`, or `too_few_spliced_reads`)
- `final_gene_body_tracking`, `final_sc_gene_{full,spliced,unspliced,ambiguous}_sparse_tar_gz`, `final_gene_body_cell_stats`, `final_gene_body_strand_check` – the same from the pooled final quantification (cluster-guided mode)
- `incl_gene_symbols_gene_{full,spliced,unspliced,ambiguous}_sparse_tar_gz` – the unfiltered final (else initial) gene-body matrices with gene symbols, when `ref_annot_gtf_source_gene_symbols` is given

Final cluster-guided results (primary deliverables):
- `final_gtf` – final transcript GTF (cluster-guided discovery + merge)
- `final_tracking` – read tracking of the pooled final quantification (.gz); all cells, including barcodes in no cluster
- `final_quant_expr` – quant.expr of the pooled final quantification
- `final_sc_*_sparse_tar_gz` – final single-cell sparse matrices (gene/iso/splice), from the pooled quantification
- `final_cluster_exprs_tar` – one quant.expr per cluster (`<sample_id>.<cluster>.LRAA.quant-only.quant.expr`), the pooled quantification summed over the cluster's cells; `RPM_total_reads` there is per million of the cluster's assigned reads, so it equals TPM
- `cluster_*_matrix` – cluster pseudobulk count/TPM matrices built from those files
- `final_cluster_read_partition` – assigned reads per cluster, and of barcodes in no cluster (counted in the single-cell matrices, in no cluster's pseudobulk)
- There is no per-cluster tracking output: a read's cluster is the cluster of the cell barcode leading its read name in `final_tracking` (`<barcode>^<UMI>^<read>`)
- With `quant_only = true` the initial pass already is the pooled quantification of `initial_annot_gtf`, so the final outputs reuse it (same tracking, quant.expr and matrices) rather than running a second identical pass. This applies when `scattering_init` equals `scattering_final_quant` and neither `init_config_update_json` nor `region` is set

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
