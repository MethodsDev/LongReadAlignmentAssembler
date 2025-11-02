#!/usr/bin/env Rscript

# CLI: Read a 10x-style gene sparse matrix, run standard Seurat filtering, clustering, UMAP,
# and write outputs: Seurat object (RDS), UMAP plot (PDF), UMAP coordinates + cluster ids (TSV),
# and simple cell->cluster assignments (TSV).

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(tidyverse)
})

parser <- ArgumentParser(description = "Seurat clustering and UMAP from gene sparse matrix (10x format)")
parser$add_argument("--sparseM_dir", required = TRUE,
                    help = "Path to gene-based sparse matrix directory (10x style: matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz))")
parser$add_argument("--output_prefix", required = TRUE,
                    help = "Prefix for output files (eg: sample.genes)")

# Basic filters and params (aligned with the Rmd defaults)
parser$add_argument("--min_cells", type = "integer", default = 10,
                    help = "Keep genes expressed in at least this many cells [default: 10]")
parser$add_argument("--min_features", type = "integer", default = 1000,
                    help = "Keep cells with at least this many detected genes [default: 1000]")
parser$add_argument("--percent_mt_max", type = "double", default = 20.0,
                    help = "Max percent.mt allowed; cells above are filtered [default: 20]")
parser$add_argument("--mt_pattern", default = "^MT-",
                    help = "Regex for mitochondrial genes (for percent.mt) [default: ^MT-]")

# Dimensionality and clustering
parser$add_argument("--npcs", type = "integer", default = 12,
                    help = "Number of PCs to use for neighbors/UMAP [default: 12]")
parser$add_argument("--resolution", type = "double", default = 0.6,
                    help = "Clustering resolution [default: 0.6]")
parser$add_argument("--n_variable_features", type = "integer", default = 2000,
                    help = "Number of variable features for HVG selection [default: 2000]")
parser$add_argument("--seed", type = "integer", default = 1,
                    help = "Random seed for UMAP/clustering reproducibility [default: 1]")

args <- parser$parse_args()

sparseM_dir   <- args$sparseM_dir
output_prefix <- args$output_prefix
min_cells     <- args$min_cells
min_features  <- args$min_features
percent_mt_max <- args$percent_mt_max
mt_pattern     <- args$mt_pattern
npcs          <- args$npcs
resolution    <- args$resolution
n_variable_features <- args$n_variable_features
seed          <- args$seed

set.seed(seed)

message("[1/7] Reading 10x gene sparse matrix from: ", sparseM_dir)
data <- Read10X(
  data.dir = sparseM_dir,
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

message("[2/7] Creating Seurat object and computing QC metrics")
seurat_obj <- CreateSeuratObject(
  counts = data,
  project = "project",
  min.cells = min_cells,
  min.features = min_features
)

# Compute percent.mt and filter cells
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)

message("- Cells before mt-filtering: ", ncol(seurat_obj))
seurat_obj <- subset(seurat_obj, subset = percent.mt < percent_mt_max)
message("- Cells after  mt-filtering: ", ncol(seurat_obj))

message("[3/7] Normalizing and selecting variable features (", n_variable_features, ")")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = n_variable_features, verbose = FALSE)

message("[4/7] Scaling data and regressing out nFeature_RNA, nCount_RNA, percent.mt")
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"), verbose = FALSE)

message("[5/7] PCA, neighbors, clustering")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:npcs, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)

message("[6/7] UMAP embedding (dims = 1:", npcs, ")")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:npcs, verbose = FALSE)

message("[7/7] Writing outputs")
# 1) Seurat object RDS
seurat_rds <- paste0(output_prefix, "-seurat_obj.rds")
saveRDS(seurat_obj, file = seurat_rds)
message("- Wrote Seurat object: ", seurat_rds)

# 2) UMAP plot (PDF)
umap_pdf <- paste0(output_prefix, "-umap.pdf")
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggplot2::theme_bw()
ggplot2::ggsave(filename = umap_pdf, plot = p, width = 7, height = 6)
message("- Wrote UMAP plot: ", umap_pdf)

# 3) UMAP coordinates + cluster id (TSV)
cell_clustering_info <- seurat_obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_barcode") %>%
  dplyr::select(cell_barcode, seurat_clusters)

umap_df <- as.data.frame(seurat_obj[["umap"]]@cell.embeddings) %>%
  tibble::rownames_to_column(var = "cell_barcode") %>%
  as_tibble() %>%
  left_join(cell_clustering_info, by = "cell_barcode")

umap_tsv <- paste0(output_prefix, "-cell_cluster_assignments.wUMAP.tsv")
readr::write_tsv(umap_df, file = umap_tsv)
message("- Wrote UMAP+clusters TSV: ", umap_tsv)

# 4) Simple cell->cluster assignments (TSV)
cluster_assignments <- tibble(
  cell_barcode = names(Idents(seurat_obj)),
  cluster = as.character(Idents(seurat_obj))
)
clusters_tsv <- paste0(output_prefix, "-cell_cluster_assignments.tsv")
readr::write_tsv(cluster_assignments, file = clusters_tsv)
message("- Wrote cluster assignments TSV: ", clusters_tsv)

message("Done.")

quit(save = "no", status = 0, runLast = FALSE)
