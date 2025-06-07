#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--sparseM_dir", help="sparse matrix directory (10x style)", required=TRUE, nargs=1)
parser$add_argument("--clusters", help="file containing cell_barcode(tab)cluster_id (ignores header row)", required=TRUE, nargs=1)
parser$add_argument("--output_matrix", help="name of output matrix file", required=TRUE, nargs=1)

args = parser$parse_args()    
sparseM_dir = args$sparseM_dir
clusters_file = args$clusters
output_matrix_filename = args$output_matrix

library(Seurat)
library(tidyverse)
library(Matrix)


message("Reading ", sparseM_dir)
data = Read10X(data.dir=sparseM_dir,
               gene.column = 1,
               cell.column = 2,
               unique.features = TRUE,
               strip.suffix = FALSE)


message("-creating pseudobulk matrix")

seurat_obj = seurat_obj <- CreateSeuratObject(counts = data)

# Read in the cluster assignments (cell_barcode and cluster_id)
cluster_df <- read.delim(clusters_file, header = FALSE, skip=1, stringsAsFactors = FALSE)
colnames(cluster_df) = c("cell_barcode", "cluster_id")

# Ensure the cell barcodes match those in the Seurat object
# Filter to matching barcodes just in case
common_barcodes <- intersect(cluster_df$cell_barcode, colnames(seurat_obj))

cluster_df <- cluster_df[cluster_df$cell_barcode %in% common_barcodes, ]

# Set cell barcodes as rownames (required for setting identities)
rownames(cluster_df) <- cluster_df$cell_barcode

# Reorder to match Seurat object's column order
cluster_df <- cluster_df[colnames(seurat_obj), , drop = FALSE]

# Add cluster assignments as metadata
seurat_obj$assigned_cluster <- as.factor(cluster_df$cluster_id)

# Set the identities (Idents) of the Seurat object to these cluster labels
Idents(seurat_obj) <- seurat_obj$assigned_cluster


###############################
## get pseudobulks per cluster

# Extract the raw counts matrix
counts <- GetAssayData(seurat_obj, slot = "counts")

# Get cluster assignments
clusters <- Idents(seurat_obj)


    
# Use Matrix package to efficiently sum counts per gene across cells within each cluster
pseudobulk_matrix <- t(sapply(levels(clusters), function(clust) {
  # Get the column indices for cells in this cluster
  cells_in_cluster <- WhichCells(seurat_obj, idents = clust)
  
  # Subset the count matrix and sum across columns (cells)
  Matrix::rowSums(counts[, cells_in_cluster, drop = FALSE])
}))

# Transpose so that rows are genes and columns are clusters
pseudobulk_matrix <- t(pseudobulk_matrix)

# Optional: add gene names as rownames
rownames(pseudobulk_matrix) <- rownames(counts)
colnames(pseudobulk_matrix) <- paste0("Cluster_", levels(clusters))

message("writing ", output_matrix_filename)
write.table(pseudobulk_matrix, output_matrix_filename, sep="\t", quote=F)

quit(save = "no", status = 0, runLast = FALSE)


