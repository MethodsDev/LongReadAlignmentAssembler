#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--seurat_obj_rds", help="seurat object in rds format post cell cluster definitions", required=TRUE, nargs=1)
parser$add_argument("--output_matrix", help="name of output matrix file", required=TRUE, nargs=1)

args = parser$parse_args()    
seurat_obj_rds_filename = args$seurat_obj_rds
output_matrix_filename = args$output_matrix

library(Seurat)
library(tidyverse)
library(Matrix)


seurat_obj = readRDS(seurat_obj_rds_filename)

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
