#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()
parser$add_argument("--sparseM_dir", help="sparse matrix directory (10x style)", required=TRUE, nargs=1)
parser$add_argument("--output_matrix", help="name of output matrix file", required=TRUE, nargs=1)

args = parser$parse_args()    
sparseM_dir = args$sparseM_dir
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

pseudobulk = AggregateExpression(seurat_obj, normalization.method = "RC")

message("-writing ", output_matrix_filename)
write.table(pseudobulk, file="test.tsv", sep="\t", quote=F, col.names=F)

