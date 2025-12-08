#!/usr/bin/env Rscript

# CLI: Read a 10x-style sparse matrix, filter to cells with cluster assignments,
# build Seurat object, and run pairwise DE analysis across all clusters.

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(tidyverse)
})

parser <- ArgumentParser(description = "Pairwise DE analysis across clusters from sparse matrix")
parser$add_argument("--sparseM_dir", required = TRUE,
                    help = "Path to sparse matrix directory (10x style: matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz))")
parser$add_argument("--cell_clusters", required = TRUE,
                    help = "Path to cell_clusters.txt file (tab-delimited, no header: cell_barcode<tab>cluster_id)")
parser$add_argument("--output_prefix", required = TRUE,
                    help = "Prefix for output files")

# DE parameters
parser$add_argument("--test_use", default = "wilcox",
                    help = "Test to use for DE (wilcox, MAST, DESeq2, etc.) [default: wilcox]")
parser$add_argument("--min_pct", type = "double", default = 0.1,
                    help = "Minimum fraction of cells expressing the gene in either group [default: 0.1]")
parser$add_argument("--logfc_threshold", type = "double", default = 0.25,
                    help = "Minimum log fold change threshold [default: 0.25]")
parser$add_argument("--only_pos", action = "store_true", default = FALSE,
                    help = "Only return positive markers (genes upregulated in cluster1 vs cluster2)")
parser$add_argument("--pval_adj_cutoff", type = "double", default = 0.05,
                    help = "Adjusted p-value cutoff for filtering final results [default: 0.05]")

# Seurat object building parameters
parser$add_argument("--min_cells", type = "integer", default = 3,
                    help = "Keep genes expressed in at least this many cells [default: 3]")
parser$add_argument("--min_features", type = "integer", default = 200,
                    help = "Keep cells with at least this many detected genes [default: 200]")

args <- parser$parse_args()

sparseM_dir      <- args$sparseM_dir
cell_clusters_file <- args$cell_clusters
output_prefix    <- args$output_prefix
test_use         <- args$test_use
min_pct          <- args$min_pct
logfc_threshold  <- args$logfc_threshold
only_pos         <- args$only_pos
pval_adj_cutoff  <- args$pval_adj_cutoff
min_cells        <- args$min_cells
min_features     <- args$min_features

# -----------------------------------------------
# 1. Load cell cluster assignments
# -----------------------------------------------
message("[1/6] Loading cell cluster assignments from: ", cell_clusters_file)
cell_clusters_df <- read.table(cell_clusters_file, header = FALSE, sep = "\t", 
                               stringsAsFactors = FALSE)

# Assign column names
colnames(cell_clusters_df) <- c("cell_barcode", "cluster_id")

message("- Found ", nrow(cell_clusters_df), " cells with cluster assignments")
message("- Number of clusters: ", length(unique(cell_clusters_df$cluster_id)))

# -----------------------------------------------
# 2. Load sparse matrix
# -----------------------------------------------
message("[2/6] Reading 10x sparse matrix from: ", sparseM_dir)
data <- Read10X(
  data.dir = sparseM_dir,
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE
)

message("- Matrix dimensions: ", nrow(data), " genes x ", ncol(data), " cells")

# -----------------------------------------------
# 3. Filter cells to match cluster assignments
# -----------------------------------------------
message("[3/6] Filtering sparse matrix to cells in cluster assignments")
cells_to_keep <- intersect(colnames(data), cell_clusters_df$cell_barcode)

if (length(cells_to_keep) == 0) {
  stop("Error: No matching cells found between sparse matrix and cell_clusters file")
}

message("- Cells in sparse matrix: ", ncol(data))
message("- Cells in cluster file: ", nrow(cell_clusters_df))
message("- Matching cells: ", length(cells_to_keep))

data_filtered <- data[, cells_to_keep, drop = FALSE]
message("- Filtered matrix dimensions: ", nrow(data_filtered), " genes x ", ncol(data_filtered), " cells")

# -----------------------------------------------
# 4. Create Seurat object
# -----------------------------------------------
message("[4/6] Creating Seurat object")
seurat_obj <- CreateSeuratObject(
  counts = data_filtered,
  project = "project",
  min.cells = min_cells,
  min.features = min_features
)

message("- Seurat object created with ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " genes")

# Add cluster assignments to metadata
cluster_vector <- cell_clusters_df$cluster_id
names(cluster_vector) <- cell_clusters_df$cell_barcode
seurat_obj@meta.data$cluster_id <- cluster_vector[colnames(seurat_obj)]

# Set cluster_id as active identity
Idents(seurat_obj) <- "cluster_id"

message("- Cluster distribution:")
print(table(Idents(seurat_obj)))

# Normalize data (required for DE analysis)
message("- Normalizing data")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                           scale.factor = 10000, verbose = FALSE)

# Save Seurat object
seurat_rds <- paste0(output_prefix, "-seurat_obj.rds")
saveRDS(seurat_obj, file = seurat_rds)
message("- Wrote Seurat object: ", seurat_rds)

# -----------------------------------------------
# 5. Perform pairwise DE analysis
# -----------------------------------------------
message("[5/6] Running pairwise DE analysis across clusters")

perform_pairwise_de <- function(seurat_obj, group_by, 
                               test_use = "wilcox", min_pct = 0.1, 
                               logfc_threshold = 0.25, only_pos = TRUE) {
  
  # Set identity
  Idents(seurat_obj) <- group_by
  clusters <- levels(Idents(seurat_obj))
  
  # Initialize results list
  all_results <- list()
  
  # Perform pairwise comparisons
  for(i in 1:length(clusters)) {
    for(j in 1:length(clusters)) {
      
      if (i == j) { next }
      
      cluster1 <- clusters[i]
      cluster2 <- clusters[j]
      
      cat("Comparing", cluster1, "vs", cluster2, "\n")
      
      de_genes <- FindMarkers(
        seurat_obj,
        ident.1 = cluster1,
        ident.2 = cluster2,
        test.use = test_use,
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        only.pos = only_pos
      )
      
      # Add cluster information
      de_genes$cluster1 <- cluster1
      de_genes$cluster2 <- cluster2
      de_genes$gene <- rownames(de_genes)
      
      # Store results
      comparison_name <- paste0(cluster1, "_vs_", cluster2)
      all_results[[comparison_name]] <- de_genes
    }
  }
  
  return(all_results)
}

# Run the analysis
pairwise_de_results <- perform_pairwise_de(
  seurat_obj, 
  group_by = "cluster_id",
  test_use = test_use,
  min_pct = min_pct,
  logfc_threshold = logfc_threshold,
  only_pos = only_pos
)

# -----------------------------------------------
# 6. Save results
# -----------------------------------------------
message("[6/6] Saving results")

# Combine all results into a single data frame
combined_results <- do.call(rbind, lapply(names(pairwise_de_results), function(x) {
  df <- pairwise_de_results[[x]]
  df$comparison <- x
  return(df)
}))

# Write all results
all_results_file <- paste0(output_prefix, ".cluster_pairwise_DE_results.all.tsv")
write.table(combined_results, file = all_results_file, sep = "\t", 
            quote = FALSE, row.names = FALSE)
message("- Wrote all DE results: ", all_results_file)

# Write filtered results (significant only)
filtered_results <- combined_results %>% filter(p_val_adj <= pval_adj_cutoff)
filtered_results_file <- paste0(output_prefix, ".cluster_pairwise_DE_results.sig.tsv")
write.table(filtered_results, file = filtered_results_file, sep = "\t", 
            quote = FALSE, row.names = FALSE)
message("- Wrote significant DE results (p_val_adj <= ", pval_adj_cutoff, "): ", 
        filtered_results_file)
message("  - ", nrow(filtered_results), " significant DE results")

# Generate summary statistics
DE_df_counts <- combined_results %>% 
  select(cluster1, cluster2, gene) %>% 
  unique() %>%
  group_by(cluster1, cluster2) %>% 
  tally()

summary_file <- paste0(output_prefix, ".cluster_pairwise_DE_summary.tsv")
write.table(DE_df_counts, file = summary_file, sep = "\t", 
            quote = FALSE, row.names = FALSE)
message("- Wrote DE summary: ", summary_file)

# Generate heatmap plot of DE counts
heatmap_pdf <- paste0(output_prefix, ".cluster_pairwise_DE_counts_heatmap.pdf")
p <- ggplot(DE_df_counts, aes(x = factor(cluster1), y = factor(cluster2), fill = n)) + 
  geom_tile() +
  geom_text(aes(label = n), color = "white", size = 3) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  theme_bw() +
  labs(
    x = "Cluster 1 (upregulated)",
    y = "Cluster 2 (downregulated)",
    fill = "# DE genes",
    title = "Pairwise DE gene counts across clusters"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(filename = heatmap_pdf, plot = p, width = 8, height = 7)
message("- Wrote DE counts heatmap: ", heatmap_pdf)

message("Done.")

quit(save = "no", status = 0, runLast = FALSE)
