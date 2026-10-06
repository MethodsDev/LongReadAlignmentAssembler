#!/usr/bin/env Rscript
# Read-track figures for the showcase events of the test run (read_tracks/manifest.tsv):
# read_tracks.<gene>.<kind>.<class>.pdf, drawn with the same functions as the
# site-usage notebook (util/sc/notebook_templates/site_usage_funcs.R).
suppressMessages({library(tidyverse); library(cowplot)})
lraa_home = Sys.getenv("LRAA_HOME", normalizePath("../../.."))
source(file.path(lraa_home, "util/sc/notebook_templates/site_usage_funcs.R"))

# cluster labels, as the notebook's get_cluster_labels() makes them from its UMAP table
cell_types = read_tsv("data/cluster_cell_types.tsv", show_col_types = FALSE, col_types = "cc")
cells = read_tsv("data/cell_clusters.tsv", show_col_types = FALSE, col_types = "cc")
get_cluster_labels = function(max_chars = 34) {
    cells %>% count(seurat_clusters = cluster, name = "n_cells") %>%
        left_join(cell_types, by = c("seurat_clusters" = "cluster")) %>%
        mutate(cluster = paste0("Cluster_", seurat_clusters),
               cluster_label = paste0(str_trunc(cell_type, max_chars), " (", seurat_clusters, "; ", n_cells, " cells)"))
}

m = read_tsv("read_tracks/manifest.tsv", show_col_types = FALSE, col_types = cols(.default = "c"))
options(site_usage.overwrite_figures = TRUE)
for (i in seq_len(nrow(m))) {
    f = paste0("read_tracks.", m$tag[i], ".pdf")
    plot_site_event_read_tracks(m[i, ], "read_tracks", "data/models.gtf", file = f)
    message("wrote ", f)
}
