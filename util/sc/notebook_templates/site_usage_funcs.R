# Helpers for the site-first TSS / PolyA usage notebooks (util/sc/site_usage/).
#
# Reads the outputs of prep_site_table.py, count_site_read_ends.py, site_usage_dexseq.R,
# annotate_site_usage_events.py and classify_site_pairs_by_splicing.py, and draws site-usage figures in the same style as
# the isoform-level DTU figures: it reuses draw_isoform_share_series and
# draw_isoform_expression_series from iso_plot_funcs.R (source that first), with the
# two sites of an event standing in for the two isoforms of a DTU pair.
#
# Expects in the global environment, as iso_plot_funcs.R's notebooks set them up:
#   umap_df              cell_barcode, umap_1, umap_2, seurat_clusters and a cell-type column
#   cell_neighbor_graph  (optional) from load_cell_neighbor_graph(), for smoothed umaps

library(tidyverse)
library(Matrix)
library(cowplot)

SITE_KINDS = c("TSS", "PolyA")

# classify_site_pairs_by_splicing.py's classes, in reporting order, with display names
SPLICING_CLASSES = c("alt_terminal_usage" = "alternative terminal usage",
                     "alt_splicing:terminal_exon" = "alt splicing: terminal exon",
                     "alt_splicing:internal" = "alt splicing: internal",
                     "unspliced_site" = "unspliced site",
                     "unresolved" = "unresolved")


# All tables of one analysis directory, by site kind.
load_site_usage_results = function(prefix = "PBMCs", dexseq_prefix = paste0(prefix, ".dexseq")) {

    read = function(f) read_tsv(f, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)

    res = list(sites = read(paste0(prefix, ".sites.tsv")),
               count_summary = read(paste0(prefix, ".summary.tsv")),
               ends_at_no_site = read(paste0(prefix, ".ends_at_no_site.tsv")))

    for (kind in SITE_KINDS) {
        p = paste0(dexseq_prefix, ".", kind)
        res[[kind]] = list(genes = read(paste0(p, ".genes.tsv")),
                           site_seeds = read(paste0(p, ".sites.tsv")),
                           usage = read(paste0(p, ".cluster_usage.tsv.gz")),
                           pairwise = read(paste0(p, ".pairwise.tsv.gz")),
                           events = read(paste0(p, ".events.tsv")) %>% mutate(flags = replace_na(flags, "")))
    }

    # each event's site pair classified by the splicing of the reads at the two sites
    split_file = paste0(dexseq_prefix, ".site_pairs.splicing.tsv")
    if (file.exists(split_file)) {
        split = read(split_file)
        for (kind in SITE_KINDS) {
            res[[kind]]$events = res[[kind]]$events %>%
                left_join(split %>% filter(kind == !!kind) %>% select(-kind, -gene_symbol),
                          by = c("gene_key", "gained_site", "lost_site")) %>%
                mutate(splicing_class = factor(SPLICING_CLASSES[splicing_class], levels = SPLICING_CLASSES),
                       splicing_group = case_when(str_starts(splicing_class, "alt splicing") ~ "alternative splicing",
                                                  splicing_class == "alternative terminal usage" ~ "alternative terminal usage",
                                                  TRUE ~ "unclassified"))
        }
    }
    res
}


# sites x cells read-end counts, as written by count_site_read_ends.py
load_site_counts = function(prefix = "PBMCs") {
    m = readMM(paste0(prefix, ".site_counts.mtx.gz"))
    m = as(m, "CsparseMatrix")
    rownames(m) = read_tsv(paste0(prefix, ".sites.tsv.gz"), col_names = FALSE, show_col_types = FALSE)$X1
    colnames(m) = read_tsv(paste0(prefix, ".barcodes.tsv.gz"), col_names = FALSE, show_col_types = FALSE)$X1
    m
}


cluster_cell_counts_from_umap = function() {
    umap_df %>% count(seurat_clusters) %>% { setNames(.$n, paste0("Cluster_", .$seurat_clusters)) }
}


# Events with both clusters of at least min_cluster_cells cells.
filter_events_by_cluster_size = function(events, min_cluster_cells = 200) {
    sizes = cluster_cell_counts_from_umap()
    events %>% filter(sizes[cluster_A] >= min_cluster_cells, sizes[cluster_B] >= min_cluster_cells)
}


# Step-by-step counts for one site kind: from the site collection to high-confidence
# events, as sites, genes and gene x cluster-pair events.
site_usage_funnel = function(res, kind, min_cluster_cells = 200) {

    sites = res$sites %>% filter(kind == !!kind)
    r = res[[kind]]
    ev = r$events
    big = filter_events_by_cluster_size(ev, min_cluster_cells)

    n_genes = function(x) n_distinct(x$gene_key)

    competing = sites %>% filter(competing)
    multi = competing %>% count(gene_key) %>% filter(n >= 2)

    tribble(
        ~step, ~sites, ~genes, ~events,
        "sites in the collection", nrow(sites), NA, NA,
        "  assigned to one gene (competing)", nrow(competing), n_distinct(competing$gene_key), NA,
        "  in genes with >= 2 sites", sum(multi$n), nrow(multi), NA,
        "tested (site filter; >= 2 sites left)", nrow(r$site_seeds), nrow(r$genes), NA,
        "significant gene, seed 1 (q < 0.05)", NA, sum(r$genes$q_seed1 < 0.05, na.rm = TRUE), NA,
        "stable gene (q < 0.05 in >= 4 of 5 seeds)", NA, sum(r$genes$stable), NA,
        "pairwise event (site padj < 0.05, |delta| >= 0.2, >= 20 gene reads per cluster)", NA, n_genes(ev), nrow(ev),
        "  reciprocal", NA, n_genes(filter(ev, switch_class == "reciprocal")), sum(ev$switch_class == "reciprocal"),
        "  reciprocal, FSM >= 5 at both sites", NA, n_genes(filter(ev, switch_class == "reciprocal", both_sites_FSM)),
            sum(ev$switch_class == "reciprocal" & ev$both_sites_FSM),
        "  high confidence (also no flags)", NA, n_genes(filter(ev, high_confidence)), sum(ev$high_confidence),
        paste0("  high confidence, both clusters >= ", min_cluster_cells, " cells"), NA,
            n_genes(filter(big, high_confidence)), sum(big$high_confidence)
    ) %>% mutate(kind = kind, .before = 1)
}


# The alt TSS / alt PolyA funnel, split by splicing: genes tested -> stable -> with
# events, then each gene's events by splicing class, and within each class the
# quality tiers (nested: reciprocal; + FSM at both sites; + no flags = high
# confidence; + both clusters of >= min_cluster_cells cells and >= min_gene_reads
# gene read ends each = showcase-eligible). Gene counts: a gene is counted in every
# class it has an event of, so the class rows overlap; `best` counts each gene once,
# by its best event (best_event_per_gene).
splicing_funnel = function(res, kind, min_cluster_cells = 200, min_gene_reads = 50) {

    r = res[[kind]]
    ev = r$events
    sizes = cluster_cell_counts_from_umap()
    ev = ev %>% mutate(
        tier_reciprocal = switch_class == "reciprocal",
        tier_FSM = tier_reciprocal & both_sites_FSM,
        tier_high_confidence = high_confidence,
        tier_showcase = high_confidence & sizes[cluster_A] >= min_cluster_cells & sizes[cluster_B] >= min_cluster_cells &
            gene_reads_A >= min_gene_reads & gene_reads_B >= min_gene_reads)

    head = tibble(step = c("genes tested (>= 2 sites)", "stable significant genes", "genes with site-switch events"),
                  genes = c(nrow(r$genes), sum(r$genes$stable), n_distinct(ev$gene_key)),
                  events = c(NA, NA, nrow(ev)))

    best = ev %>% best_event_per_gene() %>% count(splicing_class, name = "best")
    by_class = ev %>% group_by(splicing_class) %>%
        summarize(genes = n_distinct(gene_key), events = n(),
                  reciprocal = n_distinct(gene_key[tier_reciprocal]),
                  `+ FSM both sites` = n_distinct(gene_key[tier_FSM]),
                  `+ no flags (high confidence)` = n_distinct(gene_key[tier_high_confidence]),
                  # in practice the same as the FSM tier: flagged sites (monoexonic, downstream TSS without
                  # FSM reads) are the ones lacking FSM support
                  `+ showcase-eligible` = n_distinct(gene_key[tier_showcase]), .groups = "drop") %>%
        left_join(best, by = "splicing_class") %>%
        relocate(best, .after = genes) %>%
        complete(splicing_class, fill = list(genes = 0L, events = 0L, best = 0L)) %>%
        mutate(across(where(is.numeric), ~ replace_na(.x, 0L)))

    list(head = head %>% mutate(kind = kind, .before = 1),
         by_class = by_class %>% mutate(kind = kind, .before = 1))
}


# One row per gene: its best event (high confidence first, then largest |delta|).
best_event_per_gene = function(events) {
    events %>%
        arrange(desc(high_confidence), desc(switch_class == "reciprocal"), desc(abs_delta)) %>%
        group_by(gene_key) %>% slice(1) %>% ungroup() %>%
        arrange(desc(high_confidence), desc(abs_delta))
}


# Names an event's two sites by where they lie along the transcript, e.g.
# "proximal PolyA (chr19:48330553)" / "distal PolyA (...)", "upstream TSS" / "downstream TSS".
site_event_labels = function(event, kind) {
    plus = str_split_fixed(event$gene_key, fixed("|"), 3)[, 3] == "+"
    downstream = if (plus) event$gained_pos > event$lost_pos else event$gained_pos < event$lost_pos
    where = if (kind == "PolyA") c("distal PolyA", "proximal PolyA") else c("downstream TSS", "upstream TSS")
    gained_where = if (downstream) where[1] else where[2]
    lost_where = if (downstream) where[2] else where[1]
    chrom = str_split_fixed(event$gene_key, fixed("|"), 3)[, 2]
    setNames(c(paste0(gained_where, " (", format(event$gained_pos, big.mark = ","), ")"),
               paste0(lost_where, " (", format(event$lost_pos, big.mark = ","), ")")),
             c(event$gained_site, event$lost_site))
}


# Per-cell share of the gene's read ends (over its tested sites) at each of site_ids;
# smoothed when smooth_graph is given by pooling counts over each cell's neighbourhood and
# dividing once (as get_isoform_usage_fraction_umap does for isoforms). Only cells with
# reads of the gene's tested sites of their own are returned (cell_evidence = "gene").
get_site_usage_cells = function(site_counts, site_ids, gene_sites, smooth_graph = NULL) {

    m = site_counts[gene_sites, , drop = FALSE]
    den = Matrix::colSums(m)
    cells = colnames(m)

    if (! is.null(smooth_graph)) {
        common = intersect(cells, rownames(smooth_graph))
        G = smooth_graph[common, common]
        m = m[, common, drop = FALSE]
        den_own = den[common]
        sm_den = as.vector(G %*% den_own)
        out = bind_rows(lapply(site_ids, function(s) {
            tibble(site_id = s, cell_barcode = common, own_reads = den_own,
                   usage = ifelse(sm_den > 0, as.vector(G %*% as.vector(m[s, ])) / sm_den, NA_real_))
        }))
    } else {
        out = bind_rows(lapply(site_ids, function(s) {
            tibble(site_id = s, cell_barcode = cells, own_reads = den,
                   usage = ifelse(den > 0, as.vector(m[s, ]) / den, NA_real_))
        }))
    }

    out %>% filter(own_reads > 0, ! is.na(usage)) %>% inner_join(umap_df, by = "cell_barcode")
}


plot_site_usage_umap = function(event, kind, res, site_counts, smooth_graph = NULL, point_size = 0.6) {

    labels = site_event_labels(event, kind)
    gene_sites = res[[kind]]$site_seeds %>% filter(gene_key == event$gene_key) %>% pull(site_id)
    df = get_site_usage_cells(site_counts, names(labels), gene_sites, smooth_graph) %>%
        mutate(site = factor(labels[site_id], levels = unname(labels))) %>%
        arrange(usage)

    n_cells = n_distinct(df$cell_barcode)
    centers = umap_df %>% group_by(seurat_clusters) %>%
        summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2), .groups = "drop")

    ggplot(umap_df, aes(umap_1, umap_2)) +
        geom_point(color = "gray", alpha = 0.1, size = point_size) +
        geom_point(data = df, aes(color = usage), size = point_size) +
        geom_text(data = centers, aes(label = seurat_clusters), size = 4, color = "purple", fontface = "bold") +
        facet_wrap(~ site, ncol = 2) +
        scale_color_viridis_c(limits = c(0, 1), name = "site usage") +
        labs(title = paste0(event$gene_symbol, ": alternative ", kind, " site usage per cell"),
             subtitle = paste0("share of the gene's read ends at its tested ", kind, " sites (ends at no site excluded)",
                               if (is.null(smooth_graph)) "" else ", SNN-smoothed",
                               " (", format(n_cells, big.mark = ","), " cells with reads of the gene)")) +
        theme_void() +
        theme(plot.margin = margin(t = 18, r = 6, b = 6, l = 18))
}


# cluster x site read counts and usage as the matrices the iso_plot_funcs.R drawing helpers
# take (rows = sites, columns = clusters)
site_cluster_matrix = function(res, kind, value = c("reads", "usage")) {
    value = match.arg(value)
    res[[kind]]$usage %>%
        select(site_id, cluster, all_of(value)) %>%
        pivot_wider(names_from = cluster, values_from = all_of(value), values_fill = 0) %>%
        column_to_rownames("site_id") %>% as.matrix()
}


# the pairwise row for a site between two clusters, in either orientation
find_site_pairwise_row = function(pairwise, site_id, cluster_x, cluster_y) {
    pairwise %>% filter(site_id == !!site_id,
                        (cluster_A == cluster_x & cluster_B == cluster_y) | (cluster_A == cluster_y & cluster_B == cluster_x)) %>%
        head(1)
}


site_event_setup = function(event, kind) {
    labels = site_event_labels(event, kind)
    levels = c(unname(labels), "other isoforms")   # draw_isoform_share_series keys "other" on this name
    list(isoform_ids = names(labels), isoform_levels = levels,
         fill_values = setNames(c(ISOFORM_PAIR_COLORS, OTHER_ISOFORMS_COLOR), levels))
}


relabel_other = function(p, setup) {
    suppressMessages(p + scale_fill_manual(values = setup$fill_values, drop = FALSE, name = NULL,
                                           labels = function(x) sub("other isoforms", "other sites", x)))
}


# The event's two sites as stacked shares of the gene's read ends across a series of
# clusters (default: the event's two), with the shift bands and each adjacent comparison's
# pairwise padj for the gained site and switch class.
plot_site_event_shares = function(event, kind, res, clusters = NULL) {

    if (is.null(clusters)) clusters = c(event$cluster_A, event$cluster_B)
    setup = site_event_setup(event, kind)
    usage = res[[kind]]$usage

    cluster_pi = usage %>% filter(site_id %in% setup$isoform_ids, cluster %in% clusters) %>%
        transmute(cluster, isoform_id = site_id, pi = usage)

    pw = lapply(seq_len(length(clusters) - 1), function(k)
        find_site_pairwise_row(res[[kind]]$pairwise, event$gained_site, clusters[k], clusters[k + 1]))
    band_p = sapply(pw, function(r) if (nrow(r)) r$padj else NA_real_)
    band_class = sapply(seq_len(length(clusters) - 1), function(k) {
        e = res[[kind]]$events %>% filter(gene_key == event$gene_key,
                                          (cluster_A == clusters[k] & cluster_B == clusters[k + 1]) |
                                          (cluster_A == clusters[k + 1] & cluster_B == clusters[k]))
        if (nrow(e)) e$switch_class[1] else NA_character_
    })

    p = draw_isoform_share_series(cluster_pi, clusters, setup, band_p, band_class) +
        labs(y = paste0("share of read ends at the\ngene's tested ", kind, " sites"),
             title = paste0(event$gene_symbol, ": ", kind, " site shares across clusters"))
    relabel_other(p, setup)
}


plot_site_event_expression = function(event, kind, res, clusters = NULL) {
    if (is.null(clusters)) clusters = c(event$cluster_A, event$cluster_B)
    setup = site_event_setup(event, kind)
    draw_isoform_expression_series(clusters, setup, site_cluster_matrix(res, kind, "reads")) +
        labs(y = paste0("read ends per million ", kind, " site read ends"),
             title = paste0(event$gene_symbol, ": read ends at the two sites"))
}


# (A) smoothed per-cell site-usage umaps, (B) site shares across the clusters, (C) read ends
# per million at the two sites in the same clusters.
plot_site_event = function(event, kind, res, site_counts, smooth_graph = NULL, clusters = NULL,
                           file = NULL, width = 10, height = 13.5) {

    p = plot_grid(plot_site_usage_umap(event, kind, res, site_counts, smooth_graph),
                  plot_site_event_shares(event, kind, res, clusters) + theme(plot.margin = margin(t = 18, r = 6, b = 6, l = 18)),
                  plot_site_event_expression(event, kind, res, clusters) + theme(plot.margin = margin(t = 18, r = 6, b = 6, l = 18)),
                  ncol = 1, rel_heights = c(1.15, 1, 0.95), labels = c("A", "B", "C"))

    if (! is.null(file)) {
        ggsave(p, file = file, width = width, height = height,
               device = if (capabilities("cairo")) cairo_pdf else "pdf")
    }
    p
}


# The events table, trimmed to the columns worth reading.
event_table_view = function(events, kind) {
    if (! "isoform_alt_termini_DTU_same_pair" %in% names(events)) {
        events$isoform_alt_termini_DTU_same_pair = NA
    }
    events %>%
        transmute(gene = gene_symbol, clusters = paste0(sub("Cluster_", "", cluster_A), " -> ", sub("Cluster_", "", cluster_B)),
                  type = event_type,
                  splicing = if ("splicing_class" %in% names(events)) as.character(splicing_class) else NA,
                  delta = round(abs_delta, 2),
                  gained = paste0(round(100 * gained_usage_A), "% -> ", round(100 * gained_usage_B), "%"),
                  lost = paste0(round(100 * lost_usage_A), "% -> ", round(100 * lost_usage_B), "%"),
                  sep_bp = separation,
                  FSM_gained = gained_max_isoform_uniq_FSM, FSM_lost = lost_max_isoform_uniq_FSM,
                  switch = switch_class, flags,
                  seeds = n_seeds_significant,
                  iso_DTU_same_pair = isoform_alt_termini_DTU_same_pair)
}
