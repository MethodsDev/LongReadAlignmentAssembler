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


# Write a figure to `file` (pdf), unless it already exists: re-knitting a notebook
# then leaves existing figure files untouched (cairo_pdf output differs byte-wise from
# run to run even when the figure does not). options(site_usage.overwrite_figures = TRUE)
# rewrites them.
save_figure = function(p, file, width, height) {
    if (is.null(file)) return(invisible(FALSE))
    if (file.exists(file) && ! isTRUE(getOption("site_usage.overwrite_figures", FALSE))) return(invisible(FALSE))
    ggsave(p, file = file, width = width, height = height,
           device = if (capabilities("cairo")) cairo_pdf else "pdf")
    invisible(TRUE)
}

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
        labs(title = paste0(event$gene_symbol, ": ", site_switch_description(kind, event$splicing_class),
                            " -- site usage per cell"),
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
# What kind of switch an event is, for titles: from its read-based splicing class
# (classify_site_pairs_by_splicing.py) or its showcase group ("terminal_usage" / "alt_splicing").
site_switch_description = function(kind, class_or_group) {
    x = as.character(class_or_group)
    if (length(x) == 0 || is.na(x)) x = ""
    terminal = x %in% c("alt_terminal_usage", "terminal_usage", "alternative terminal usage")
    splicing = startsWith(x, "alt_splicing") || startsWith(x, "alt splicing")
    if (terminal) return(if (kind == "TSS") "tandem TSSs" else "tandem 3' UTR PolyA sites")
    if (splicing) return(paste0("alternative ", kind, " with alternative splicing"))
    paste0("alternative ", kind)
}

plot_site_event = function(event, kind, res, site_counts, smooth_graph = NULL, clusters = NULL,
                           file = NULL, width = 10, height = 13.5) {

    p = plot_grid(plot_site_usage_umap(event, kind, res, site_counts, smooth_graph),
                  plot_site_event_shares(event, kind, res, clusters) + theme(plot.margin = margin(t = 18, r = 6, b = 6, l = 18)),
                  plot_site_event_expression(event, kind, res, clusters) + theme(plot.margin = margin(t = 18, r = 6, b = 6, l = 18)),
                  ncol = 1, rel_heights = c(1.15, 1, 0.95), labels = c("A", "B", "C"))

    save_figure(p, file, width, height)
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


# ---------------------------------------------------------------------------------
# Read tracks: isoform structures over a sample of the reads supporting them, per
# cell cluster (reads from util/sc/site_usage/extract_isoform_read_tracks.py).
# ---------------------------------------------------------------------------------

# exons of the given transcripts (ids as in the gtf, e.g. "SELENOH^t:chr11:+:comp-569:iso-10")
read_transcript_exons = function(gtf, transcript_ids) {
    pattern = paste(sprintf('transcript_id "%s"', transcript_ids), collapse = "|")
    lines = system2("grep", c("-E", shQuote(paste0("\t(exon)\t.*(", gsub("([\\^.+])", "\\\\\\1", pattern), ")")), gtf),
                    stdout = TRUE)
    tibble(raw = lines) %>%
        separate(raw, into = c("chrom", "src", "feature", "start", "end", "score", "strand", "frame", "attr"),
                 sep = "\t") %>%
        transmute(transcript_id = str_match(attr, 'transcript_id "([^"]+)"')[, 2],
                  start = as.integer(start), end = as.integer(end), strand)
}


# One read-track panel. `transcripts`: named vector, names = labels, values = gtf
# transcript ids, in drawing order; `track_ids` maps them to the tracking-file ids used
# in `reads` (default: the part after "^"). `clusters`: named vector, names = labels,
# values = cluster numbers, in drawing order. `sites`: named positions to mark (e.g. the
# two TSSs). `xlim`: region shown.
#
# `read_totals` (optional): cluster, transcript_id (tracking id), n -- all reads of each
# isoform in each cluster, shown in the cluster's header so the sample can be read
# against them (count_isoform_reads_by_cluster()).
#
# If `reads` has a `read_class` column, reads other than "uniq_FSM" (e.g. "compatible":
# partial reads sharing the isoform's terminus) are drawn lighter.
plot_read_track_panel = function(exons, reads, transcripts, clusters, sites = NULL, xlim = NULL,
                                 colors = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"),
                                 read_height = 0.7, title = NULL, show_legend = TRUE, read_totals = NULL,
                                 highlight = NULL, base_size = 10, model_scale = 1, site_linewidth = 0.3) {

    tx_label = setNames(names(transcripts), transcripts)
    tx_col = setNames(colors[seq_along(transcripts)], names(transcripts))
    track_id = sub("^[^^]*\\^", "", transcripts)
    reads = reads %>% mutate(isoform = names(transcripts)[match(transcript_id, track_id)],
                             cluster_label = names(clusters)[match(as.character(cluster), as.character(clusters))]) %>%
        filter(! is.na(isoform), ! is.na(cluster_label))

    # rows: isoform models on top, then each cluster's reads, grouped by isoform and 5' start
    plus = exons$strand[1] == "+"
    read_order = reads %>% distinct(cluster_label, read_name, isoform, read_start, read_end) %>%
        mutate(cluster_label = factor(cluster_label, levels = names(clusters)),
               isoform = factor(isoform, levels = names(transcripts)),
               five_prime = if (plus) read_start else -read_end) %>%
        arrange(cluster_label, isoform, five_prime)
    gap = 3
    # model_scale: thickness (and row spacing) of the isoform models relative to a read row
    model_rows = tibble(isoform = names(transcripts), y = -(seq_along(transcripts)) * model_scale)
    # below the bottom model's lower edge (models are model_scale thick), so thick models
    # don't run into the first cluster's header
    y = min(model_rows$y) - 0.5 * model_scale - gap
    row_y = numeric(nrow(read_order))
    header = list()
    for (cl in names(clusters)) {
        idx = which(read_order$cluster_label == cl)
        header[[cl]] = y
        row_y[idx] = y - seq_along(idx)
        y = y - length(idx) - gap
    }
    read_order$y = row_y
    # the y range from all rows, fixed before any trimming to a zoom window: a zoom that
    # loses the top model row or the bottom reads must not rescale, or its rows drift from
    # the other panels' rows
    ylim = c(min(c(row_y, model_rows$y)) - 0.6, max(model_rows$y) + 0.5 * model_scale + 0.1)
    header = tibble(cluster_label = names(header), y = unlist(header))
    if (! is.null(read_totals)) {
        tot = read_totals %>% mutate(isoform = names(transcripts)[match(transcript_id, track_id)],
                                     cluster_label = names(clusters)[match(as.character(cluster), as.character(clusters))]) %>%
            filter(! is.na(isoform), ! is.na(cluster_label)) %>%
            mutate(isoform = factor(isoform, levels = names(transcripts))) %>% arrange(isoform) %>%
            group_by(cluster_label) %>%
            summarize(totals = paste0(sub(" .*", "", isoform), " ", format(n, big.mark = ","), collapse = " / "), .groups = "drop")
        header = header %>% left_join(tot, by = "cluster_label") %>%
            mutate(cluster_label = if_else(is.na(totals), cluster_label, paste0(cluster_label, "   (reads: ", totals, ")")))
    }

    blocks = reads %>% left_join(read_order %>% select(read_name, y), by = "read_name") %>%
        mutate(fsm = if ("read_class" %in% names(.)) read_class == "uniq_FSM" else TRUE)
    spans = read_order
    model_ex = exons %>% mutate(isoform = tx_label[transcript_id]) %>% left_join(model_rows, by = "isoform")
    model_span = model_ex %>% group_by(isoform, y) %>% summarize(start = min(start), end = max(end), .groups = "drop")
    if (! is.null(xlim)) {
        # the panel is unclipped (so the cluster headers can run past it): trim the reads
        # and models to the window here instead, or a zoom's exons spill into the next panel
        trim = function(d, s, e) {
            d %>% filter(.data[[e]] >= xlim[1], .data[[s]] <= xlim[2]) %>%
                mutate("{s}" := pmax(.data[[s]], xlim[1]), "{e}" := pmin(.data[[e]], xlim[2]))
        }
        blocks = trim(blocks, "block_start", "block_end")
        spans = trim(spans, "read_start", "read_end")
        model_ex = trim(model_ex, "start", "end")
        model_span = trim(model_span, "start", "end")
    }

    if (! is.null(highlight) && ! is.null(xlim)) {   # unclipped panel: keep the bands inside it
        highlight = lapply(highlight, function(w) c(max(w[1], xlim[1]), min(w[2], xlim[2])))
        highlight = highlight[vapply(highlight, function(w) w[1] < w[2], logical(1))]
        if (! length(highlight)) highlight = NULL
    }
    p = ggplot() + highlight_layer(highlight) +
        geom_segment(data = spans, aes(x = read_start, xend = read_end, y = y, yend = y), color = "#b8b7b1", linewidth = 0.25) +
        geom_rect(data = blocks, aes(xmin = block_start, xmax = block_end, ymin = y - read_height / 2, ymax = y + read_height / 2,
                                     fill = isoform, alpha = fsm), color = NA) +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
        geom_segment(data = model_span, aes(x = start, xend = end, y = y, yend = y), color = "#0b0b0b",
                     linewidth = 0.4 * sqrt(model_scale)) +
        geom_rect(data = model_ex, aes(xmin = start, xmax = end, ymin = y - 0.38 * model_scale, ymax = y + 0.38 * model_scale,
                                       fill = isoform),
                  color = "#0b0b0b", linewidth = 0.2) +
        geom_text(data = header, aes(x = -Inf, y = y - 0.2, label = cluster_label), hjust = -0.02, vjust = 0,
                  size = base_size * 0.3, color = "#0b0b0b") +
        scale_fill_manual(values = tx_col, breaks = names(transcripts), name = NULL) +
        scale_y_continuous(breaks = NULL) +
        labs(x = NULL, y = NULL, title = title) +
        theme_minimal(base_size = base_size) +
        theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
              legend.position = if (show_legend) "bottom" else "none",
              legend.key.size = unit(9, "pt"), legend.spacing.x = unit(4, "pt"),
              legend.text = element_text(margin = margin(r = 14)))
    if (! is.null(sites)) {
        # the panel is unclipped: a site outside the window would be drawn on the next panel
        in_view = if (is.null(xlim)) sites else sites[sites >= xlim[1] & sites <= xlim[2]]
        if (length(in_view)) {
            p = p + geom_vline(xintercept = in_view, linetype = "dashed", color = "#5a5954", linewidth = site_linewidth)
        }
    }
    p = p + coord_cartesian(xlim = xlim, ylim = ylim, expand = is.null(xlim), clip = if (is.null(xlim)) "on" else "off")
    p + scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                           n.breaks = if (! is.null(xlim) && diff(xlim) < 1000) 3 else 5)
}


# light band behind the data marking the windows (list of c(start, end)) shown zoomed
highlight_layer = function(highlight, fill = "#e4e3dc") {
    if (is.null(highlight)) return(NULL)
    annotate("rect", xmin = sapply(highlight, `[`, 1), xmax = sapply(highlight, `[`, 2),
             ymin = -Inf, ymax = Inf, fill = fill)
}


# Read-end density per cluster: each cluster's read 5' (TSS) or 3' (PolyA) ends in
# `binwidth`-bp bins, as a share of that cluster's ends in the region, one row per
# cluster. `ends`: cluster, pos, reads (extract_isoform_read_tracks.py --ends_output).
plot_read_end_density = function(ends, clusters, sites = NULL, xlim, binwidth = 5, kind = "TSS",
                                 fill = "#5a5954", highlight = NULL, axis_side = "right", ends_at_sites = FALSE,
                                 base_size = 9, show_subtitle = TRUE, site_linewidth = 0.3) {
    d = ends %>% mutate(cluster_label = names(clusters)[match(as.character(cluster), as.character(clusters))]) %>%
        filter(! is.na(cluster_label)) %>%
        group_by(cluster_label) %>% mutate(share = reads / sum(reads)) %>% ungroup() %>%
        filter(pos >= xlim[1], pos <= xlim[2]) %>%
        mutate(bin = floor((pos - xlim[1]) / binwidth) * binwidth + xlim[1] + binwidth / 2) %>%
        group_by(cluster_label, bin) %>% summarize(share = sum(share), .groups = "drop") %>%
        mutate(cluster_label = factor(cluster_label, levels = names(clusters)))
    labels = tibble(cluster_label = factor(names(clusters), levels = names(clusters)))
    # cluster names in the corner away from the sites, so they don't sit on the peaks
    label_right = ! is.null(sites) && mean(sites) < mean(xlim)
    p = ggplot(d, aes(x = bin, y = share)) + highlight_layer(highlight) +
        geom_col(width = binwidth * 0.9, fill = fill) +
        geom_text(data = labels, aes(x = if (label_right) Inf else -Inf, y = Inf, label = cluster_label), inherit.aes = FALSE,
                  hjust = if (label_right) 1.03 else -0.03, vjust = 1.3, size = base_size * 0.3, color = "#0b0b0b") +
        facet_grid(cluster_label ~ .) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1), n.breaks = 3,
                           position = axis_side, expand = expansion(mult = c(0, 0.35))) +
        labs(x = NULL, y = NULL,
             subtitle = if (show_subtitle) paste0("read ", if (kind == "TSS") "5'" else "3'", " ends",
                               if (ends_at_sites) paste0(" at the gene's ", kind, " sites") else "",
                               " per ", binwidth, " bp") else paste0(binwidth, " bp bins")) +
        theme_minimal(base_size = base_size) +
        theme(panel.grid.minor = element_blank(), strip.text = element_blank(), panel.spacing.y = unit(8, "pt"),
              axis.text.x = element_blank(), plot.subtitle = element_text(size = base_size - 0.5))
    if (! is.null(sites)) {
        p = p + geom_vline(xintercept = sites, linetype = "dashed", color = "#5a5954", linewidth = site_linewidth)
    }
    p + coord_cartesian(xlim = xlim, expand = FALSE)
}


# Full-gene view beside a zoom on the varying terminus (TSS: 5' end; PolyA: 3' end),
# optionally with the read-end density of each cluster above the reads (`ends`). Sites
# up to `max_joint_zoom` apart share one zoom; farther apart (alternative first / last
# exons) each site gets its own. The whole-gene density bins scale with the region.
plot_isoform_read_tracks = function(exons, reads, transcripts, clusters, sites, kind = c("TSS", "PolyA"),
                                    zoom_flank = 120, title = NULL, file = NULL, width = 11, height = 8,
                                    read_totals = NULL, ends = NULL, zoom_bin = 4,
                                    density_height = 0.28, max_joint_zoom = 400, ends_at_sites = FALSE,
                                    base_size = 10, zoom = TRUE, model_scale = 1, site_linewidth = 0.3) {
    kind = match.arg(kind)
    full_xlim = range(c(exons$start, exons$end, reads$read_start, reads$read_end))
    # pad by 2% so a read-end peak at the gene's terminus isn't drawn on the panel edge
    full_xlim = full_xlim + c(-1, 1) * max(20, round(0.02 * diff(full_xlim)))
    full_bin = max(10, round(diff(full_xlim) / 250))
    zooms = if (! zoom) {
        list()   # whole gene only
    } else if (diff(range(sites)) <= max_joint_zoom) {
        list(c(min(sites) - zoom_flank, max(sites) + zoom_flank))
    } else {
        lapply(sort(sites), function(p) c(p - zoom_flank, p + zoom_flank))
    }

    xlims = c(list(full_xlim), zooms)
    bins = c(full_bin, rep(zoom_bin, length(zooms)))
    # no zooms: one column, no column title
    col_titles = c(if (length(zooms)) "whole gene" else "", if (length(zooms) == 1) paste0(kind, " region") else if (length(zooms) > 1)
        paste0(kind, " at ", format(sort(sites), big.mark = ",")))
    # the whole-gene column shades the zoomed windows and keeps its density axis on its
    # outer (left) side, with some space before the zooms, so the column boundary is clear
    gap = function(i) if (i > 1) theme(plot.margin = margin(5.5, 5.5, 5.5, 16)) else NULL
    tracks = lapply(seq_along(xlims), function(i)
        plot_read_track_panel(exons, reads, transcripts, clusters, sites = sites, xlim = xlims[[i]],
                              title = if (is.null(ends) && nzchar(col_titles[i])) col_titles[i] else NULL, show_legend = FALSE,
                              read_totals = if (i == 1) read_totals else NULL,
                              highlight = if (i == 1 && length(zooms)) zooms else NULL, base_size = base_size,
                              model_scale = model_scale, site_linewidth = site_linewidth) + gap(i))
    # one grid, so every column's read rows share the same height and line up across
    # panels; the legend goes under the whole figure, not under one column
    widths = c(1.6, rep(if (length(zooms) == 1) 1 else 0.7, length(zooms)))
    p = if (is.null(ends)) {
        plot_grid(plotlist = tracks, nrow = 1, rel_widths = widths, align = "h", axis = "tb")
    } else {
        dens = lapply(seq_along(xlims), function(i)
            plot_read_end_density(ends, clusters, sites, xlims[[i]], binwidth = bins[i], kind = kind,
                                  ends_at_sites = ends_at_sites, base_size = base_size - 1, show_subtitle = i == 1,
                                  highlight = if (i == 1 && length(zooms)) zooms else NULL,
                                  axis_side = if (i == 1) "left" else "right", site_linewidth = site_linewidth) +
                labs(title = if (nzchar(col_titles[i])) col_titles[i] else NULL) + gap(i))
        # patchwork aligns the panels across the grid with each row's axis space sized to
        # that row (cowplot's align = "hv", axis = "tblr" gave the density row the height of
        # the read tracks' x-axis labels as empty space)
        patchwork::patchworkGrob(patchwork::wrap_plots(c(dens, tracks), nrow = 2, widths = widths,
                                                       heights = c(density_height, 1)))
    }
    legend = get_plot_component(tracks[[1]] + theme(legend.position = "bottom"), "guide-box-bottom")
    p = plot_grid(p, legend, ncol = 1, rel_heights = c(1, 0.03))
    if (! is.null(title)) {
        p = plot_grid(ggdraw() + draw_label(title, x = 0.01, hjust = 0, size = base_size + 2), p, ncol = 1, rel_heights = c(0.04, 1))
    }
    save_figure(p, file, if (length(zooms) > 1) width * 1.15 else width, height)
    p
}


# reads per isoform and cluster from an LRAA quant.tracking file (optionally only unique
# full-splice-match reads), for plot_isoform_read_tracks(read_totals = ...)
count_isoform_reads_by_cluster = function(tracking, cell_clusters, unique_FSM_only = TRUE) {
    tr = read_tsv(tracking, comment = "#", show_col_types = FALSE,
                  col_types = cols(.default = "c")) %>% filter(transcript_id != "transcript_id")
    if (unique_FSM_only) tr = tr %>% filter(is_unique == "1", is_FSM == "1")
    # cell_barcode <tab> cluster; any header line (with or without a tab) is dropped
    f = str_split_fixed(readLines(cell_clusters), "\t", 3)
    cl = tibble(cell_barcode = f[, 1], cluster = str_trim(f[, 2])) %>% filter(str_detect(cluster, "^-?[0-9]+$"))
    tr %>% transmute(transcript_id = sub(".*@", "", transcript_id), cell_barcode = sub("\\^.*", "", read_name)) %>%
        inner_join(cl, by = "cell_barcode") %>% count(cluster, transcript_id)
}


# Path of a file in the LRAA checkout: under $LRAA_HOME if set; else found through this
# helper's own symlink in the working directory (notebooks link site_usage_funcs.R from
# LRAA's util/sc/notebook_templates/), which holds whichever user runs R -- RStudio may
# run as a different user, whose ~ is not the checkout owner's; else ~/GITHUB/MDL/.
lraa_util_path = function(rel) {
    home = Sys.getenv("LRAA_HOME")
    if (home == "") {
        link = Sys.readlink("site_usage_funcs.R")
        home = if (! is.na(link) && nzchar(link)) normalizePath(file.path(dirname(link), "../../.."), mustWork = FALSE)
               else path.expand("~/GITHUB/MDL/LongReadAlignmentAssembler")
    }
    file.path(home, rel)
}


# Read-track data for a set of site-switch events (util/sc/site_usage/
# build_site_event_read_tracks.py), cached in `outdir`: rerun only when the events or the
# script change. `events`: tag, gene_symbol, kind, gained_site, lost_site, cluster_A,
# cluster_B. Returns the manifest (one row per event).
build_site_event_read_tracks = function(events, outdir, sites, gtf, cluster_quant_tar, tracking, bam,
                                        cell_clusters, max_reads = 30,
                                        script = lraa_util_path("util/sc/site_usage/build_site_event_read_tracks.py")) {
    if (! file.exists(script)) stop("cannot find ", script, "; set LRAA_HOME to the LRAA checkout")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    events_tsv = file.path(outdir, "events.tsv")
    new_tsv = file.path(outdir, "events.tsv.new")
    write_tsv(events %>% select(tag, gene_symbol, kind, gained_site, lost_site, cluster_A, cluster_B), new_tsv)
    key = paste(unname(tools::md5sum(new_tsv)), unname(tools::md5sum(script)), max_reads)
    key_file = file.path(outdir, "cache_key")
    manifest = file.path(outdir, "manifest.tsv")
    if (! (file.exists(manifest) && file.exists(key_file) && readLines(key_file, n = 1) == key)) {
        file.rename(new_tsv, events_tsv)
        status = system2("python3", c(script, "--events", events_tsv, "--sites", sites, "--gtf", gtf,
                                      "--cluster_quant_tar", cluster_quant_tar, "--tracking", tracking,
                                      "--bam", bam, "--cell_clusters", cell_clusters,
                                      "--max_reads", max_reads, "--outdir", outdir),
                         stdout = file.path(outdir, "build.log"), stderr = file.path(outdir, "build.log"))
        if (status != 0) stop("build_site_event_read_tracks.py failed; see ", file.path(outdir, "build.log"))
        writeLines(key, key_file)
    } else {
        unlink(new_tsv)
    }
    read_tsv(manifest, show_col_types = FALSE, col_types = cols(.default = "c"))
}


# One event's read-track figure from a manifest row: the gained-site isoform (blue) and
# the lost-site isoform (orange) over a sample of each cluster's reads -- unique FSM
# reads first, then (lighter) partial reads compatible with the isoform that share its
# terminus (build_site_event_read_tracks.py) -- with each cluster's read-end density
# above; labelled as in plot_site_event.
plot_site_event_read_tracks = function(m, outdir, gtf, file = NULL, max_chars = 34, ...) {
    kind = m$kind
    gained_pos = as.numeric(m$gained_pos)
    lost_pos = as.numeric(m$lost_pos)
    ev = tibble(gene_key = paste(m$gene_symbol, m$chrom, m$strand, sep = "|"),
                gained_site = m$gained_site, lost_site = m$lost_site,
                gained_pos = gained_pos, lost_pos = lost_pos)
    site_lab = unname(site_event_labels(ev, kind))
    iso = function(x) sub(".*:", "", x)
    transcripts = setNames(c(m$gained_gtf_id, m$lost_gtf_id),
                           c(paste0(site_lab[1], ", ", iso(m$gained_tx)), paste0(site_lab[2], ", ", iso(m$lost_tx))))
    cl_labels = get_cluster_labels(max_chars = max_chars)
    cl = c(m$cluster_A, m$cluster_B)
    clusters = setNames(cl, cl_labels$cluster_label[match(paste0("Cluster_", cl), cl_labels$cluster)])
    exons = read_transcript_exons(gtf, unname(transcripts))
    reads = read_tsv(file.path(outdir, paste0(m$tag, ".reads.tsv")), show_col_types = FALSE)
    ends = read_tsv(file.path(outdir, paste0(m$tag, ".ends.tsv")), show_col_types = FALSE)
    totals = read_tsv(file.path(outdir, paste0(m$tag, ".totals.tsv")), show_col_types = FALSE)
    plot_isoform_read_tracks(exons, reads, transcripts, clusters, sites = c(gained_pos, lost_pos), kind = kind,
                             read_totals = totals, ends = ends, ends_at_sites = TRUE,
                             title = paste0(m$gene_symbol, ": ", site_switch_description(kind, sub(".*\\.", "", m$tag)), ", ",
                                            if (any(reads$read_class %in% "compatible"))
                                                paste0("unique FSM reads (solid) and partial reads sharing the ", kind, " (light), ")
                                            else "unique FSM reads, ",
                                            max(table(reads$cluster[!duplicated(reads$read_name)])), " sampled per cluster"),
                             file = file, height = 10, ...)
}
