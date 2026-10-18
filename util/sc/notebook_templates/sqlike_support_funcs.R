# Read-support views for the SQLIKE gtf analysis notebooks (PBMC and ALS share this file, so
# the two reports partition and draw support the same way). Two metrics:
#
#   uniq_reads metric (quant-based, util/misc/collect_sqlike_read_support.py):
#                       uniquely supported / no unique reads
#   FSM_reads metric (same-reads, SQLIKE read mode of the original aligned BAM against each GTF;
#                     the paper repo's __util/read_fsm_support_summary.py):
#                       FSM read in a clustered cell         >= 1 read reproduces the intron chain,
#                                                            from a barcode in a cell cluster
#                       FSM read, unclustered barcodes only  such reads, none from a clustered barcode
#                       no FSM read                          no read reproduces the chain
#                       (monoexonic features are n/a: no intron chain to reproduce)
#
# The collector also writes a quant-based FSM_reads metric (unique / shared FSM read, from the
# quant.expr). read_sqlike_support() replaces it with the same-reads one when given that table:
# scg is quantified per cell cluster, so its quant never sees reads from unclustered barcodes and
# quant-based FSM support is not comparable between basic and scg runs.
#
# The notebook defines `ordered_cats` (the SQANTI-like category order) before using these.

SQLIKE_SUPPORT_LEVELS = c("uniquely supported", "no unique reads",
                          "FSM read in a clustered cell", "FSM read, unclustered barcodes only",
                          "unique FSM read", "shared FSM read only",
                          "no FSM read", "FSM unmeasured (NA)",
                          "monoexonic (n/a)", "no quant record")

# supported is opaque, unsupported is faint, and the middle level sits between them so the
# three-way FSM stack reads as a gradient from attested to unobserved
SQLIKE_SUPPORT_ALPHAS = c("uniquely supported"=1, "no unique reads"=0.2,
                          "FSM read in a clustered cell"=1, "FSM read, unclustered barcodes only"=0.5,
                          "unique FSM read"=1, "shared FSM read only"=0.5, "no FSM read"=0.15,
                          "FSM unmeasured (NA)"=0.35,
                          "monoexonic (n/a)"=0.35, "no quant record"=0.35)

model_types = c("denovo-basic", "denovo-scg", "refGuided-basic", "refGuided-scg")
spC_types = paste0(model_types, "-spC")

read_sqlike_support = function(summary_tsv, ordered_cats, read_mode_fsm_tsv=NULL) {
    df = read.csv(summary_tsv, header=T, sep="\t")
    if (! is.null(read_mode_fsm_tsv)) {
        df = bind_rows(df %>% filter(metric != "FSM_reads"),
                       read.csv(read_mode_fsm_tsv, header=T, sep="\t"))
    }
    df$Category = factor(df$Category, levels=ordered_cats)
    df$support = factor(df$support, levels=SQLIKE_SUPPORT_LEVELS)
    df
}

# The legend of each view lists only the levels its metric can take (the alpha mapping itself is
# fixed, so a level is shaded the same in every view)
SQLIKE_METRIC_LEVELS = list(
  uniq_reads = c("uniquely supported", "no unique reads"),
  FSM_reads  = c("FSM read in a clustered cell", "FSM read, unclustered barcodes only",
                 "unique FSM read", "shared FSM read only", "no FSM read", "monoexonic (n/a)"))

SQLIKE_METRIC_TITLES = c(uniq_reads = "unique-read support", FSM_reads = "FSM-read support")

SQLIKE_LEVEL_UNITS = c(model = "transcript models", splice_chain = "splice-collapsed chains")

# The three views of a class differ only in how each bar is partitioned, and the two classes
# differ only in which categories they keep, so one function emits a class and nothing can drift
# between levels or classes. In the split views, supported sits at the base of the stack and
# unsupported rides on top, so the unsupported slice is read off the top of the bar rather than
# against the axis. Bars are the four LRAA runs, in the same order and colour in every view, and
# the facet strip names the category. Within a panel the runs sit as two abutting pairs, de novo
# (basic, scg) and ref-guided (basic, scg), with a gap between the pairs, so basic vs scg is read
# within a pair.

SQLIKE_BAR_WIDTH = 0.95
SQLIKE_PAIR_GAP = 0.6

# x position of each run: consecutive within a pair, offset by the gap between pairs (the pair is
# the run name up to its first '-', i.e. denovo or refGuided)
run_x = function(type, type_levels) {
  pair = sub("-.*", "", type_levels)
  pos = seq_along(type_levels) + SQLIKE_PAIR_GAP * (match(pair, unique(pair)) - 1)
  pos[match(type, type_levels)]
}

sqlike_bar_theme = function() {
  list(theme_bw(),
       theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
             panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
             strip.text = element_text(size=8)))
}

counts_plot = function(df, type_levels, y_label, caption=NULL) {
  df %>%
    mutate(x = run_x(type, type_levels), type = factor(type, levels=type_levels)) %>%
    ggplot(aes(x=x, y=Count, fill=type)) +
    geom_bar(stat='identity', width=SQLIKE_BAR_WIDTH, color='black', linewidth=0.1) +
    facet_wrap(~Category, nrow=1, scale='free_x') +
    sqlike_bar_theme() +
    labs(x=NULL, y=y_label, fill="LRAA run", caption=caption)
}

support_split_plot = function(df, type_levels, metric, y_label, caption=NULL) {
  shown_levels = intersect(SQLIKE_METRIC_LEVELS[[metric]], as.character(unique(df$support)))
  df %>%
    mutate(x = run_x(type, type_levels), type = factor(type, levels=type_levels)) %>%
    ggplot(aes(x=x, y=Count, fill=type, alpha=support)) +
    geom_bar(stat='identity', width=SQLIKE_BAR_WIDTH, position=position_stack(reverse=TRUE),
             color='black', linewidth=0.1) +
    facet_wrap(~Category, nrow=1, scale='free_x') +
    scale_alpha_manual(values=SQLIKE_SUPPORT_ALPHAS, breaks=shown_levels, drop=FALSE) +
    guides(fill=guide_legend(order=1),
           alpha=guide_legend(order=2, override.aes=list(fill="grey20", color="black"))) +
    sqlike_bar_theme() +
    labs(x=NULL, y=y_label, fill="LRAA run", alpha=SQLIKE_METRIC_TITLES[[metric]],
         caption=caption)
}

# In the all-categories class the FSM view keeps the monoexonic features, shaded as n/a: they
# have no intron chain to reproduce, so they are unmeasurable here rather than unsupported.
# The spliced-only class drops them, which is the same thing said by omission.
make_three_views = function(counts_df, support_df, model_level, type_levels, file_prefix,
                            spliced_only=FALSE) {
  spliced_cats = grep("^se_", ordered_cats, value=TRUE, invert=TRUE)
  keep_cats = if (spliced_only) spliced_cats else ordered_cats
  caption = if (spliced_only) "spliced categories only" else NULL
  width = if (spliced_only) 11 else 15
  prefix = paste0(file_prefix, if (spliced_only) "_spliced" else "")
  y_label = SQLIKE_LEVEL_UNITS[[model_level]]

  split_view = function(which_metric) {
    support_df %>%
      filter(level == model_level, metric == which_metric, Category %in% keep_cats) %>%
      support_split_plot(type_levels, which_metric, y_label, caption=caption)
  }

  views = list(counts = counts_df %>% filter(Category %in% keep_cats) %>%
                 counts_plot(type_levels, y_label, caption),
               uniq   = split_view("uniq_reads"),
               FSM    = split_view("FSM_reads"))

  ggsave(views$counts, file=paste0(prefix, "_plot.pdf"), width=width, height=6)
  ggsave(views$uniq, file=paste0(prefix, "_uniq_read_support_plot.pdf"), width=width, height=6)
  ggsave(views$FSM, file=paste0(prefix, "_FSM_read_support_plot.pdf"), width=width, height=6)

  views
}

# each metric must re-partition exactly the features the counts view plots
check_support_partitions = function(support_df, counts_df, model_level) {
  stopifnot(
    support_df %>%
      filter(level == model_level) %>%
      group_by(type, metric, Category = as.character(Category)) %>%
      summarize(Count = sum(Count), .groups="drop") %>%
      inner_join(counts_df %>%
                   transmute(type, Category = as.character(Category), total = Count),
                 by=c("type", "Category")) %>%
      summarize(all_equal = all(Count == total)) %>%
      pull(all_equal)
  )
}

# FSM-category features unsupported under each metric (unique reads; an FSM read, or for the
# quant-based FSM metric a unique FSM read)
fsm_unsupported_table = function(support_df) {
  support_df %>%
    filter(Category == "FSM", ! support %in% c("monoexonic (n/a)", "no quant record")) %>%
    group_by(type, level, metric) %>%
    summarize(features = sum(Count),
              supported = sum(Count[support %in% c("uniquely supported", "unique FSM read",
                                                    "FSM read in a clustered cell",
                                                    "FSM read, unclustered barcodes only")]),
              .groups="drop") %>%
    mutate(pct_unsupported = (features - supported) / features * 100) %>%
    arrange(level, type, metric) %>%
    as.data.frame()
}

# which explanation the FSM-category features without a unique FSM read have: the chain was
# never traversed by a whole read (no FSM read), or its whole reads are shared with another
# model (shared FSM read only)
fsm_never_traversed_table = function(support_df) {
  support_df %>%
    filter(metric == "FSM_reads", Category == "FSM",
           support %in% c("unique FSM read", "shared FSM read only", "no FSM read")) %>%
    select(type, level, support, Count) %>%
    pivot_wider(names_from=support, values_from=Count, values_fill=0) %>%
    mutate(pct_never_traversed = `no FSM read` /
             (`unique FSM read` + `shared FSM read only` + `no FSM read`) * 100) %>%
    arrange(level, type) %>%
    as.data.frame()
}

# where the FSM reads of the FSM-category features come from (same-reads FSM metric): a barcode in
# a cell cluster, only barcodes outside the clusters, or no read reproduces the chain at all
fsm_read_origin_table = function(support_df) {
  support_df %>%
    filter(metric == "FSM_reads", Category == "FSM",
           support %in% c("FSM read in a clustered cell", "FSM read, unclustered barcodes only",
                          "no FSM read")) %>%
    select(type, level, support, Count) %>%
    pivot_wider(names_from=support, values_from=Count, values_fill=0) %>%
    mutate(features = `FSM read in a clustered cell` + `FSM read, unclustered barcodes only` +
             `no FSM read`,
           pct_unclustered_only = `FSM read, unclustered barcodes only` / features * 100,
           pct_no_FSM_read = `no FSM read` / features * 100) %>%
    arrange(level, type) %>%
    as.data.frame()
}
