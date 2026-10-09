# Read-support views for the SQLIKE gtf analysis notebooks (PBMC and ALS share this file, so
# the two reports partition and draw support the same way). The tables come from
# util/misc/collect_sqlike_read_support.py:
#
#   uniq_reads metric:  uniquely supported / no unique reads
#   FSM_reads metric:   unique FSM read       uniq_FSM_reads >= 1
#                       shared FSM read only  uniq_FSM_reads == 0, has_FSM_read == 1
#                       no FSM read           has_FSM_read == 0
#                       (monoexonic features are n/a: no intron chain to reproduce)
#
# The notebook defines `ordered_cats` (the SQANTI-like category order) before using these.

SQLIKE_SUPPORT_LEVELS = c("uniquely supported", "no unique reads",
                          "unique FSM read", "shared FSM read only",
                          "no FSM read", "FSM unmeasured (NA)",
                          "monoexonic (n/a)", "no quant record")

# supported is opaque, unsupported is faint, and the middle level sits between them so the
# three-way FSM stack reads as a gradient from attested to unobserved
SQLIKE_SUPPORT_ALPHAS = c("uniquely supported"=1, "no unique reads"=0.2,
                          "unique FSM read"=1, "shared FSM read only"=0.5, "no FSM read"=0.15,
                          "FSM unmeasured (NA)"=0.35,
                          "monoexonic (n/a)"=0.35, "no quant record"=0.35)

model_types = c("denovo-basic", "denovo-scg", "refGuided-basic", "refGuided-scg")
spC_types = paste0(model_types, "-spC")

read_sqlike_support = function(summary_tsv, ordered_cats) {
    df = read.csv(summary_tsv, header=T, sep="\t")
    df$Category = factor(df$Category, levels=ordered_cats)
    df$support = factor(df$support, levels=SQLIKE_SUPPORT_LEVELS)
    df
}

# The three views of a class differ only in how each bar is partitioned, and the two classes
# differ only in which categories they keep, so one function emits a class and nothing can drift
# between levels or classes. In the split views, supported sits at the base of the stack and
# unsupported rides on top, so the unsupported slice is read off the top of the bar rather than
# against the axis.

counts_plot = function(df, caption=NULL) {
  df %>%
    ggplot(aes(x=Category, y=Count, fill=type)) + geom_bar(stat='identity', position='dodge') +
    theme_bw() +
    facet_wrap(~Category, nrow=1, scale='free_x') +
    labs(caption=caption) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

support_split_plot = function(df, type_levels, caption=NULL) {
  df %>%
    mutate(type = factor(type, levels=type_levels)) %>%
    ggplot(aes(x=type, y=Count, fill=type, alpha=support)) +
    geom_bar(stat='identity', position=position_stack(reverse=TRUE), color='black', linewidth=0.1) +
    theme_bw() +
    facet_wrap(~Category, nrow=1, scale='free_x') +
    scale_alpha_manual(values=SQLIKE_SUPPORT_ALPHAS, drop=FALSE) +
    labs(caption=caption) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
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

  split_view = function(which_metric) {
    support_df %>%
      filter(level == model_level, metric == which_metric, Category %in% keep_cats) %>%
      support_split_plot(type_levels, caption=caption)
  }

  views = list(counts = counts_df %>% filter(Category %in% keep_cats) %>% counts_plot(caption),
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

# FSM-category features unsupported under each metric (unique reads; unique FSM read)
fsm_unsupported_table = function(support_df) {
  support_df %>%
    filter(Category == "FSM", ! support %in% c("monoexonic (n/a)", "no quant record")) %>%
    group_by(type, level, metric) %>%
    summarize(features = sum(Count),
              supported = sum(Count[support %in% c("uniquely supported", "unique FSM read")]),
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
