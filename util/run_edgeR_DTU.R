#!/usr/bin/env Rscript

# Differential transcript usage (DTU) analysis with edgeR + diffSplice.
#
# Inputs:
#   1. Transcript-level count matrix with required columns:
#        gene_id, transcript_id, <sample1>, <sample2>, ...
#   2. samples.txt with two columns and no header:
#        <condition> <sample>
#   3. contrasts.txt with two columns and no header:
#        <numerator_condition> <denominator_condition>
#
# Statistical workflow:
#   - Subset the matrix to the requested samples.
#   - Precompute transcript fractions within each gene from the full
#     sample-subset matrix BEFORE transcript filtering. These fractions are
#     used only for reporting columns in the transcript-level output.
#   - Build an edgeR DGEList from transcript counts.
#   - Filter low-information transcripts using:
#       * edgeR::filterByExpr()
#       * a user-tunable CPM prevalence rule
#       * retention of genes with at least 2 surviving transcripts
#   - Normalize library sizes with edgeR's normLibSizes().
#   - Fit transcript-level quasi-likelihood GLMs.
#   - Run edgeR/limma diffSplice() grouped by the configured gene grouping key
#     (full gene_id by default; gene_symbol via --gene-grouping) so testing
#     asks whether transcript usage within a gene differs across conditions.
#
# Output logic:
#   - gene_F.tsv: primary gene-level DTU table from the quasi-F test.
#   - gene_simes.tsv: alternate gene-level ranking from Simes aggregation.
#   - transcript_t.tsv: transcript-level table including:
#       * diffSplice transcript statistics
#       * mean raw counts and normalized expression summaries by contrast group
#       * mean isoform fractions and delta_isoform_fraction computed from the
#         pre-filter full-gene denominators.
#
# Important interpretation note:
#   The DTU hypothesis test is performed on the FILTERED transcript set.
#   The reported isoform_fraction columns are computed from the PRE-FILTER
#   sample-subset matrix so they reflect each transcript's fraction of the
#   full observed gene output, not its fraction within only the tested subset.

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "run_edgeR_DTU.R --matrix counts.tsv[.gz] --samples samples.txt --contrasts contrasts.txt --outdir outdir",
  "[--min-cpm 1] [--min-samples 2] [--min-gene-transcripts 2] [--fdr 1] [--gene-grouping gene_id|gene_symbol] [--verbose]",
  sep = "\n"
)

parse_args <- function(args) {
  opts <- list(
    `--matrix` = NULL,
    `--samples` = NULL,
    `--contrasts` = NULL,
    `--outdir` = NULL,
    `--min-cpm` = 1,
    `--min-samples` = 2,
    `--min-gene-transcripts` = 2,
    `--fdr` = 1,
    `--gene-grouping` = "gene_id",
    `--verbose` = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% names(opts)) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    if (identical(key, "--verbose")) {
      opts[[key]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    opts[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  for (required in c("--matrix", "--samples", "--contrasts", "--outdir")) {
    if (is.null(opts[[required]])) {
      stop("Missing required argument: ", required, "\n", usage, call. = FALSE)
    }
  }

  opts[["--min-cpm"]] <- as.numeric(opts[["--min-cpm"]])
  opts[["--min-samples"]] <- as.integer(opts[["--min-samples"]])
  opts[["--min-gene-transcripts"]] <- as.integer(opts[["--min-gene-transcripts"]])
  opts[["--fdr"]] <- as.numeric(opts[["--fdr"]])

  valid_groupings <- c("gene_id", "gene_symbol")
  if (!opts[["--gene-grouping"]] %in% valid_groupings) {
    stop(
      "--gene-grouping must be one of: ", paste(valid_groupings, collapse = ", "),
      " (got '", opts[["--gene-grouping"]], "')",
      call. = FALSE
    )
  }

  opts
}

opts <- parse_args(args)

log_msg <- function(...) {
  if (isTRUE(opts[["--verbose"]])) {
    message(...)
  }
}

read_no_header_tsv <- function(path, colnames) {
  df <- read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    comment.char = "",
    quote = ""
  )
  if (ncol(df) != length(colnames)) {
    stop("Expected ", length(colnames), " columns in ", path, " but found ", ncol(df), call. = FALSE)
  }
  names(df) <- colnames
  df
}

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

make_compact_feature_label <- function(transcript_id, max_chars = 15) {
  first_component <- strsplit(transcript_id, "[\\^|:;,_-]+")[[1]][1]
  if (is.na(first_component) || !nzchar(first_component)) {
    first_component <- transcript_id
  }
  substr(first_component, 1, max_chars)
}

make_transcript_volcano_plot <- function(df, contrast_name, output_file, top_n_each = 25) {
  plot_df <- df
  plot_df$neg_log10_pvalue <- -log10(pmax(plot_df$P.Value, .Machine$double.xmin))
  plot_df$direction <- ifelse(plot_df$logFC > 0, "up", ifelse(plot_df$logFC < 0, "down", "flat"))

  # Rank labels using a combined effect-size / significance priority score.
  # This favors features that are both strongly shifted and well supported.
  plot_df$priority_score <- abs(plot_df$logFC) * plot_df$neg_log10_pvalue
  plot_df$label <- vapply(plot_df$transcript_id, make_compact_feature_label, character(1))

  up_df <- plot_df[plot_df$logFC > 0, , drop = FALSE]
  down_df <- plot_df[plot_df$logFC < 0, , drop = FALSE]
  up_labels <- head(up_df[order(-up_df$priority_score, up_df$P.Value, -abs(up_df$logFC)), , drop = FALSE], top_n_each)
  down_labels <- head(down_df[order(-down_df$priority_score, down_df$P.Value, -abs(down_df$logFC)), , drop = FALSE], top_n_each)
  label_df <- rbind(up_labels, down_labels)

  p <- ggplot(plot_df, aes(x = logFC, y = neg_log10_pvalue)) +
    geom_point(aes(color = direction), alpha = 0.55, size = 1.2) +
    scale_color_manual(values = c(up = "#c0392b", down = "#1f618d", flat = "grey60")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    theme_bw(base_size = 11) +
    labs(
      title = paste0("Transcript DTU Volcano: ", contrast_name),
      subtitle = "Labels prioritize combined effect size and statistical significance",
      x = "logFC",
      y = expression(-log[10](P.Value)),
      color = "Direction"
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  if (nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 3,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf,
      seed = 1
    )
  }

  ggplot2::ggsave(output_file, plot = p, width = 11, height = 8.5, units = "in")
}

matrix_file <- normalizePath(opts[["--matrix"]], mustWork = TRUE)
samples_file <- normalizePath(opts[["--samples"]], mustWork = TRUE)
contrasts_file <- normalizePath(opts[["--contrasts"]], mustWork = TRUE)
outdir <- opts[["--outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)

log_msg("Reading samples from ", samples_file)
samples_df <- read_no_header_tsv(samples_file, c("condition", "sample"))
if (anyDuplicated(samples_df$sample)) {
  dupes <- unique(samples_df$sample[duplicated(samples_df$sample)])
  stop("Duplicate sample names in samples file: ", paste(dupes, collapse = ", "), call. = FALSE)
}

condition_levels <- unique(samples_df$condition)
samples_df$condition <- factor(samples_df$condition, levels = condition_levels)

log_msg("Reading contrasts from ", contrasts_file)
contrasts_df <- read_no_header_tsv(contrasts_file, c("numerator", "denominator"))
if (nrow(contrasts_df) < 1) {
  stop("No contrasts found in ", contrasts_file, call. = FALSE)
}

unknown_conditions <- setdiff(unique(unlist(contrasts_df)), condition_levels)
if (length(unknown_conditions) > 0) {
  stop(
    "Contrast file contains conditions absent from samples file: ",
    paste(unknown_conditions, collapse = ", "),
    call. = FALSE
  )
}

log_msg("Reading count matrix from ", matrix_file)
count_df <- read.delim(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("gene_id", "transcript_id")
missing_required <- setdiff(required_cols, names(count_df))
if (length(missing_required) > 0) {
  stop("Matrix is missing required columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

missing_samples <- setdiff(samples_df$sample, names(count_df))
if (length(missing_samples) > 0) {
  stop("Samples missing from matrix: ", paste(missing_samples, collapse = ", "), call. = FALSE)
}

count_df <- count_df[, c(required_cols, samples_df$sample), drop = FALSE]
counts <- as.matrix(count_df[, samples_df$sample, drop = FALSE])
mode(counts) <- "numeric"
rownames(counts) <- count_df$transcript_id

# Determine the grouping key that defines a "gene" for DTU testing.
# Default groups by the full LRAA gene_id. With --gene-grouping gene_symbol,
# transcripts are grouped by gene symbol, i.e. the substring of the gene_id
# preceding the first '^' (LRAA encodes ids as
# "<symbol>^g:<chr>:<strand>:comp-<n>"). Grouping by symbol pools transcripts
# from loci that share a symbol but were split across separate gene_id
# components.
if (identical(opts[["--gene-grouping"]], "gene_symbol")) {
  count_df$group_id <- sub("\\^.*$", "", count_df$gene_id)
} else {
  count_df$group_id <- count_df$gene_id
}
log_msg(
  "Gene grouping key: ", opts[["--gene-grouping"]],
  " (", length(unique(count_df$group_id)), " groups across ",
  length(unique(count_df$gene_id)), " gene_ids)"
)

# Reporting-only fractions are computed before filtering so downstream
# summaries reflect the full observed transcript mixture for each gene.
full_gene_totals <- rowsum(counts, group = count_df$group_id, reorder = FALSE)
full_gene_index <- match(count_df$group_id, rownames(full_gene_totals))
full_gene_totals_by_tx <- full_gene_totals[full_gene_index, , drop = FALSE]
full_isoform_fraction <- counts / full_gene_totals_by_tx
full_isoform_fraction[!is.finite(full_isoform_fraction)] <- 0

design <- model.matrix(~0 + condition, data = samples_df)
colnames(design) <- sub("^condition", "", colnames(design))

y <- DGEList(
  counts = counts,
  genes = data.frame(
    group_id = count_df$group_id,
    gene_id = count_df$gene_id,
    transcript_id = count_df$transcript_id,
    stringsAsFactors = FALSE,
    row.names = count_df$transcript_id
  )
)

lib_sizes_raw <- colSums(y$counts)
sample_summary <- data.frame(
  sample = samples_df$sample,
  condition = as.character(samples_df$condition),
  lib_size_raw = as.numeric(lib_sizes_raw[samples_df$sample]),
  stringsAsFactors = FALSE
)

# Transcript filtering happens in three stages:
#   1. edgeR::filterByExpr for design-aware low-count removal
#   2. explicit CPM prevalence threshold
#   3. retention only of genes that still have >= min-gene-transcripts
# The DTU model is fit only on this filtered feature space.
keep_expr <- filterByExpr(y, design = design)
keep_cpm <- rowSums(edgeR::cpm(y) >= opts[["--min-cpm"]]) >= opts[["--min-samples"]]
keep_expr <- keep_expr & keep_cpm
y <- y[keep_expr, , keep.lib.sizes = FALSE]

tx_per_gene <- table(y$genes$group_id)
keep_gene <- y$genes$group_id %in% names(tx_per_gene[tx_per_gene >= opts[["--min-gene-transcripts"]]])
y <- y[keep_gene, , keep.lib.sizes = FALSE]

if (nrow(y) == 0) {
  stop("No transcripts remain after filtering", call. = FALSE)
}

y <- normLibSizes(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

# These summaries are used only to annotate transcript-level result tables.
# They are not inputs to the statistical test itself.
observed_cpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = FALSE)
observed_logcpm <- edgeR::cpm(y, normalized.lib.sizes = TRUE, log = TRUE)
fitted_values <- fit$fitted.values
rownames(fitted_values) <- rownames(y$counts)
colnames(fitted_values) <- colnames(y$counts)
fitted_logcpm <- edgeR::cpm(fitted_values, normalized.lib.sizes = TRUE, log = TRUE)

sample_summary$lib_size_filtered <- y$samples$lib.size[match(sample_summary$sample, rownames(y$samples))]
sample_summary$norm_factor <- y$samples$norm.factors[match(sample_summary$sample, rownames(y$samples))]
sample_summary$lib_size_effective <- sample_summary$lib_size_filtered * sample_summary$norm_factor

filtering_summary <- data.frame(
  metric = c(
    "transcripts_input",
    "genes_input",
    "transcripts_after_expr_filter",
    "genes_after_expr_filter",
    "transcripts_after_gene_tx_filter",
    "genes_after_gene_tx_filter"
  ),
  value = c(
    nrow(count_df),
    length(unique(count_df$group_id)),
    sum(keep_expr),
    length(unique(count_df$group_id[keep_expr])),
    nrow(y),
    length(unique(y$genes$group_id))
  ),
  stringsAsFactors = FALSE
)

write.table(samples_df, file.path(outdir, "samples.used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(contrasts_df, file.path(outdir, "contrasts.used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sample_summary, file.path(outdir, "sample.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtering_summary, file.path(outdir, "filtering.summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(
  data.frame(transcript_id = y$genes$transcript_id, gene_id = y$genes$gene_id, group_id = y$genes$group_id, stringsAsFactors = FALSE),
  file.path(outdir, "features.retained.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(sample = rownames(design), design, check.names = FALSE),
  file.path(outdir, "design.matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

for (i in seq_len(nrow(contrasts_df))) {
  numerator <- contrasts_df$numerator[[i]]
  denominator <- contrasts_df$denominator[[i]]
  contrast_name <- paste0(sanitize_name(numerator), "_vs_", sanitize_name(denominator))
  log_msg("Running DTU contrast: ", numerator, " vs ", denominator)

  numerator_samples <- samples_df$sample[samples_df$condition == numerator]
  denominator_samples <- samples_df$sample[samples_df$condition == denominator]

  contrast_matrix <- limma::makeContrasts(contrasts = paste0(numerator, "-", denominator), levels = design)
  ds <- diffSplice(
    fit,
    contrast = contrast_matrix,
    geneid = "group_id",
    exonid = "transcript_id",
    verbose = isTRUE(opts[["--verbose"]])
  )

  # edgeR still uses exon-oriented naming in diffSplice because the method was
  # originally developed for exon usage. Here, "exons" correspond to transcript
  # features grouped by the gene grouping key (gene_id or gene_symbol).
  gene_f <- topSplice(ds, test = "F", number = Inf, FDR = opts[["--fdr"]], sort.by = "p")
  gene_simes <- topSplice(ds, test = "simes", number = Inf, FDR = opts[["--fdr"]], sort.by = "p")
  transcript_t <- topSplice(ds, test = "t", number = Inf, FDR = opts[["--fdr"]], sort.by = "p")

  # topSplice carries the group_id grouping column (constant within each group)
  # into the gene-level tables; make it the leading column and drop rownames.
  gene_f <- gene_f[, c("group_id", setdiff(names(gene_f), "group_id")), drop = FALSE]
  gene_simes <- gene_simes[, c("group_id", setdiff(names(gene_simes), "group_id")), drop = FALSE]
  rownames(gene_f) <- NULL
  rownames(gene_simes) <- NULL
  transcript_t <- cbind(feature_id = rownames(transcript_t), transcript_t, row.names = NULL)

  tx_ids <- transcript_t$transcript_id
  tx_index <- match(tx_ids, rownames(y$counts))

  # Add contrast-specific summary columns similar in spirit to the intron-level
  # Diff-Splice-Finder output, but adapted for transcript DTU.
  transcript_t$contrast <- contrast_name
  transcript_t$contrast_group1_mean_count <- rowMeans(y$counts[tx_index, numerator_samples, drop = FALSE])
  transcript_t$contrast_group2_mean_count <- rowMeans(y$counts[tx_index, denominator_samples, drop = FALSE])
  transcript_t$contrast_group1_mean_CPM <- rowMeans(observed_cpm[tx_index, numerator_samples, drop = FALSE])
  transcript_t$contrast_group2_mean_CPM <- rowMeans(observed_cpm[tx_index, denominator_samples, drop = FALSE])
  transcript_t$contrast_group1_mean_logCPM <- rowMeans(observed_logcpm[tx_index, numerator_samples, drop = FALSE])
  transcript_t$contrast_group2_mean_logCPM <- rowMeans(observed_logcpm[tx_index, denominator_samples, drop = FALSE])
  transcript_t$contrast_group1_mean_fitted_logCPM <- rowMeans(fitted_logcpm[tx_index, numerator_samples, drop = FALSE])
  transcript_t$contrast_group2_mean_fitted_logCPM <- rowMeans(fitted_logcpm[tx_index, denominator_samples, drop = FALSE])
  full_tx_index <- match(tx_ids, rownames(full_isoform_fraction))
  # These isoform fractions use full-gene denominators from the sample-subset
  # matrix before filtering. They are descriptive summaries only.
  transcript_t[[paste0(numerator, "_mean_isoform_fraction")]] <- rowMeans(full_isoform_fraction[full_tx_index, numerator_samples, drop = FALSE])
  transcript_t[[paste0(denominator, "_mean_isoform_fraction")]] <- rowMeans(full_isoform_fraction[full_tx_index, denominator_samples, drop = FALSE])
  transcript_t$delta_isoform_fraction <-
    transcript_t[[paste0(numerator, "_mean_isoform_fraction")]] -
    transcript_t[[paste0(denominator, "_mean_isoform_fraction")]]

  prefix <- file.path(outdir, paste0("edgeR.DTU.", contrast_name))
  write.table(gene_f, paste0(prefix, ".gene_F.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(gene_simes, paste0(prefix, ".gene_simes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(transcript_t, paste0(prefix, ".transcript_t.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  make_transcript_volcano_plot(
    transcript_t,
    contrast_name = contrast_name,
    output_file = paste0(prefix, ".transcript_volcano.pdf")
  )
}

writeLines(capture.output(sessionInfo()), con = file.path(outdir, "sessionInfo.txt"))
