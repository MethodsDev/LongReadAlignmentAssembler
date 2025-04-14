#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(tidyr)
  library(dplyr)
  library(Matrix)
  library(DropletUtils)
  library(ggplot2)
  library(scales)
  library(gridExtra)  # For table grob and arranging plots
})

# Argument parser for Linux command-line execution.
parser <- ArgumentParser(description='Sparse Matrix Outputs Pipeline')
parser$add_argument("--working_dir", required=TRUE, help="Working directory path")
parser$add_argument("--track_file", required=TRUE, help="Tracking file (e.g., BrCa_NeoMET_26b.LRAA.quant.tracking)")
parser$add_argument("--bstats_file", required=FALSE, default=NULL, help="BCStats file path")
parser$add_argument("--sample_name", required=TRUE, help="Sample name")
parser$add_argument("--create_filtered", type="integer", default=1,
                    help="1 to create filtered matrix and related QC plots, 0 to generate raw matrix only")
parser$add_argument("--fdr_threshold", type="double", default=0.01, help="FDR threshold")
parser$add_argument("--umistats_file", required=FALSE, default=NULL, help="UMI stats file path")
args <- parser$parse_args()

# Working directory and assign variables from parsed arguments.
setwd(args$working_dir)
dat_filename_track  <- args$track_file
dat_filename_bstats <- args$bstats_file
sample_name         <- args$sample_name
create_filtered     <- (args$create_filtered == 1)
fdr_threshold       <- args$fdr_threshold
umistats_file       <- args$umistats_file

# Read the tracking data.
message("-reading in ", dat_filename_track)
data_track <- fread(dat_filename_track, sep="\t", header=TRUE)

# Extract cell barcode, UMI, and core read name.
data_track_mod <- data_track %>% 
  separate_wider_delim(read_name, "^", names=c("cell_barcode", "UMI", "core_read_name"))

# Get gene/cell counts.
gene_cell_counts <- data_track_mod %>% 
  group_by(gene_id, cell_barcode) %>% 
  summarize(sum_UMIs = sum(frac_assigned)) %>% 
  ungroup()
gene_cell_counts_tsv <- paste0(sample_name, ".gene_cell_counts.tsv")
fwrite(gene_cell_counts, file=gene_cell_counts_tsv, sep="\t", row.names=FALSE, quote=FALSE)
message("-writing gene cell counts table: ", gene_cell_counts_tsv)

# Get isoform/cell counts.
isoform_cell_counts <- data_track_mod %>% 
  group_by(transcript_id, cell_barcode) %>% 
  summarize(sum_UMIs = sum(frac_assigned)) %>% 
  ungroup()
isoform_cell_counts_tsv <- paste0(sample_name, ".isoform_cell_counts.tsv")
fwrite(isoform_cell_counts, file=isoform_cell_counts_tsv, sep="\t", row.names=FALSE, quote=FALSE)
message("-writing isoform cell counts table: ", isoform_cell_counts_tsv)

# Function to create sparse matrix outputs with an adjustable FDR threshold.
# Additional parameters:
#   - bstats_file_path: if provided, additional reads-based summary stats are computed.
#   - umistats_df: if provided (as a data frame or file path), sequencing saturation is calculated.
#   - create_filtered: if TRUE, the filtered matrix and related QC plots are created;
#                      if FALSE, only the raw sparse matrix and its associated QC plots/summary are generated.
make_sparse_matrix_outputs <- function(counts_data, outdirname, fdr_threshold = 0.01,
                                       bstats_file_path = NULL, umistats_df = NULL,
                                       create_filtered = TRUE) {
  message("-making sparse matrix outputs for: ", outdirname)
  
  # If umistats_df is provided as a file path, load it.
  if (!is.null(umistats_df) && is.character(umistats_df)) {
    umistats_df <- fread(umistats_df, sep="\t", header=TRUE)
  }
  
  # Create output directories.
  raw_dir <- file.path(outdirname, "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
  if (create_filtered) {
    filtered_dir <- file.path(outdirname, "filtered_feature_bc_matrix")
    if (!dir.exists(filtered_dir)) dir.create(filtered_dir, recursive = TRUE)
  }
  
  # Ensure counts_data has expected column names.
  colnames(counts_data) <- c("feature_id", "cell_barcode", "UMI_counts")
  feature_ids <- counts_data %>% select(feature_id) %>% unique() %>% pull(feature_id)
  cell_barcodes <- counts_data %>% select(cell_barcode) %>% unique() %>% pull(cell_barcode)
  
  # Factorize identifiers.
  feature_factor <- factor(counts_data$feature_id, levels = feature_ids)
  cell_factor <- factor(counts_data$cell_barcode, levels = cell_barcodes)
  i <- as.numeric(feature_factor)
  j <- as.numeric(cell_factor)
  
  # Create the raw sparse matrix.
  sparseM <- Matrix::sparseMatrix(i = i, j = j, x = counts_data$UMI_counts,
                                  dims = c(length(feature_ids), length(cell_barcodes)),
                                  dimnames = list(feature_ids, cell_barcodes))
  
  # Write raw matrix and identifiers.
  Matrix::writeMM(obj = sparseM, file = file.path(raw_dir, "matrix.mtx"))
  write(feature_ids, file = file.path(raw_dir, "features.tsv"))
  write(cell_barcodes, file = file.path(raw_dir, "barcodes.tsv"))
  
  # Run emptyDrops with error handling.
  droplet_measure <- tryCatch({
    emptyDrops(sparseM)
  }, error = function(e) {
    if (grepl("no counts available", e$message)) {
      message("Default emptyDrops failed; retrying with lower = 0.")
      emptyDrops(sparseM, lower = 0)
    } else stop(e)
  })
  
  # If filtering is requested, perform filtering.
  if (create_filtered) {
    real_droplets <- which(droplet_measure$FDR <= fdr_threshold)
    filtered_sparseM <- sparseM[, real_droplets]
    filtered_cell_barcodes <- colnames(filtered_sparseM)
    
    # Write filtered matrix.
    Matrix::writeMM(obj = filtered_sparseM, file = file.path(filtered_dir, "matrix.mtx"))
    write(feature_ids, file = file.path(filtered_dir, "features.tsv"))
    write(filtered_cell_barcodes, file = file.path(filtered_dir, "barcodes.tsv"))
  } else {
    filtered_sparseM <- NA
    filtered_cell_barcodes <- cell_barcodes
  }
  
  # Compress output files.
  system(paste0("gzip ", file.path(raw_dir, "*")))
  if (create_filtered) {
    system(paste0("gzip ", file.path(filtered_dir, "*")))
  }
  
  # Create QC plots directory.
  qc_dir <- file.path(outdirname, paste0(sample_name, "^QC_plots"))
  if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
  
  ## Knee Plot Calculation ##
  br.out <- barcodeRanks(sparseM)
  if (create_filtered) {
    valid_idx <- which(!is.na(droplet_measure$FDR) & droplet_measure$FDR <= fdr_threshold)
    valid_br <- br.out[rownames(br.out) %in% colnames(sparseM)[valid_idx], ]
  } else {
    valid_idx <- seq_along(cell_barcodes)
    valid_br <- br.out
  }
  
  # Compute basic summary statistics.
  all_barcodes_count <- length(cell_barcodes)
  total_umi <- sum(sparseM)
  fdr_count <- length(valid_idx)
  valid_total_umi <- sum(droplet_measure$Total[valid_idx])
  median_umi_all <- median(droplet_measure$Total, na.rm = TRUE)
  median_umi_filtered <- median(droplet_measure$Total[valid_idx], na.rm = TRUE)
  
  if (!is.null(bstats_file_path)) {
    bstats_file <- fread(bstats_file_path, sep="\t", header=TRUE)
    all_bstats <- bstats_file
    if (create_filtered) {
      filtered_bstats <- bstats_file[`#BarcodeSequence` %in% filtered_cell_barcodes]
    } else {
      filtered_bstats <- bstats_file[`#BarcodeSequence` %in% cell_barcodes]
    }
    
    mean_reads_all <- mean(all_bstats$NumberOfReads, na.rm = TRUE)
    total_reads_all <- sum(all_bstats$NumberOfReads, na.rm = TRUE)
    mean_reads_filtered <- mean(filtered_bstats$NumberOfReads, na.rm = TRUE)
    total_reads_filtered <- sum(filtered_bstats$NumberOfReads, na.rm = TRUE)
    fraction_reads <- total_reads_filtered / total_reads_all
    
    if (create_filtered) {
      summary_df <- data.frame(
        Metric = c("Total barcodes #", "Total UMI #", "Real Droplet Barcodes #", "Real Droplet UMI #",
                   "Total Reads (All)", "Mean Reads (All)", "Total Reads (Filtered)", "Mean Reads (Filtered)",
                   "Fraction Reads (Filtered/All)", "Median UMI per Cell (All)", "Median UMI per Cell (Filtered)"),
        Value = c(all_barcodes_count, total_umi, fdr_count, valid_total_umi,
                  total_reads_all, round(mean_reads_all, 2), total_reads_filtered, round(mean_reads_filtered, 2),
                  round(fraction_reads, 4), median_umi_all, median_umi_filtered)
      )
    } else {
      summary_df <- data.frame(
        Metric = c("Total barcodes #", "Total UMI #", "Total Reads (All)", "Mean Reads (All)", 
                   "Median UMI per Cell (All)"),
        Value = c(all_barcodes_count, total_umi, total_reads_all, round(mean_reads_all, 2), median_umi_all)
      )
    }
  } else {
    if (create_filtered) {
      summary_df <- data.frame(
        Metric = c("Total barcodes #", "Total UMI #", "Real Droplet Barcodes #", "Real Droplet UMI #",
                   "Median UMI per Cell (All)", "Median UMI per Cell (Filtered)"),
        Value = c(all_barcodes_count, total_umi, fdr_count, valid_total_umi,
                  median_umi_all, median_umi_filtered)
      )
    } else {
      summary_df <- data.frame(
        Metric = c("Total barcodes #", "Total UMI #", "Median UMI per Cell (All)"),
        Value = c(all_barcodes_count, total_umi, median_umi_all)
      )
    }
  }
  
  ### Sequencing Saturation Calculation and Plot (if umistats_df is provided) ###
  if (!is.null(umistats_df)) {
    umistats_dt <- as.data.table(umistats_df)
    # Global saturation for all cells.
    all_cell_summary <- umistats_dt[, .(total_reads = sum(Count), unique_UMIs = .N), by = cell_barcode]
    global_saturation_all <- 1 - (sum(all_cell_summary$unique_UMIs) / sum(all_cell_summary$total_reads))
    # Global saturation for filtered cells.
    filtered_cell_summary <- umistats_dt[cell_barcode %in% (if(create_filtered) filtered_cell_barcodes else cell_barcodes), 
                                         .(total_reads = sum(Count), unique_UMIs = .N), by = cell_barcode]
    global_saturation_filtered <- 1 - (sum(filtered_cell_summary$unique_UMIs) / sum(filtered_cell_summary$total_reads))
    # Per-cell saturation for filtered cells.
    cell_summary <- filtered_cell_summary
    cell_summary[, saturation := 1 - (unique_UMIs / total_reads)]
    # Fit Michaelis-Menten model.
    fit <- nls(saturation ~ 1 - (Vmax/(K + total_reads)), data = cell_summary,
               start = list(Vmax = max(cell_summary$unique_UMIs), K = 1e4))
    newdata <- data.frame(total_reads = seq(min(cell_summary$total_reads), max(cell_summary$total_reads), length.out = 100))
    newdata$sat_fit <- predict(fit, newdata)
    # Create saturation plot.
    sat_plot <- ggplot(cell_summary, aes(x = total_reads, y = saturation)) +
      geom_point() +
      geom_line(data = newdata, aes(x = total_reads, y = sat_fit), color = "red", size = 1) +
      labs(
        x = "Total Reads per Cell",
        y = "Sequencing Saturation",
        title = "Sequencing Saturation vs Total Reads (Michaelis-Menten Fit)"
      ) +
      theme_classic() +
      scale_y_continuous(
        limits = c(0, 1),
        labels = percent_format()          # Format 0..1 as 0%..100%
      ) +
      scale_x_continuous(
        limits = c(0, max(cell_summary$total_reads, na.rm = TRUE)),
        labels = label_number(scale_cut = cut_short_scale()), 
        breaks = pretty_breaks(5)          # About 5 major tick marks
      ) +
      annotate(
        "text",
        x = Inf, 
        y = Inf, 
        label = paste("Global Saturation (All):", round(global_saturation_all, 4)),
        hjust = 1.1, vjust = 2, size = 4
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste("Global Saturation (Filtered):", round(global_saturation_filtered, 4)),
        hjust = 1.1, vjust = 3.5, size = 4
      )
    sat_plot_file <- file.path(qc_dir, "sequencing_saturation.png")
    ggsave(filename = sat_plot_file, plot = sat_plot, width = 8, height = 6, dpi = 300)
    message("Sequencing saturation plot saved to: ", sat_plot_file)
    
    summary_df <- rbind(summary_df,
                        data.frame(Metric = "Global Saturation (All)", Value = global_saturation_all),
                        data.frame(Metric = "Global Saturation (Filtered)", Value = global_saturation_filtered))
  }
  
  ### Export Summary Statistics ###
  summary_csv_file <- file.path(qc_dir, paste0(basename(outdirname), "_summary_statistics.csv"))
  write.csv(summary_df, file = summary_csv_file, row.names = FALSE)
  summary_df$Value <- sapply(summary_df$Value, function(x) format(x, scientific = FALSE, trim = TRUE))
  summary_pdf_file <- file.path(qc_dir, paste0(basename(outdirname), "_summary_statistics.pdf"))
  pdf(summary_pdf_file, width = 8, height = 6)
  grid.arrange(tableGrob(summary_df))
  dev.off()
  
  ### Knee Plot Creation ###
  knee_plot_filename <- file.path(qc_dir, paste0(basename(outdirname), "_barcoderank_realdroplet_threshold_FDR", fdr_threshold, ".png"))
  png(filename = knee_plot_filename, width = 800, height = 600)
  plot(br.out$rank, br.out$total, log = "xy", 
       xlab = "Barcode Rank", ylab = "Total UMI count",
       main = paste("Knee Plot:", basename(outdirname)), pch = 20)
  if(nrow(valid_br) > 0) {
    points(valid_br$rank, valid_br$total, col = "red", pch = 20)
  }
  abline(h = metadata(br.out)$knee, col = "dodgerblue", lty = 2)
  abline(h = metadata(br.out)$inflection, col = "forestgreen", lty = 2)
  valid_threshold_total <- min(droplet_measure$Total[droplet_measure$FDR <= fdr_threshold], na.rm = TRUE)
  abline(h = valid_threshold_total, col = "black", lty = 2)
  valid_threshold_rank <- if(nrow(valid_br) > 0) max(valid_br$rank, na.rm = TRUE) else NA
  if(!is.na(valid_threshold_rank)) abline(v = valid_threshold_rank, col = "black", lty = 2)
  
  # Add legend for abline annotations in the bottom left.
  legend("bottomleft", 
         legend = c("knee", "inflection", "FDR <= 0.01"), 
         col = c("dodgerblue", "forestgreen", "red"), 
         lty = c(2, 2, NA), 
         pch = c(NA, NA, 20))
  
  if (!is.null(bstats_file_path)) {
    legend_text <- c("Summary",
                     paste("Total barcodes #:", all_barcodes_count),
                     paste("Total UMI #:", total_umi),
                     paste("Real Droplet Barcodes #:", fdr_count),
                     paste("Real Droplet UMI #:", valid_total_umi),
                     paste("Total Reads (All):", total_reads_all),
                     paste("Mean Reads (All):", round(mean_reads_all, 2)),
                     paste("Total Reads (Filtered):", total_reads_filtered),
                     paste("Mean Reads (Filtered):", round(mean_reads_filtered, 2)),
                     paste("Fraction Reads (Filtered/All):", round(fraction_reads, 4)),
                     paste("Median UMI per Cell (All):", median_umi_all),
                     paste("Median UMI per Cell (Filtered):", median_umi_filtered))
  } else {
    legend_text <- c("Summary",
                     paste("Total barcodes #:", all_barcodes_count),
                     paste("Total UMI #:", total_umi),
                     paste("Real Droplet Barcodes #:", fdr_count),
                     paste("Real Droplet UMI #:", valid_total_umi),
                     paste("Median UMI per Cell (All):", median_umi_all),
                     paste("Median UMI per Cell (Filtered):", median_umi_filtered))
  }
  legend("topright", legend = legend_text, bty = "n")
  dev.off()
  
  # Create emptyDrops probability plot only if filtering is enabled.
  if (create_filtered) {
    empty_drop_plot <- file.path(qc_dir, paste0(basename(outdirname), "_emptydrops_probability_UMI.png"))
    png(filename = empty_drop_plot, width = 800, height = 600)
    plot(droplet_measure$Total, -droplet_measure$LogProb,
         col = ifelse(droplet_measure$FDR <= fdr_threshold, "red", "black"),
         xlab = "Total UMI count", ylab = "-Log Probability")
    title(basename(outdirname))
    dev.off()
  }
  
  message("Summary Statistics:")
  print(summary_df)
  message("done with ", outdirname)
  
  return(list(raw = sparseM,
              filtered = if (create_filtered) filtered_sparseM else NA,
              droplet_info = droplet_measure,
              summary_stats = summary_df))
}

# Example function calls.
gene_count_list <- make_sparse_matrix_outputs(gene_cell_counts, 
                                               paste0(sample_name, "^gene-sparseM"),
                                               fdr_threshold = 0.01,
                                               bstats_file_path = dat_filename_bstats,
                                               umistats_df = umistats_file,
                                               create_filtered = TRUE)
isoform_count_list <- make_sparse_matrix_outputs(isoform_cell_counts, 
                                                  paste0(sample_name, "^isoform-sparseM"),
                                                  fdr_threshold = 0.01,
                                                  bstats_file_path = dat_filename_bstats,
                                                  umistats_df = umistats_file,
                                                  create_filtered = TRUE)

message("Gene and Isoform matrices Created.")
