#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(Matrix)
  library(DropletUtils)
})

# Argument parser for command-line execution.
parser <- ArgumentParser(description='Filter sparse matrix to retain only good cells using emptyDrops')
parser$add_argument("--matrix_dir", required=TRUE, 
                    help="Input directory containing matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz")
parser$add_argument("--output_dir", required=TRUE, 
                    help="Output directory for filtered matrix")
parser$add_argument("--fdr_threshold", type="double", default=0.01, 
                    help="FDR threshold for emptyDrops (default: 0.01)")
parser$add_argument("--lower", type="integer", default=NULL, 
                    help="Lower UMI count threshold for emptyDrops (default: auto)")
parser$add_argument("--isoform_matrix_dir", default=NULL, 
                    help="Optional: Input directory for isoform sparse matrix to filter with same cells")
parser$add_argument("--isoform_output_dir", default=NULL, 
                    help="Optional: Output directory for filtered isoform matrix")
parser$add_argument("--splice_pattern_matrix_dir", default=NULL, 
                    help="Optional: Input directory for splice pattern sparse matrix to filter with same cells")
parser$add_argument("--splice_pattern_output_dir", default=NULL, 
                    help="Optional: Output directory for filtered splice pattern matrix")
args <- parser$parse_args()

# Extract arguments
matrix_dir <- args$matrix_dir
output_dir <- args$output_dir
fdr_threshold <- args$fdr_threshold
lower_threshold <- args$lower
isoform_matrix_dir <- args$isoform_matrix_dir
isoform_output_dir <- args$isoform_output_dir
splice_pattern_matrix_dir <- args$splice_pattern_matrix_dir
splice_pattern_output_dir <- args$splice_pattern_output_dir

# Validate input directory
if (!dir.exists(matrix_dir)) {
  stop("Input directory does not exist: ", matrix_dir)
}

# Create output directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

# Read sparse matrix files
message("Reading sparse matrix from: ", matrix_dir)
matrix_file <- file.path(matrix_dir, "matrix.mtx.gz")
features_file <- file.path(matrix_dir, "features.tsv.gz")
barcodes_file <- file.path(matrix_dir, "barcodes.tsv.gz")

# Check if files exist
if (!file.exists(matrix_file)) {
  # Try without .gz extension
  matrix_file <- file.path(matrix_dir, "matrix.mtx")
  if (!file.exists(matrix_file)) {
    stop("matrix.mtx or matrix.mtx.gz not found in: ", matrix_dir)
  }
}
if (!file.exists(features_file)) {
  features_file <- file.path(matrix_dir, "features.tsv")
  if (!file.exists(features_file)) {
    stop("features.tsv or features.tsv.gz not found in: ", matrix_dir)
  }
}
if (!file.exists(barcodes_file)) {
  barcodes_file <- file.path(matrix_dir, "barcodes.tsv")
  if (!file.exists(barcodes_file)) {
    stop("barcodes.tsv or barcodes.tsv.gz not found in: ", matrix_dir)
  }
}

# Read the sparse matrix
sparse_matrix <- readMM(matrix_file)
features <- readLines(features_file)
barcodes <- readLines(barcodes_file)

# Set dimension names
rownames(sparse_matrix) <- features
colnames(sparse_matrix) <- barcodes

message("Loaded matrix: ", nrow(sparse_matrix), " features x ", ncol(sparse_matrix), " barcodes")
message("Total UMI count: ", sum(sparse_matrix))

# Run emptyDrops with error handling
message("Running emptyDrops with FDR threshold: ", fdr_threshold)
droplet_measure <- tryCatch({
  if (!is.null(lower_threshold)) {
    message("Using lower threshold: ", lower_threshold)
    emptyDrops(sparse_matrix, lower = lower_threshold)
  } else {
    emptyDrops(sparse_matrix)
  }
}, error = function(e) {
  if (grepl("no counts available", e$message)) {
    message("Default emptyDrops failed; retrying with lower = 0.")
    emptyDrops(sparse_matrix, lower = 0)
  } else {
    stop(e)
  }
})

# Identify good cells based on FDR threshold
good_cells_idx <- which(droplet_measure$FDR <= fdr_threshold)
n_good_cells <- length(good_cells_idx)

if (n_good_cells == 0) {
  stop("No cells passed the FDR threshold of ", fdr_threshold, 
       ". Consider adjusting --fdr_threshold or --lower parameters.")
}

message("Found ", n_good_cells, " good cells (", 
        round(100 * n_good_cells / ncol(sparse_matrix), 2), "% of total)")

# Filter the sparse matrix
filtered_matrix <- sparse_matrix[, good_cells_idx]
filtered_barcodes <- barcodes[good_cells_idx]

message("Filtered matrix: ", nrow(filtered_matrix), " features x ", 
        ncol(filtered_matrix), " barcodes")
message("Filtered UMI count: ", sum(filtered_matrix))

# Write filtered matrix
message("Writing filtered matrix to: ", output_dir)
writeMM(obj = filtered_matrix, file = file.path(output_dir, "matrix.mtx"))
write(features, file = file.path(output_dir, "features.tsv"))
write(filtered_barcodes, file = file.path(output_dir, "barcodes.tsv"))

# Compress output files
message("Compressing output files...")
system(paste0("gzip -f ", file.path(output_dir, "matrix.mtx")))
system(paste0("gzip -f ", file.path(output_dir, "features.tsv")))
system(paste0("gzip -f ", file.path(output_dir, "barcodes.tsv")))

# Write summary statistics
summary_stats <- data.frame(
  Metric = c("Input Barcodes", "Input UMI Count", 
             "Filtered Barcodes", "Filtered UMI Count",
             "Fraction Barcodes Retained", "Fraction UMI Retained",
             "FDR Threshold"),
  Value = c(ncol(sparse_matrix), sum(sparse_matrix),
            ncol(filtered_matrix), sum(filtered_matrix),
            round(ncol(filtered_matrix) / ncol(sparse_matrix), 4),
            round(sum(filtered_matrix) / sum(sparse_matrix), 4),
            fdr_threshold)
)

summary_file <- file.path(output_dir, "filtering_summary.tsv")
write.table(summary_stats, file = summary_file, sep = "\t", 
            row.names = FALSE, quote = FALSE)
message("Summary statistics written to: ", summary_file)

# Print summary to console
message("\nFiltering Summary:")
print(summary_stats)

message("\nFiltering complete!")

# Helper function to filter additional sparse matrices using the same good cell barcodes
filter_additional_matrix <- function(input_dir, output_dir, good_barcodes) {
  message("\n--- Filtering additional sparse matrix ---")
  message("Input directory: ", input_dir)
  message("Output directory: ", output_dir)
  
  # Read matrix files
  matrix_file <- file.path(input_dir, "matrix.mtx.gz")
  features_file <- file.path(input_dir, "features.tsv.gz")
  barcodes_file <- file.path(input_dir, "barcodes.tsv.gz")
  
  # Check if files exist (try without .gz extension if not found)
  if (!file.exists(matrix_file)) {
    matrix_file <- file.path(input_dir, "matrix.mtx")
    if (!file.exists(matrix_file)) {
      stop("matrix.mtx or matrix.mtx.gz not found in: ", input_dir)
    }
  }
  if (!file.exists(features_file)) {
    features_file <- file.path(input_dir, "features.tsv")
    if (!file.exists(features_file)) {
      stop("features.tsv or features.tsv.gz not found in: ", input_dir)
    }
  }
  if (!file.exists(barcodes_file)) {
    barcodes_file <- file.path(input_dir, "barcodes.tsv")
    if (!file.exists(barcodes_file)) {
      stop("barcodes.tsv or barcodes.tsv.gz not found in: ", input_dir)
    }
  }
  
  # Read the sparse matrix
  sparse_matrix <- readMM(matrix_file)
  features <- readLines(features_file)
  barcodes <- readLines(barcodes_file)
  
  # Set dimension names
  rownames(sparse_matrix) <- features
  colnames(sparse_matrix) <- barcodes
  
  message("Loaded matrix: ", nrow(sparse_matrix), " features x ", ncol(sparse_matrix), " barcodes")
  
  # Find indices of good cells in this matrix
  good_cells_idx <- which(barcodes %in% good_barcodes)
  
  if (length(good_cells_idx) == 0) {
    warning("No matching good cell barcodes found in this matrix!")
    return(NULL)
  }
  
  message("Found ", length(good_cells_idx), " matching good cells in this matrix")
  
  # Filter the sparse matrix to keep only good cells
  filtered_matrix <- sparse_matrix[, good_cells_idx]
  filtered_barcodes <- barcodes[good_cells_idx]
  
  message("Filtered matrix: ", nrow(filtered_matrix), " features x ", ncol(filtered_matrix), " barcodes")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  # Write filtered matrix
  message("Writing filtered matrix to: ", output_dir)
  writeMM(obj = filtered_matrix, file = file.path(output_dir, "matrix.mtx"))
  write(features, file = file.path(output_dir, "features.tsv"))
  write(filtered_barcodes, file = file.path(output_dir, "barcodes.tsv"))
  
  # Compress output files
  message("Compressing output files...")
  system(paste0("gzip -f ", file.path(output_dir, "matrix.mtx")))
  system(paste0("gzip -f ", file.path(output_dir, "features.tsv")))
  system(paste0("gzip -f ", file.path(output_dir, "barcodes.tsv")))
  
  message("Additional matrix filtering complete!")
}

# Filter isoform matrix if provided
if (!is.null(isoform_matrix_dir) && !is.null(isoform_output_dir)) {
  message("\n========================================")
  message("Filtering isoform sparse matrix")
  message("========================================")
  filter_additional_matrix(isoform_matrix_dir, isoform_output_dir, filtered_barcodes)
}

# Filter splice pattern matrix if provided
if (!is.null(splice_pattern_matrix_dir) && !is.null(splice_pattern_output_dir)) {
  message("\n========================================")
  message("Filtering splice pattern sparse matrix")
  message("========================================")
  filter_additional_matrix(splice_pattern_matrix_dir, splice_pattern_output_dir, filtered_barcodes)
}

message("\n========================================")
message("All filtering operations complete!")
message("========================================")
