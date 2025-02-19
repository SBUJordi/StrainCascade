#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# genome_annotation_aggregation.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(qs)
  library(optparse)
  library(dplyr)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character", default=NULL, help="Output directory for results"),
  make_option(c("--bakta"), type="character", default=NULL, help="Path to the Bakta results file"),
  make_option(c("--prokka"), type="character", default=NULL, help="Path to the Prokka results file"),
  make_option(c("--microbeannotator"), type="character", default=NULL, help="Path to the MicrobeAnnotator results file"),
  make_option(c("--version"), type="character", default=NULL, help="StrainCascade version")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to merge data frames by common columns
merge_by_common_columns <- function(df1, df2) {
  common_columns <- intersect(colnames(df1), colnames(df2))
  
  if (length(common_columns) == 0) {
    stop("No common columns to merge by.")
  }
  
  merged_df <- merge(df1, df2, by = common_columns, all = TRUE)
  return(merged_df)
}

# Function to remove rows where all specified columns are NA
remove_na_rows <- function(df, pattern) {
  source_cols <- grep(pattern, colnames(df), value = TRUE)
  df <- df[!apply(df[source_cols], 1, function(x) all(is.na(x))), ]
  return(df)
}

# Main execution
tryCatch({
  # Check which input files are available and import them
  results <- list()
  if (!is.null(opt$bakta)) results$bakta <- qread(file.path(opt$output_dir, opt$bakta))
  if (!is.null(opt$prokka)) results$prokka <- qread(file.path(opt$output_dir, opt$prokka))
  if (!is.null(opt$microbeannotator)) results$microbeannotator <- qread(file.path(opt$output_dir, opt$microbeannotator))
  
  # Check if any files were imported
  if (length(results) == 0) {
    cat("INFO: No input files available. Exiting script.\n")
    quit(save = "no", status = 0)
  }
  
  # Merge available data frames
  if (length(results) == 1) {
    annotation_results_aggregated <- results[[1]]  # Only one file available
  } else {
    # Merge all available data frames
    annotation_results_aggregated <- Reduce(merge_by_common_columns, results)
    
    # Remove rows where all columns containing 'source_file' in their name are NA
    annotation_results_aggregated <- remove_na_rows(annotation_results_aggregated, "source_file")
    
    # Remove duplicate rows
    annotation_results_aggregated <- unique(annotation_results_aggregated)
  }
  
  # Add StrainCascade version if not already present
  if (!"software_version" %in% colnames(annotation_results_aggregated)) {
    annotation_results_aggregated$software_version <- paste("StrainCascade v", opt$version, sep="")
  }
  
  # Save the final result
  output_file <- file.path(opt$output_dir, "annotation_results_aggregated.qs")
  qsave(annotation_results_aggregated, output_file, preset = "archive")
  
  cat("SUCCESS: Annotation results processed and saved successfully to", output_file, "\n")
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})