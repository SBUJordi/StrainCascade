#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_circlator.R

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(optparse)
})

# Function to replace "NA" strings with actual NA values, excluding date-time columns
replace_na_strings <- function(df) {
  # Identify date-time columns
  date_cols <- sapply(df, inherits, what = c("POSIXt", "Date"))
  
  # Define a pattern to match various representations of NA
  na_pattern <- "^(NA|N/A|N/a|n/A|n/a|)$"
  
  # Perform replacement on non-date-time columns
  df[!date_cols] <- df[!date_cols] %>%
    mutate(across(everything(), ~ replace(.x, grepl(na_pattern, .x), NA)))
  
  return(df)
}

# Define command line options
option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="CHARACTER"),
  make_option(c("-l", "--log_file"), type="character", default=NULL, 
              help="Input Circlator log file path", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$log_file) || is.null(opt$version)) {
  stop("All arguments (output_dir, log_file, and version) must be supplied.")
}

# Read and process the Circlator log file
circlator_log <- read.delim(opt$log_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
circlator_log <- replace_na_strings(circlator_log)

# Rename columns to match the desired format
colnames(circlator_log) <- c(
  "modus_circlator",
  "contig_tag_fasta",
  "repetitive_deleted_circlator",
  "circularised_using_nucmer_circlator",
  "circularised_using_spades_circlator",
  "circularisation_status"
)

# Create the circlator dataframe
circlator <- circlator_log %>%
  mutate(
    repetitive_deleted_circlator = as.integer(repetitive_deleted_circlator),
    circularised_using_nucmer_circlator = as.integer(circularised_using_nucmer_circlator),
    circularised_using_spades_circlator = as.integer(circularised_using_spades_circlator),
    circularisation_status = as.integer(circularisation_status)
  )

circlator$software_version <- paste0("StrainCascade v", opt$version)
circlator$source_file_circlator <- basename(opt$log_file)

# Save as qs
output_file <- file.path(opt$output_dir, "circlator_results.qs")
qsave(circlator, output_file, preset = "archive")

cat("Circlator results processed and saved to:", output_file, "\n")