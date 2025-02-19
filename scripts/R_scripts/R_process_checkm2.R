#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_checkm2.r

set.seed(42)

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="CHARACTER"),
  make_option(c("-t", "--tsv"), type="character", default=NULL, 
              help="Input TSV file path", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$tsv) || is.null(opt$version)) {
  stop("All arguments (output_dir, tsv, and version) must be supplied.")
}

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

# Read the TSV file
tsv <- read.delim(opt$tsv, header = TRUE, sep = "\t", fill = TRUE)
tsv <- replace_na_strings(tsv)
tsv$source_file_checkm2 <- basename(opt$tsv)

# Process the data
checkm2 <- data.frame(
  completeness = ifelse("Completeness" %in% names(tsv), tsv$Completeness, NA_real_),
  contamination = ifelse("Contamination" %in% names(tsv), tsv$Contamination, NA_real_),
  coding_density = ifelse("Coding_Density" %in% names(tsv), tsv$Coding_Density, NA_real_),
  total_coding_sequences = ifelse("Total_Coding_Sequences" %in% names(tsv), tsv$Total_Coding_Sequences, NA_real_),
  average_gene_length = ifelse("Average_Gene_Length" %in% names(tsv), tsv$Average_Gene_Length, NA_real_),
  max_contig_length = ifelse("Max_Contig_Length" %in% names(tsv), tsv$Max_Contig_Length, NA_real_),
  number_of_contigs = ifelse("Total_Contigs" %in% names(tsv), tsv$Total_Contigs, NA_real_),
  source_file_checkm2 = ifelse("source_file_checkm2" %in% names(tsv), tsv$source_file_checkm2, NA_character_)
)
colnames(checkm2) <- paste0(colnames(checkm2), "_checkm2")
checkm2$software_version <- paste0("StrainCascade v", opt$version)

# Save as qs
output_file <- file.path(opt$output_dir, "checkm2_results.qs")
qsave(checkm2, output_file, preset = "archive")

cat("CheckM2 results processed and saved to:", output_file, "\n")