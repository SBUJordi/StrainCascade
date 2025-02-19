#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_quast.R

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

# Read and process the TSV file
tsv <- read.delim(opt$tsv, header = F, sep = "\t", fill = TRUE)
tsv <- replace_na_strings(tsv)
tsv$source_file_checkm2 <- basename(opt$tsv)

# Process column names
tsv <- tsv %>%
  mutate(across(everything(), ~ ifelse(row_number() == 1, gsub("[()]", "", .), .))) %>%
  mutate(across(everything(), ~ ifelse(row_number() == 1, gsub("'", "", .), .))) %>%
  mutate(across(everything(), ~ ifelse(row_number() == 1, gsub(">= ", "_greater_or_equal_", .), .))) %>%
  mutate(across(everything(), ~ ifelse(row_number() == 1, gsub("^# ", "number_of_", .), .))) %>%
  mutate(across(everything(), ~ ifelse(row_number() == 1, gsub(" ", "", .), .)))

colnames(tsv) <- tsv[1, ]
tsv <- tsv[-1, ]
colnames(tsv) <- paste0(colnames(tsv), "_quast")
note_regarding_file_naming <- tsv[2,1]
tsv <- tsv[1,]

# Create the quast dataframe
quast <- data.frame(
  assessed_file_quast = tsv$Assembly_quast,
  number_of_contigs_quast = as.numeric(tsv$number_of_contigs_quast),  
  largest_contig_quast = as.numeric(tsv$Largestcontig_quast),  
  total_length_quast = as.numeric(tsv$Totallength_quast),  
  number_of_contigs_greater_or_equal_1000bp_quast = as.numeric(tsv$number_of_contigs_greater_or_equal_1000bp_quast),
  number_of_contigs_greater_or_equal_5000bp_quast = as.numeric(tsv$number_of_contigs_greater_or_equal_5000bp_quast),
  number_of_contigs_greater_or_equal_10000bp_quast = as.numeric(tsv$number_of_contigs_greater_or_equal_10000bp_quast),
  number_of_contigs_greater_or_equal_25000bp_quast = as.numeric(tsv$number_of_contigs_greater_or_equal_25000bp_quast),
  number_of_contigs_greater_or_equal_50000bp_quast = as.numeric(tsv$number_of_contigs_greater_or_equal_50000bp_quast),
  total_length_greater_or_equal_0bp_quast = as.numeric(tsv$Totallength_greater_or_equal_0bp_quast),
  total_length_greater_or_equal_1000bp_quast = as.numeric(tsv$Totallength_greater_or_equal_1000bp_quast),
  total_length_greater_or_equal_5000bp_quast = as.numeric(tsv$Totallength_greater_or_equal_5000bp_quast),
  total_length_greater_or_equal_10000bp_quast = as.numeric(tsv$Totallength_greater_or_equal_10000bp_quast),
  total_length_greater_or_equal_25000bp_quast = as.numeric(tsv$Totallength_greater_or_equal_25000bp_quast),
  total_length_greater_or_equal_50000bp_quast = as.numeric(tsv$Totallength_greater_or_equal_50000bp_quast),
  GC_percentage_quast = as.numeric(tsv$`GC%_quast`),
  N50_quast = as.numeric(tsv$N50_quast),
  N90_quast = as.numeric(tsv$N90_quast),
  number_of_N_per100kbp_quast = as.numeric(tsv$number_of_Nsper100kbp_quast)
)

quast$software_version <- paste0("StrainCascade v", opt$version)
quast$source_file_quast <- basename(opt$tsv)
quast$note_regarding_file_naming <- note_regarding_file_naming

# Save as qs
output_file <- file.path(opt$output_dir, "quast_results.qs")
qsave(quast, output_file, preset = "archive")

cat("QUAST results processed and saved to:", output_file, "\n")