#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_microbeannotator_KEGG_modules.R

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
  make_option(c("-t", "--tab"), type="character", default=NULL, 
              help="Input tab file path", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$tab) || is.null(opt$version)) {
  stop("All arguments (output_dir, tab, and version) must be supplied.")
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

# Read and process the tab file
tab <- read.delim(opt$tab, header = TRUE, sep = "\t", fill = TRUE)
names(tab) <- c("M_number_microbeannotator", "KEGG_module_microbeannotator", "KEGG_module_group_microbeannotator", "KEGG_module_completeness_microbeannotator")

tab <- replace_na_strings(tab)

tab$source_file_microbeannotator <- basename(opt$tab)
tab$software_version <- paste0("StrainCascade v", opt$version)

# Save as qs
library(qs)
output_file <- file.path(opt$output_dir, "microbeannotator_KEGG_modules_results.qs")
qsave(tab, output_file, preset = "archive")

cat("MicrobeAnnotator KEGG module results processed and saved to:", output_file, "\n")