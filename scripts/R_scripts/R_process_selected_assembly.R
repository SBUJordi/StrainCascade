#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_selected_assembly.R

set.seed(42)

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(optparse)
  library(Biostrings)
  library(stringr)
})

# Define command line options
option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="CHARACTER"),
  make_option(c("-f", "--fasta"), type="character", default=NULL, 
              help="Input FASTA file path", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$fasta) || is.null(opt$version)) {
  stop("All arguments must be supplied.")
}

# Read the FASTA file
assembly <- readDNAStringSet(opt$fasta)

# Process the assembly
assembly_data <- data.frame(
  contig_tag_fasta = names(assembly),
  contig_length_bp = width(assembly),
  nucleotide_code = as.character(assembly)
) %>%
  mutate(
    A_percentage = str_count(nucleotide_code, "A") / contig_length_bp * 100,
    C_percentage = str_count(nucleotide_code, "C") / contig_length_bp * 100,
    G_percentage = str_count(nucleotide_code, "G") / contig_length_bp * 100,
    T_percentage = str_count(nucleotide_code, "T") / contig_length_bp * 100,
    N_percentage = str_count(nucleotide_code, "N") / contig_length_bp * 100
  )

# Add metadata
assembly_data$source_file_fasta <- basename(opt$fasta)
assembly_data$software_version <- paste0("StrainCascade v", opt$version)

# Save as qs
output_file <- file.path(opt$output_dir, "selected_assembly_results.qs")
qsave(assembly_data, output_file, preset = "archive")

cat("Selected assembly processed and saved to:", output_file, "\n")