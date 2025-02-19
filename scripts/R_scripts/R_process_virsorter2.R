#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_virsorter2.R

set.seed(42)

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(optparse)
  library(Biostrings)
})

# Define command line options
option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="CHARACTER"),
  make_option(c("-f", "--combined_fa"), type="character", default=NULL, 
              help="Input combined.fa file path", metavar="CHARACTER"),
  make_option(c("-s", "--score_tsv"), type="character", default=NULL, 
              help="Input score.tsv file path", metavar="CHARACTER"),
  make_option(c("-b", "--boundary_tsv"), type="character", default=NULL, 
              help="Input boundary.tsv file path", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$combined_fa) || is.null(opt$score_tsv) || is.null(opt$boundary_tsv) || is.null(opt$version)) {
  stop("All arguments (output_dir, combined_fa, score_tsv, boundary_tsv, and version) must be supplied.")
}

# Function to safely read files
safe_read <- function(file_path, read_func, ...) {
  if (!file.exists(file_path)) {
    warning(paste("File does not exist:", file_path))
    return(data.frame())
  }
  tryCatch({
    df <- read_func(file_path, ...)
    if (nrow(df) == 0) {
      warning(paste("File is empty or contains only headers:", file_path))
    }
    return(df)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(data.frame())
  })
}

# Function to create empty data frame with specified columns
create_empty_df <- function(columns) {
  df <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(df) <- columns
  return(df)
}

# Read and process the .fa file
fa_file <- tryCatch({
  readDNAStringSet(opt$combined_fa)
}, error = function(e) {
  warning(paste("Error reading combined.fa file:", e$message))
  DNAStringSet()
})

if (length(fa_file) > 0) {
  fa_file <- data.frame(
    contig_tag_virsorter2 = names(fa_file),
    nucleotide_code = as.character(fa_file),
    length_bp = width(fa_file),
    source_file_virsorter2_fa = basename(opt$combined_fa),
    software_version = paste0("StrainCascade v", opt$version)
  )
} else {
  fa_file <- create_empty_df(c("contig_tag_virsorter2", "nucleotide_code", "length_bp", "source_file_virsorter2_fa", "software_version"))
}

# Read and process the score_tsv file
score <- safe_read(opt$score_tsv, read.delim, header = TRUE, sep = "\t", fill = TRUE)
if (nrow(score) > 0) {
  score <- score %>%
    dplyr::rename(
      contig_tag_virsorter2 = seqname,
      length_bp = length,
      score_virsorter2 = max_score,
      max_score_group_virsorter2 = max_score_group,
      hallmark_gene_count = hallmark
    ) %>%
    mutate(
      source_file_virsorter2_tsv1 = basename(opt$score_tsv),
      software_version = paste0("StrainCascade v", opt$version)
    )
} else {
  score <- create_empty_df(c("contig_tag_virsorter2", "length_bp", "score_virsorter2", "max_score_group_virsorter2", "hallmark_gene_count", "source_file_virsorter2_tsv1", "software_version"))
}

# Read and process the boundary_tsv file
boundary <- safe_read(opt$boundary_tsv, read.delim, header = TRUE, sep = "\t", fill = TRUE)
if (nrow(boundary) > 0) {
  boundary <- boundary %>%
    select(seqname_new, trim_bp_start, trim_bp_end, arc, bac, euk, vir, mix, unaligned) %>%
    dplyr::rename(
      contig_tag_virsorter2 = seqname_new,
      start_position = trim_bp_start,
      end_position = trim_bp_end,
      archeal = arc,
      bacterial = bac,
      eukaryotic = euk,
      viral = vir,
      mixed = mix,
      unaligned = unaligned
    ) %>%
    dplyr::mutate(
      length_bp = end_position - start_position + 1,
      source_file_virsorter2_tsv2 = basename(opt$boundary_tsv),
      software_version = paste0("StrainCascade v", opt$version)
    )
} else {
  boundary <- create_empty_df(c("contig_tag_virsorter2", "start_position", "end_position", "archeal", "bacterial", "eukaryotic", "viral", "mixed", "unaligned", "length_bp", "source_file_virsorter2_tsv2", "software_version"))
}

# Define the desired column order
desired_order <- c(
  "contig_tag_fasta",
  "contig_tag_virsorter2",
  "start_position",
  "end_position",
  "length_bp",
  "nucleotide_code",
  "dsDNAphage",
  "ssDNA",
  "score_virsorter2",
  "max_score_group_virsorter2",
  "hallmark_gene_count", 
  "viral",
  "cellular", 
  "archeal",
  "bacterial",
  "eukaryotic",
  "mixed",
  "unaligned",
  "software_version",
  "source_file_virsorter2_tsv1",
  "source_file_virsorter2_tsv2",
  "source_file_virsorter2_fa"
)

# Create an empty data frame with the desired structure
virsorter2_results <- create_empty_df(desired_order)

# Merge the data frames if they are not empty
if (nrow(fa_file) > 0 || nrow(score) > 0 || nrow(boundary) > 0) {
  virsorter2_results <- full_join(score, boundary, by = c("contig_tag_virsorter2", "length_bp", "software_version"))
  virsorter2_results <- full_join(virsorter2_results, fa_file, by = c("contig_tag_virsorter2", "length_bp", "software_version"))
  
  # Add contig_tag_fasta column
  virsorter2_results$contig_tag_fasta <- sub("\\|\\|.*", "", virsorter2_results$contig_tag_virsorter2)
  
  # Ensure all desired columns exist, add missing ones as NA
  for (col in desired_order) {
    if (!(col %in% names(virsorter2_results))) {
      virsorter2_results[[col]] <- NA
    }
  }
  
  # Reorder the columns
  virsorter2_results <- virsorter2_results[, desired_order]
}

# Save as qs
output_file <- file.path(opt$output_dir, "virsorter2_results.qs")
qsave(virsorter2_results, output_file, preset = "archive")

cat("VirSorter2 results processed and saved to:", output_file, "\n")
cat("Number of rows in final virsorter2_results:", nrow(virsorter2_results), "\n")