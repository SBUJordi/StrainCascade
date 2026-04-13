#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_isescan.R

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(optparse)
  library(Biostrings)
})

option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="CHARACTER"),
  make_option(c("-t", "--tsv"), type="character", default=NULL, 
              help="Input TSV file path", metavar="CHARACTER"),
  make_option(c("-i", "--is_fna"), type="character", default=NULL, 
              help="Input IS DNA sequences file path", metavar="CHARACTER"),
  make_option(c("-f", "--orf_fna"), type="character", default=NULL, 
              help="Input ORF DNA sequences file path", metavar="CHARACTER"),
  make_option(c("-a", "--orf_faa"), type="character", default=NULL, 
              help="Input ORF amino acid sequences file path", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to safely read files
safe_read <- function(file_path, read_function, ...) {
  if (is.null(file_path) || !file.exists(file_path) || file.size(file_path) == 0) {
    warning(paste("File not found or empty:", file_path))
    return(data.frame())
  }
  tryCatch({
    read_function(file_path, ...)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(data.frame())
  })
}

# Function to safely convert to integer (silently returns NA for non-numeric values)
safe_as_integer <- function(x) {
  suppressWarnings(as.integer(x))
}

# Read files safely
tsv_data <- safe_read(opt$tsv, read.delim, header = TRUE, sep = "\t", fill = TRUE)
is_sequences <- safe_read(opt$is_fna, function(x) as.data.frame(readDNAStringSet(x)))
orf_dna_sequences <- safe_read(opt$orf_fna, function(x) as.data.frame(readDNAStringSet(x)))
orf_aa_sequences <- safe_read(opt$orf_faa, function(x) as.data.frame(readAAStringSet(x)))

# Create an empty dataframe with the expected structure
empty_isescan_results <- data.frame(
  contig_tag_fasta = character(),
  insertion_sequence_family = character(),
  tpase_cluster = character(),
  insertion_sequence_start_position = integer(),
  insertion_sequence_end_position = integer(),
  strand = character(),
  insertion_sequence_length_bp = integer(),
  insertion_sequence_copy_number = integer(),
  first_inverted_repeat_start_position = integer(),
  first_inverted_repeat_end_position = integer(),
  second_inverted_repeat_start_position = integer(),
  second_inverted_repeat_end_position = integer(),
  inverted_repeat_score_isescan = numeric(),
  inverted_repeat_identity_matches_bp = integer(),
  inverted_repeat_length_bp = integer(),
  inverted_repeat_alignment_gaps_number = integer(),
  orf_start_position = integer(),
  orf_end_position = integer(),
  orf_length_bp = integer(),
  terminal_inverted_repeat_sequence = character(),
  insertion_sequence_element_completness = character(),
  insertion_sequence_E_value_hmmer = numeric(),
  insertion_sequence_best_E_value_hmmer_among_copies = numeric(),
  insertion_sequence_OV_value_hmmer = numeric(),
  orf_nucleotide_code = character(),
  orf_amino_acid_code = character(),
  insertion_sequence_nucleotide_code = character(),
  stringsAsFactors = FALSE
)

# Process TSV data if it exists
if (nrow(tsv_data) > 0) {
  isescan_results <- tsv_data %>%
    dplyr::rename(
      contig_tag_fasta = seqID,
      insertion_sequence_family = family,
      tpase_cluster = cluster,
      insertion_sequence_start_position = isBegin,
      insertion_sequence_end_position = isEnd,
      insertion_sequence_length_bp = isLen,
      strand = strand,
      insertion_sequence_copy_number = ncopy4is,
      first_inverted_repeat_start_position = start1,
      first_inverted_repeat_end_position = end1,
      second_inverted_repeat_start_position = start2,
      second_inverted_repeat_end_position = end2,
      inverted_repeat_score_isescan = score,
      inverted_repeat_identity_matches_bp = irId,
      inverted_repeat_length_bp = irLen,
      inverted_repeat_alignment_gaps_number = nGaps,
      orf_start_position = orfBegin,
      orf_end_position = orfEnd,
      orf_length_bp = orfLen,
      terminal_inverted_repeat_sequence = tir,
      insertion_sequence_element_completness = type,
      insertion_sequence_E_value_hmmer = E.value,
      insertion_sequence_best_E_value_hmmer_among_copies = E.value4copy,
      insertion_sequence_OV_value_hmmer = ov
    ) %>%
    mutate(
      insertion_sequence_element_completness = case_when(
        insertion_sequence_element_completness == "c" ~ "complete_insertion_sequence_element",
        insertion_sequence_element_completness == "p" ~ "partial_insertion_sequence_element",
        TRUE ~ insertion_sequence_element_completness
      ),
      terminal_inverted_repeat_sequence = ifelse(
        strand == "-",
        paste0(sub("^[^:]*:", "", terminal_inverted_repeat_sequence), ":", sub(":.*$", "", terminal_inverted_repeat_sequence)),
        terminal_inverted_repeat_sequence
      )
    )
} else {
  isescan_results <- empty_isescan_results
}

# Process sequence data if it exists
if (nrow(is_sequences) > 0) {
  is_sequences <- is_sequences %>%
    mutate(
      full_header = rownames(is_sequences),
      insertion_sequence_nucleotide_code = x,
      # Split header: "contig_1_314476_316624_- IS21_259" or "contig_1_684983_686492 ISL3_126|..."
      # First part (before space) contains contig and positions, second part is cluster
      header_parts = strsplit(full_header, " ", fixed = TRUE),
      contig_part = sapply(header_parts, `[`, 1),
      tpase_cluster = sapply(header_parts, function(x) if(length(x) > 1) paste(x[-1], collapse = " ") else NA_character_),
      # Check if strand is present at the end of contig_part (e.g., "_-" or "_+")
      has_strand = grepl("_[+-]$", contig_part),
      strand = ifelse(has_strand, sub(".*_([+-])$", "\\1", contig_part), NA_character_),
      # Remove strand suffix if present
      contig_part = ifelse(has_strand, sub("_[+-]$", "", contig_part), contig_part),
      # Now contig_part is like "contig_1_314476_316624"
      # Extract end position (last number after underscore)
      insertion_sequence_end_position = safe_as_integer(sub(".*_([0-9]+)$", "\\1", contig_part)),
      contig_part = sub("_[0-9]+$", "", contig_part),
      # Extract start position
      insertion_sequence_start_position = safe_as_integer(sub(".*_([0-9]+)$", "\\1", contig_part)),
      contig_tag_fasta = sub("_[0-9]+$", "", contig_part)
    ) %>%
    select(-x, -full_header, -header_parts, -contig_part, -has_strand)
  
  isescan_results <- merge(isescan_results, is_sequences, 
                           by = c("contig_tag_fasta", "tpase_cluster", "insertion_sequence_start_position", 
                                  "insertion_sequence_end_position", "strand"), 
                           all = TRUE)
}

# Process ORF DNA sequences if they exist
# ORF headers: "contig_1_314943_316304_-" (with strand) or "contig_1_314943_316304" (without)
if (nrow(orf_dna_sequences) > 0) {
  orf_dna_sequences <- orf_dna_sequences %>%
    mutate(
      contig_part = rownames(orf_dna_sequences),
      orf_nucleotide_code = x,
      # Check if strand is present at the end (e.g., "_-" or "_+")
      has_strand = grepl("_[+-]$", contig_part),
      strand = ifelse(has_strand, sub(".*_([+-])$", "\\1", contig_part), NA_character_),
      # Remove strand suffix if present
      contig_part = ifelse(has_strand, sub("_[+-]$", "", contig_part), contig_part),
      # Extract end position (last number after underscore)
      orf_end_position = safe_as_integer(sub(".*_([0-9]+)$", "\\1", contig_part)),
      contig_part = sub("_[0-9]+$", "", contig_part),
      # Extract start position
      orf_start_position = safe_as_integer(sub(".*_([0-9]+)$", "\\1", contig_part)),
      contig_tag_fasta = sub("_[0-9]+$", "", contig_part)
    ) %>%
    select(-x, -contig_part, -has_strand)
  
  isescan_results <- merge(isescan_results, orf_dna_sequences, 
                           by = c("contig_tag_fasta", "orf_start_position", 
                                  "orf_end_position", "strand"), 
                           all = TRUE)
}

# Process ORF amino acid sequences if they exist
# ORF headers: "contig_1_314943_316304_-" (with strand) or "contig_1_314943_316304" (without)
if (nrow(orf_aa_sequences) > 0) {
  orf_aa_sequences <- orf_aa_sequences %>%
    mutate(
      contig_part = rownames(orf_aa_sequences),
      orf_amino_acid_code = x,
      # Check if strand is present at the end (e.g., "_-" or "_+")
      has_strand = grepl("_[+-]$", contig_part),
      strand = ifelse(has_strand, sub(".*_([+-])$", "\\1", contig_part), NA_character_),
      # Remove strand suffix if present
      contig_part = ifelse(has_strand, sub("_[+-]$", "", contig_part), contig_part),
      # Extract end position (last number after underscore)
      orf_end_position = safe_as_integer(sub(".*_([0-9]+)$", "\\1", contig_part)),
      contig_part = sub("_[0-9]+$", "", contig_part),
      # Extract start position
      orf_start_position = safe_as_integer(sub(".*_([0-9]+)$", "\\1", contig_part)),
      contig_tag_fasta = sub("_[0-9]+$", "", contig_part)
    ) %>%
    select(-x, -contig_part, -has_strand)
  
  isescan_results <- merge(isescan_results, orf_aa_sequences, 
                           by = c("contig_tag_fasta", "orf_start_position", 
                                  "orf_end_position", "strand"), 
                           all = TRUE)
}

# Ensure all expected columns are present, even if empty
for (col in names(empty_isescan_results)) {
  if (!(col %in% names(isescan_results))) {
    isescan_results[[col]] <- NA
  }
}

# Add metadata (even if the dataframe is empty)
isescan_results$source_file_tsv <- ifelse(is.null(opt$tsv), NA, basename(opt$tsv))
isescan_results$source_file_is_fna <- ifelse(is.null(opt$is_fna), NA, basename(opt$is_fna))
isescan_results$source_file_orf_fna <- ifelse(is.null(opt$orf_fna), NA, basename(opt$orf_fna))
isescan_results$source_file_orf_faa <- ifelse(is.null(opt$orf_faa), NA, basename(opt$orf_faa))
isescan_results$software_version <- paste0("StrainCascade v", opt$version)

# Define the desired column order
desired_order <- c(
  "contig_tag_fasta",
  "tpase_cluster",
  "insertion_sequence_family",
  "insertion_sequence_start_position",
  "insertion_sequence_end_position",
  "strand",
  "insertion_sequence_length_bp",
  "insertion_sequence_copy_number",
  "first_inverted_repeat_start_position",
  "first_inverted_repeat_end_position",
  "second_inverted_repeat_start_position",
  "second_inverted_repeat_end_position",
  "inverted_repeat_score_isescan",
  "inverted_repeat_identity_matches_bp",
  "inverted_repeat_length_bp",
  "inverted_repeat_alignment_gaps_number",
  "orf_start_position",
  "orf_end_position",
  "orf_length_bp",
  "orf_nucleotide_code",
  "orf_amino_acid_code",
  "terminal_inverted_repeat_sequence",
  "insertion_sequence_nucleotide_code",
  "insertion_sequence_E_value_hmmer",
  "insertion_sequence_best_E_value_hmmer_among_copies",
  "insertion_sequence_element_completness",
  "insertion_sequence_OV_value_hmmer",
  "source_file_tsv",
  "source_file_is_fna",
  "source_file_orf_fna",
  "source_file_orf_faa",
  "software_version"
)

# Reorder the columns based on availability
isescan_results <- isescan_results[, intersect(desired_order, names(isescan_results))]

# Ensure at least one row exists
if (nrow(isescan_results) == 0) {
  isescan_results <- isescan_results[1, ]
}

# Save as qs
output_file <- file.path(opt$output_dir, "isescan_results.qs")
qsave(isescan_results, output_file, preset = "archive")

cat("ISEScan results processed and saved to:", output_file, "\n")