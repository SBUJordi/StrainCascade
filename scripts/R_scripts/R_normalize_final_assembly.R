#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_normalize_final_assembly.R
# Description: Normalizes contig names in the final assembly for reproducibility and consistency.
#              Sorts contigs by length (descending), renames to contig_1, contig_2, etc.,
#              updates the FASTA file, and creates a mapping file for traceability.

set.seed(42)

suppressPackageStartupMessages({
  library(Biostrings)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-f", "--fasta"), type = "character", default = NULL,
              help = "Path to the final assembly FASTA file", metavar = "CHARACTER"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory for mapping file", metavar = "CHARACTER"),
  make_option(c("-c", "--circlator_qs"), type = "character", default = NULL,
              help = "Path to Circlator results QS file (optional)", metavar = "CHARACTER")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$fasta) || is.null(opt$output_dir)) {
  stop("Both --fasta and --output_dir arguments are required.")
}

if (!file.exists(opt$fasta)) {
  stop(sprintf("FASTA file not found: %s", opt$fasta))
}

tryCatch({
  # Read the assembly
  assembly <- readDNAStringSet(opt$fasta)
  
  if (length(assembly) == 0) {
    stop("Assembly file is empty.")
  }
  
  # Get original names and lengths
  original_names <- names(assembly)
  contig_lengths <- width(assembly)
  
  # Sort by length (descending) for reproducible ordering
  size_order <- order(contig_lengths, decreasing = TRUE)
  assembly <- assembly[size_order]
  original_names_sorted <- original_names[size_order]
  contig_lengths_sorted <- contig_lengths[size_order]
  
  # Create new standardized names
  normalized_names <- paste0("contig_", seq_along(assembly))
  names(assembly) <- normalized_names
  
  # Create mapping dataframe
  contig_mapping <- data.frame(
    original_name = original_names_sorted,
    normalized_name = normalized_names,
    contig_length_bp = contig_lengths_sorted,
    stringsAsFactors = FALSE
  )
  
  # Also create a column with stripped polishing suffixes for matching with pre-polishing tools
  strip_polishing_suffix <- function(name) {
    name <- gsub(" (polypolish|medaka|racon|arrow)$", "", name)
    name <- gsub("_(polypolish|medaka|racon|arrow)$", "", name)
    return(name)
  }
  contig_mapping$original_name_stripped <- sapply(contig_mapping$original_name, strip_polishing_suffix)
  
  # Write the normalized FASTA (overwrite original)
  writeXStringSet(assembly, opt$fasta)
  cat(sprintf("SUCCESS: Normalized FASTA written to %s\n", opt$fasta))
  cat(sprintf("INFO: %d contigs sorted by length and renamed (contig_1 to contig_%d)\n", 
              length(assembly), length(assembly)))
  
  # Write the mapping file
  mapping_file <- file.path(opt$output_dir, "contig_name_mapping.tsv")
  write.table(contig_mapping, mapping_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("SUCCESS: Contig name mapping saved to %s\n", mapping_file))
  
  # Update Circlator results if provided
  if (!is.null(opt$circlator_qs) && file.exists(opt$circlator_qs)) {
    library(qs)
    
    circlator_results <- qread(opt$circlator_qs)
    
    if ("contig_tag_fasta" %in% colnames(circlator_results)) {
      # Match Circlator contig names to mapping (Circlator has pre-polish names)
      # First try exact match, then try stripped match
      circlator_results$contig_tag_original <- circlator_results$contig_tag_fasta
      
      # Create lookup from stripped names to normalized names
      stripped_to_normalized <- setNames(
        contig_mapping$normalized_name,
        contig_mapping$original_name_stripped
      )
      
      # Also create lookup from original names (in case no polishing happened)
      original_to_normalized <- setNames(
        contig_mapping$normalized_name,
        contig_mapping$original_name
      )
      
      # Try to match - first exact original, then stripped
      circlator_results$contig_tag_fasta <- sapply(circlator_results$contig_tag_fasta, function(name) {
        if (name %in% names(original_to_normalized)) {
          return(original_to_normalized[name])
        } else if (name %in% names(stripped_to_normalized)) {
          return(stripped_to_normalized[name])
        } else {
          # No match found, keep original
          return(name)
        }
      })
      
      # Save updated Circlator results
      qsave(circlator_results, opt$circlator_qs, preset = "archive")
      cat(sprintf("SUCCESS: Circlator results updated with normalized contig names\n"))
    }
  }
  
  cat("SUCCESS: Assembly normalization completed\n")
  
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  quit(save = "no", status = 1)
})
