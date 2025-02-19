#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_islandpath.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(optparse)
  library(Biostrings)
})

# Parse command-line arguments
option_list <- list(
  make_option("--output_dir", type="character", help="Output directory for qs files"),
  make_option("--gff_file", type="character", help="Path to islandpath.gff3 file"),
  make_option("--fasta", type="character", help="Path to FASTA file"),
  make_option("--version", type="character", help="StrainCascade version")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to import islandpath data
import_islandpath_data <- function(gff_file) {
  tryCatch({
    if (!file.exists(gff_file)) {
      stop(paste("GFF file does not exist:", gff_file))
    }
    islandpath_data <- readGFF(gff_file, version=3, columns=NULL, tags=NULL, filter=NULL, nrows=-1, raw_data=FALSE)
    if (nrow(islandpath_data) == 0) {
      warning("Empty GFF file. Returning empty data frame.")
      return(data.frame())
    }
    return(islandpath_data)
  }, error = function(e) {
    warning(paste("Error reading GFF file:", e$message, "Returning empty data frame."))
    return(data.frame())
  })
}

# Function to import FASTA file
import_fasta_file <- function(fasta_file) {
  tryCatch({
    if (!file.exists(fasta_file)) {
      stop(paste("FASTA file does not exist:", fasta_file))
    }
    fasta_data <- readDNAStringSet(fasta_file)
    if (length(fasta_data) == 0) {
      warning("Empty FASTA file.")
    }
    return(fasta_data)
  }, error = function(e) {
    stop(paste("Error reading FASTA file:", e$message))
  })
}

import_fasta_file_as_data_frame <- function(fasta_file) {
  fasta_data <- import_fasta_file(fasta_file)
  return(as.data.frame(fasta_data))
}

# Function to process islandpath data
process_islandpath_data <- function(islandpath_data, fasta_data, gff_file, version) {
  if (nrow(islandpath_data) == 0) {
    warning("No data in islandpath_data. Returning empty data frame.")
    return(data.frame(
      GI_id_islandpath = character(),
      contig_tag_islandpath = character(),
      start_position = numeric(),
      end_position = numeric(),
      length_bp = numeric(),
      strand = character(),
      type_islandpath = character(),
      phase_islandpath = character(),
      score_islandpath = numeric(),
      source_file_islandpath = character(),
      software_version = character(),
      nucleotide_code = character(),
      contig_tag_SC = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  islandpath_processed <- data.frame(
    GI_id_islandpath = islandpath_data$ID,
    contig_tag_islandpath = islandpath_data$seqid,  
    start_position = as.numeric(islandpath_data$start),
    end_position = as.numeric(islandpath_data$end),
    length_bp = as.numeric(islandpath_data$end) - as.numeric(islandpath_data$start) + 1,
    strand = islandpath_data$strand,
    type_islandpath = islandpath_data$type,
    phase_islandpath = islandpath_data$phase,
    score_islandpath = islandpath_data$score,
    source_file_islandpath = basename(gff_file),
    software_version = paste("StrainCascade v", version, sep = ""),
    stringsAsFactors = FALSE
  )
  
  # Function to extract nucleotide sequence
  extract_nucleotide_sequence <- function(contig, start, end, strand) {
    tryCatch({
      contig <- as.character(contig)
      names(fasta_data) <- as.character(names(fasta_data))
      
      if (contig %in% names(fasta_data)) {
        seq <- as.character(subseq(fasta_data[[contig]], start, end))
        if (strand == "-") {
          seq <- as.character(reverseComplement(DNAString(seq)))
        }
        return(seq)
      } else {
        warning(paste("Contig", contig, "not found in FASTA file"))
        return(NA)
      }
    }, error = function(e) {
      warning(paste("Error extracting nucleotide sequence:", e$message))
      return(NA)
    })
  }
  
  # Add nucleotide_code column
  islandpath_processed <- islandpath_processed %>%
    rowwise() %>%
    mutate(nucleotide_code = extract_nucleotide_sequence(contig_tag_islandpath, start_position, end_position, strand)) %>%
    ungroup()
  
  islandpath_processed$contig_tag_SC <- paste0("SC_", islandpath_processed$contig_tag_islandpath)
  
  return(islandpath_processed)
}

# Main function
islandpath_main <- function(gff_file, fasta_file, version) {
  cat("Importing islandpath data...\n")
  islandpath_data <- import_islandpath_data(gff_file)
  
  cat("Importing FASTA data...\n")
  fasta_data <- import_fasta_file(fasta_file)
  
  cat("Processing islandpath data...\n")
  processed_data <- process_islandpath_data(islandpath_data, fasta_data, gff_file, version)
  
  cat("Importing FASTA data as data frame...\n")
  fasta_frame <- import_fasta_file_as_data_frame(fasta_file)
  
  # Add contig_tag_fasta as a new column 
  fasta_frame$contig_tag_fasta <- rownames(fasta_frame)
  
  # Generate contig_tag_islandpath
  fasta_frame <- fasta_frame %>%
    mutate(
      contig_tag_islandpath = paste("contig_", row_number(), sep = ""),
      rownames = NULL
    ) 
  
  contig_tags <- fasta_frame %>%
    select(contig_tag_islandpath, contig_tag_fasta)
  
  processed_data <- processed_data %>%
    left_join(contig_tags, by = "contig_tag_islandpath")
  
  return(processed_data)
}

# Main execution
tryCatch({
  cat("Starting islandpath processing...\n")
  islandpath_results <- islandpath_main(opt$gff_file, opt$fasta, opt$version)
  
  # Define the desired order of columns for the islandpath_results dataframe
  new_column_order <- c(
    "GI_id_islandpath", 
    "contig_tag_fasta",
    "contig_tag_islandpath",
    "contig_tag_SC", 
    "type_islandpath", 
    "start_position", 
    "end_position",
    "length_bp",
    "strand",
    "nucleotide_code",
    "phase_islandpath", 
    "score_islandpath", 
    "source_file_islandpath",
    "software_version"
  )
  
  # Reorder the columns in the dataframe
  islandpath_results <- islandpath_results[, intersect(new_column_order, names(islandpath_results))]
  
  # If you want to save the results
  if (!is.null(islandpath_results) && nrow(islandpath_results) > 0) {
    output_file <- file.path(opt$output_dir, "islandpath_results.qs")
    qsave(islandpath_results, output_file, preset = "archive")
    cat(sprintf("SUCCESS: Results saved to %s\n", output_file))
  } else {
    warning("No results to save. Check your input files.")
  }
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})