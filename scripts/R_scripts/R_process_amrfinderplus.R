#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_amrfinderplus.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option("--output_dir", type="character", help="Output directory for qs files"),
  make_option("--tsv_file1", type="character", help="Path to identification_amrfinderplus.tsv file"),
  make_option("--tsv_file2", type="character", help="Path to mutations_amrfinderplus.tsv file"),
  make_option("--fasta", type="character", help="Path to FASTA file"),
  make_option("--version", type="character", help="StrainCascade version")
)

opt <- parse_args(OptionParser(option_list=option_list))

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

# Function to check file status
check_file_status <- function(file_path, file_type) {
  if (!file.exists(file_path)) {
    cat(sprintf("INFO: %s file not found at %s\n", file_type, file_path))
    return("missing")
  } else if (file.size(file_path) == 0) {
    cat(sprintf("INFO: %s file is empty at %s\n", file_type, file_path))
    return("empty")
  }
  return("valid")
}

# Function to process AMRFinderPlus identification data
process_amrfinderplus_data <- function(tsv_file1, tsv_file2, fasta_file, version) {
  # Initialize empty dataframe with required structure
  empty_df <- data.frame(
    protein_id = character(),
    contig_tag_fasta = character(),
    start_position = integer(),
    end_position = integer(),
    length_bp = integer(),
    strand = character(),
    gene_symbol = character(),
    sequence_name = character(),
    scope = character(),
    element_type = character(),
    element_subtype = character(),
    resistance_class = character(),
    resistance_subclass = character(),
    method = character(),
    target_length = integer(),
    reference_length = integer(),
    coverage_percentage = numeric(),
    identity_percentage = numeric(),
    alignment_length = integer(),
    accession = character(),
    reference_name = character(),
    hmm_id = logical(),
    hmm_description = logical(),
    nucleotide_code = character(),
    mutation_details = character(),
    source_file = character(),
    software_version = character(),
    stringsAsFactors = FALSE
  )

  # Check identification file status
  id_file_status <- check_file_status(tsv_file1, "Identification")
  if (id_file_status == "missing" || id_file_status == "empty") {
    cat("WARNING: No valid identification data. Returning empty dataframe.\n")
    return(empty_df)
  }

  # Import FASTA data
  fasta_data <- import_fasta_file(fasta_file)

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

  # Read and process identification file
  tryCatch({
    df1 <- read.delim(tsv_file1, stringsAsFactors = FALSE)
    
    processed_df <- df1 %>%
      transmute(
        protein_id = Protein.id,
        contig_tag_fasta = Contig.id,
        start_position = Start,
        end_position = Stop,
        length_bp = Stop - Start + 1,
        strand = Strand,
        gene_symbol = Element.symbol,
        sequence_name = Element.name,
        scope = Scope,
        element_type = Type,
        element_subtype = Subtype,
        resistance_class = Class,
        resistance_subclass = Subclass,
        method = Method,
        target_length = Target.length,
        reference_length = Reference.sequence.length,
        coverage_percentage = X..Coverage.of.reference,
        identity_percentage = X..Identity.to.reference,
        alignment_length = Alignment.length,
        accession = Closest.reference.accession,
        reference_name = Closest.reference.name,
        hmm_id = HMM.accession,
        hmm_description = HMM.description,
        mutation_details = NA_character_,  # Initialize with NA
        source_file = basename(tsv_file1),
        software_version = paste0("StrainCascade v", version)
      )

    # Add nucleotide sequences
    processed_df <- processed_df %>%
      rowwise() %>%
      mutate(nucleotide_code = extract_nucleotide_sequence(
        contig_tag_fasta, 
        start_position, 
        end_position, 
        strand
      )) %>%
      ungroup()

    # Check mutations file status
    mutations_file_status <- check_file_status(tsv_file2, "Mutations")
    
    if (mutations_file_status == "valid") {
      tryCatch({
        mutations_df <- read.delim(tsv_file2, stringsAsFactors = FALSE)
        if (nrow(mutations_df) > 0) {
          cat(sprintf("INFO: Processing %d mutations from mutations file\n", nrow(mutations_df)))
          # Aggregate mutations by position
          mutations_summary <- mutations_df %>%
            group_by(Position) %>%
            summarize(
              mutation_details = paste(Mutation, collapse = "; "),
              .groups = 'drop'
            )
          
          # Match mutations with genes based on position
          for (i in seq_len(nrow(mutations_summary))) {
            pos <- mutations_summary$Position[i]
            matching_rows <- which(processed_df$start_position <= pos & 
                                 processed_df$end_position >= pos)
            if (length(matching_rows) > 0) {
              processed_df$mutation_details[matching_rows] <- 
                mutations_summary$mutation_details[i]
            }
          }
        } else {
          cat("INFO: Mutations file is empty (no rows)\n")
        }
      }, error = function(e) {
        cat(sprintf("WARNING: Error processing mutations file: %s\n", e$message))
      })
    } else {
      cat("INFO: No valid mutations file available. Mutation details will be NA.\n")
    }

    return(processed_df)
    
  }, error = function(e) {
    cat(sprintf("ERROR: Failed to process identification file: %s\n", e$message))
    return(empty_df)
  })
}

# Main execution
tryCatch({
  cat("Starting AMRFinderPlus data processing...\n")
  results <- process_amrfinderplus_data(opt$tsv_file1, opt$tsv_file2, opt$fasta, opt$version)
  
  # Define column order
  column_order <- c(
    "protein_id",
    "contig_tag_fasta", 
    "start_position", 
    "end_position",
    "length_bp",
    "strand",
    "gene_symbol",
    "sequence_name",
    "scope",
    "element_type",
    "element_subtype",
    "resistance_class",
    "resistance_subclass",
    "method",
    "target_length",
    "reference_length",
    "coverage_percentage",
    "identity_percentage",
    "alignment_length",
    "accession",
    "reference_name",
    "hmm_id",
    "hmm_description",
    "nucleotide_code",
    "mutation_details",
    "source_file",
    "software_version"
  )
  
  # Reorder columns
  results <- results[, intersect(column_order, names(results))]
  
  # Add information about data source
  cat(sprintf("INFO: Processed %d rows of AMR data\n", nrow(results)))
  
  # Save results
  output_file <- file.path(opt$output_dir, "amrfinderplus_results.qs")
  qsave(results, output_file, preset = "archive")
  cat(sprintf("SUCCESS: Results saved to %s\n", output_file))
  
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})