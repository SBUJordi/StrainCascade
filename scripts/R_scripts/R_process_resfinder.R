#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_resfinder.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(Biostrings)
  library(qs)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="PATH"),
  make_option(c("--tab_file"), type="character", default=NULL, 
              help="Path to ResFinder_results_tab file", metavar="FILE"),
  make_option(c("--pheno_file"), type="character", default=NULL, 
              help="Path to pheno_table file", metavar="FILE"),
  make_option(c("--hit_fsa"), type="character", default=NULL, 
              help="Path to ResFinder_Hit_in_genome_seq file", metavar="FILE"),
  make_option(c("--ref_fsa"), type="character", default=NULL, 
              help="Path to ResFinder_Resistance_gene_seq file", metavar="FILE"),
  make_option(c("--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$output_dir) || is.null(opt$tab_file) || is.null(opt$pheno_file) || 
    is.null(opt$hit_fsa) || is.null(opt$ref_fsa) || is.null(opt$version)) {
  print_help(opt_parser)
  stop("All arguments must be supplied.", call.=FALSE)
}

# Function to import resfinder data
import_resfinder_data <- function(tab_file, pheno_file, hit_fsa, ref_fsa) {
  # Import .tsv file (handle empty dataframe)
  summary_data <- if(file.size(tab_file) > 0) {
    read.delim(tab_file, header = TRUE, sep = "\t", fill = TRUE)
  } else {
    data.frame()
  }
  
  # Import pheno_table file (unchanged, as it should never be empty)
  pheno_table_data <- read.delim(pheno_file, skip = 15, header = TRUE, sep = "\t", fill = TRUE)
  
  # Import .fsa file hit (handle empty dataframe)
  resfinder_hit_nucleotide_seq <- if(file.size(hit_fsa) > 0) {
    seq_data <- readDNAStringSet(hit_fsa)
    as.data.frame(seq_data) %>%
      dplyr::rename(nucleotide_code_resfinder_hit = x) %>%
      mutate(
        accession_number_resfinder = rownames(.),
        rownames = NULL
      ) %>%
      separate(col = 2, into = c("AMR_gene_resfinder", "nucleotide_identity_resfinder", "alignment_length_to_reference_resfinder", "percentage_length_of_reference_sequence_resfinder", "position_in_reference_resfinder", "contig_tag_fasta", "position_in_contig_resfinder"), sep = ", ", fill = "right") %>%
      mutate(across(everything(), ~ sub(".*: ", "", .))) %>%
      select(-nucleotide_identity_resfinder, -percentage_length_of_reference_sequence_resfinder)
  } else {
    data.frame(
      nucleotide_code_resfinder_hit = character(),
      accession_number_resfinder = character(),
      AMR_gene_resfinder = character(),
      alignment_length_to_reference_resfinder = character(),
      position_in_reference_resfinder = character(),
      contig_tag_fasta = character(),
      position_in_contig_resfinder = character()
    )
  }
  
  # Import .fsa file ref (handle empty dataframe)
  resfinder_reference_nucleotide_seq <- if(file.size(ref_fsa) > 0) {
    seq_data <- readDNAStringSet(ref_fsa)
    as.data.frame(seq_data) %>%
      dplyr::rename(nucleotide_code_resfinder_reference = x) %>%
      mutate(
        info = rownames(.),
        rownames = NULL
      ) %>%
      separate(col = info, into = c("AMR_gene_resfinder", "accession_number_resfinder"), sep = "_")
  } else {
    data.frame(
      nucleotide_code_resfinder_reference = character(),
      AMR_gene_resfinder = character(),
      accession_number_resfinder = character()
    )
  }
  
  return(list(summary = summary_data, pheno_table = pheno_table_data, hit_nucleotide_seq = resfinder_hit_nucleotide_seq, reference_nucleotide_seq = resfinder_reference_nucleotide_seq))
}

# Function to process resfinder data
process_resfinder_data <- function(resfinder_data, tab_file, pheno_file, hit_fsa, ref_fsa, version) {
  # Process results data (handle empty dataframe)
  results_resfinder_processed <- if(nrow(resfinder_data$summary) > 0) {
    resfinder_data$summary %>%
      dplyr::rename(
        AMR_gene_resfinder = Resistance.gene,
        nucleotide_identity_resfinder = Identity,
        alignment_length_to_reference_resfinder = Alignment.Length.Gene.Length,
        percentage_length_of_reference_sequence_resfinder = Coverage,
        position_in_reference_resfinder = Position.in.reference,
        contig_tag_fasta = Contig,
        position_in_contig_resfinder = Position.in.contig,
        AMR_phenotype_resfinder = Phenotype,
        accession_number_resfinder = Accession.no.,
      ) %>%
      mutate(
        source_file_results_resfinder = basename(tab_file),
        rownames = NULL
      ) %>%
      separate(position_in_contig_resfinder, into = c("start_position", "end_position"), sep = "\\.\\.") %>%
      separate(position_in_reference_resfinder, into = c("start_position_resfinder_reference", "end_position_resfinder_reference"), sep = "\\.\\.") %>%
      mutate(
        start_position = as.numeric(start_position),
        end_position = as.numeric(end_position),
        start_position_resfinder_reference = as.numeric(start_position_resfinder_reference),
        end_position_resfinder_reference = as.numeric(end_position_resfinder_reference),
        source_file_results_resfinder = basename(tab_file),
        rownames = NULL
      )
  } else {
    data.frame(
      AMR_gene_resfinder = character(),
      nucleotide_identity_resfinder = numeric(),
      alignment_length_to_reference_resfinder = character(),
      percentage_length_of_reference_sequence_resfinder = numeric(),
      contig_tag_fasta = character(),
      AMR_phenotype_resfinder = character(),
      accession_number_resfinder = character(),
      source_file_results_resfinder = character(),
      start_position = numeric(),
      end_position = numeric(),
      start_position_resfinder_reference = numeric(),
      end_position_resfinder_reference = numeric()
    )
  }
  
  # Process AMR_profile (phenotype) - unchanged as pheno_file should never be empty
  AMR_profile_resfinder_processed <- as.data.frame(resfinder_data$pheno_table) %>%
    dplyr::rename(
      antimicrobial_agent_resfinder = X..Antimicrobial,
      antimicrobial_agent_class_resfinder = Class,
      genomic_AMR_phenotype_resfinder = WGS.predicted.phenotype,
      AMR_phenotype_confidence_resfinder = Match,
      AMR_phenotype_genetic_background_resfinder = Genetic.background
    ) %>%
    mutate(
      source_file_AMR_profile_resfinder = basename(pheno_file),
      rownames = NULL
    )
  
  # Process nucleotide code of resfinder hit (handle empty dataframe)
  nucleotide_code_AMR_hit_resfinder_processed <- if(nrow(resfinder_data$hit_nucleotide_seq) > 0) {
    as.data.frame(resfinder_data$hit_nucleotide_seq) %>%
      mutate(
        source_file_nucleotide_code_AMR_hit_resfinder_processed = basename(hit_fsa),
        rownames = NULL
      ) %>%
      separate(position_in_contig_resfinder, into = c("start_position", "end_position"), sep = "\\.\\.") %>%
      separate(position_in_reference_resfinder, into = c("start_position_resfinder_reference", "end_position_resfinder_reference"), sep = "\\.\\.") %>%
      mutate(
        start_position = as.numeric(start_position),
        end_position = as.numeric(end_position),
        start_position_resfinder_reference = as.numeric(start_position_resfinder_reference),
        end_position_resfinder_reference = as.numeric(end_position_resfinder_reference)
      )
  } else {
    data.frame(
      nucleotide_code_resfinder_hit = character(),
      AMR_gene_resfinder = character(),
      contig_tag_fasta = character(),
      start_position = numeric(),
      end_position = numeric(),
      source_file_nucleotide_code_AMR_hit_resfinder_processed = character()
    )
  }
  
  # Process nucleotide code of resfinder reference (handle empty dataframe)
  nucleotide_code_AMR_reference_resfinder_processed <- if(nrow(resfinder_data$reference_nucleotide_seq) > 0) {
    as.data.frame(resfinder_data$reference_nucleotide_seq) %>%
      mutate(
        source_file_nucleotide_code_AMR_reference_resfinder_processed = basename(ref_fsa),
        rownames = NULL
      )
  } else {
    data.frame(
      nucleotide_code_resfinder_reference = character(),
      AMR_gene_resfinder = character(),
      accession_number_resfinder = character(),
      source_file_nucleotide_code_AMR_reference_resfinder_processed = character()
    )
  }
  
  return(list(AMR_results_resfinder = results_resfinder_processed, AMR_profile_resfinder = AMR_profile_resfinder_processed, AMR_gene_nucleotide_code_resfinder = nucleotide_code_AMR_hit_resfinder_processed, AMR_gene_reference_nucleotide_code_resfinder=nucleotide_code_AMR_reference_resfinder_processed))
}

# Main execution
tryCatch({
  resfinder_data <- import_resfinder_data(opt$tab_file, opt$pheno_file, opt$hit_fsa, opt$ref_fsa)
  processed_data <- process_resfinder_data(resfinder_data, opt$tab_file, opt$pheno_file, opt$hit_fsa, opt$ref_fsa, opt$version)
  
  AMR_profile_resfinder <- processed_data$AMR_profile_resfinder
  AMR_gene_nucleotide_code_resfinder <- processed_data$AMR_gene_nucleotide_code_resfinder
  AMR_gene_reference_nucleotide_code_resfinder <- processed_data$AMR_gene_reference_nucleotide_code_resfinder
  AMR_results_resfinder <- processed_data$AMR_results_resfinder
  
  # Add additional error handling here
  if(nrow(AMR_results_resfinder) == 0) {
    cat("WARNING: AMR_results_resfinder is empty. Skipping merges and saving empty dataframes.\n")
  } else {
    # Perform merges only if AMR_results_resfinder is not empty
    if(nrow(AMR_gene_nucleotide_code_resfinder) > 0) {
      AMR_results_resfinder <- merge(AMR_results_resfinder, AMR_gene_nucleotide_code_resfinder, 
                                     by = c("AMR_gene_resfinder", "start_position", "end_position", "contig_tag_fasta"), 
                                     all = TRUE)
    }
    
    if(nrow(AMR_gene_reference_nucleotide_code_resfinder) > 0) {
      AMR_results_resfinder <- merge(AMR_results_resfinder, AMR_gene_reference_nucleotide_code_resfinder, 
                                     by = c("AMR_gene_resfinder", "accession_number_resfinder"), 
                                     all = TRUE)
    }
    
    AMR_results_resfinder$software_version <- paste("StrainCascade v", opt$version, sep = "")
  }

  # Define the desired order of columns for the AMR_results_resfinder dataframe
  new_column_order <- c(
    # AMR gene and phenotype information
    "AMR_gene_resfinder",
    "AMR_phenotype_resfinder",
    
    # Sequence alignment and coverage information
    "start_position",
    "end_position",
    "start_position_resfinder_reference",
    "end_position_resfinder_reference",
    "nucleotide_identity_resfinder",
    "alignment_length_to_reference_resfinder",
    "percentage_length_of_reference_sequence_resfinder",
    
    # Accession and contig information
    "accession_number_resfinder",
    "contig_tag_fasta",
    
    # Nucleotide code and source file information
    "nucleotide_code_resfinder_hit",
    "nucleotide_code_resfinder_reference",
    "source_file_results_resfinder",
    "source_file_nucleotide_code_AMR_hit_resfinder_processed",
    "source_file_nucleotide_code_AMR_reference_resfinder_processed",
    
    # Software version information
    "software_version"
  )

  # Reorder columns if AMR_results_resfinder is not empty
  if(nrow(AMR_results_resfinder) > 0) {
    AMR_results_resfinder <- AMR_results_resfinder[, new_column_order[new_column_order %in% names(AMR_results_resfinder)]]
  }

  # Save the results
  qsave(AMR_results_resfinder, file.path(opt$output_dir, "resfinder_results.qs"), preset = "archive")
  qsave(AMR_profile_resfinder, file.path(opt$output_dir, "resfinder_AMR_profile.qs"), preset = "archive")
  cat("SUCCESS: Results saved to resfinder_results.qs and resfinder_AMR_profile.qs\n")
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})