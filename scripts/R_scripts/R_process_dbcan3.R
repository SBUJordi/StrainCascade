#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_dbcan3.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(qs)
  library(optparse)
  library(Biostrings)
})

# Parse command-line arguments
option_list <- list(
  make_option("--output_dir", type="character", help="Output directory for qs files"),
  make_option("--out_file", type="character", help="Path to dbcan3.out file"),
  make_option("--txt_file", type="character", help="Path to dbcan3.txt file"),
  make_option("--fasta", type="character", help="Path to FASTA file"),
  make_option("--version", type="character", help="StrainCascade version")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to import dbcan3 data
import_dbcan3_data <- function(out_file, txt_file) {
  dbcan3_data_1 <- read.delim(out_file, header = TRUE, sep = "\t", fill = TRUE)
  dbcan3_data_2 <- read.delim(txt_file, header = TRUE, sep = "\t", fill = TRUE)
  return(list(out = dbcan3_data_1, txt = dbcan3_data_2))
}

# Function to import FASTA file
import_fasta_file <- function(fasta_file) {
  fasta_data <- readDNAStringSet(fasta_file)
  return(fasta_data)
}

# Function to process dbcan3 data
process_dbcan3_data <- function(dbcan3_data, out_file, txt_file, fasta_data) {

  # Process .out data
  dbcan3_out_processed <- dbcan3_data$out %>%
    dplyr::rename(
      CAZyme_gene_clusters_dbcan3 = CGC.,
      gene_type_dbcan3 = Gene.Type,
      contig_tag_fasta = Contig.ID,
      locus_tag_dbcan3 = Protein.ID,
      start_position = Gene.Start,
      end_position = Gene.Stop,
      strand = Direction,
      protein_family_dbcan3 = Protein.Family
    ) %>%
    mutate(
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position),
      source_file_dbcan3_out = basename(out_file)
    )

  dbcan3_out_processed[names(dbcan3_out_processed) != "strand"] <- 
    lapply(dbcan3_out_processed[names(dbcan3_out_processed) != "strand"], function(x) replace(x, x == "-", NA))
  dbcan3_out_processed[dbcan3_out_processed == "null"] <- NA

  dbcan3_out_processed$gene_type_dbcan3[dbcan3_out_processed$gene_type_dbcan3 == "TF"] <- "transcription factor"
  dbcan3_out_processed$gene_type_dbcan3[dbcan3_out_processed$gene_type_dbcan3 == "TC"] <- "transporter"
  dbcan3_out_processed$gene_type_dbcan3[dbcan3_out_processed$gene_type_dbcan3 == "STP"] <- "signal_transduction_proteins"

  # Process .txt data
  dbcan3_txt_processed <- dbcan3_data$txt %>%
    dplyr::rename(
      locus_tag_dbcan3 = Gene.ID,
      EC_number_dbcan3 = EC.,
      HMMER_dbcan3 = HMMER,
      dbCAN_sub = dbCAN_sub,
      DIAMOND_dbcan3 = DIAMOND,
      number_of_tools_dbcan3 = X.ofTools
    ) %>%
    mutate(
      source_file_dbcan3_txt = basename(txt_file)
    )

  dbcan3_txt_processed[dbcan3_txt_processed == "-"] <- NA
  dbcan3_txt_processed[dbcan3_txt_processed == "null"] <- NA

  # Merge processed data
  dbcan3_results <- merge(dbcan3_out_processed, dbcan3_txt_processed, by = "locus_tag_dbcan3", all = TRUE)

  # Function to complement bases
  complement_bases <- function(sequence) {
    chartr("ATCG", "TAGC", sequence)
  }

  # Function to extract nucleotide sequence
  extract_nucleotide_sequence <- function(contig, start, end, strand) {
    contig <- as.character(contig)  # Ensure character type because LJA creates numeric contig tags
    names(fasta_data) <- as.character(names(fasta_data))
    
    if (contig %in% names(fasta_data)) {
      seq <- as.character(subseq(fasta_data[[contig]], start, end))
      if (strand == "-") {
        seq <- as.character(reverseComplement(DNAString(seq)))  # Efficient reverse complement
      }
      return(seq)
    }
    return(NA)
  }

  # Add nucleotide_code column
  dbcan3_results <- dbcan3_results %>%
    rowwise() %>%
    mutate(nucleotide_code = extract_nucleotide_sequence(contig_tag_fasta, start_position, end_position, strand)) %>%
    ungroup()

  return(dbcan3_results)
}

# Main function
dbcan3_main <- function(out_file, txt_file, fasta_file, version) {
  dbcan3_data <- import_dbcan3_data(out_file, txt_file)
  fasta_data <- import_fasta_file(fasta_file)
  processed_data <- process_dbcan3_data(dbcan3_data, out_file, txt_file, fasta_data)

  processed_data <- processed_data[!is.na(processed_data$start_position),]

  processed_data$software_version <- paste("StrainCascade v", version, sep = "")

  # Define the desired order of columns for the dbcan3_results dataframe
  new_column_order <- c(
    "locus_tag_dbcan3",
    "gene_type_dbcan3",
    "CAZyme_gene_clusters_dbcan3",
    "EC_number_dbcan3",
    "protein_family_dbcan3",
    "contig_tag_fasta",
    "start_position",
    "end_position",
    "strand",
    "nucleotide_code",
    "HMMER_dbcan3",
    "DIAMOND_dbcan3",
    "dbCAN_sub",
    "number_of_tools_dbcan3",
    "source_file_dbcan3_out",
    "source_file_dbcan3_txt",
    "software_version"
  )

  # Reorder the columns in the dataframe
  processed_data <- processed_data[, new_column_order]

  return(processed_data)
}

# Main execution
tryCatch({
  dbcan3_results <- dbcan3_main(opt$out_file, opt$txt_file, opt$fasta, opt$version)

  # If you want to save the results
  if (!is.null(dbcan3_results)) {
    output_file <- file.path(opt$output_dir, "dbcan3_results.qs")
    qsave(dbcan3_results, output_file, preset = "archive")
    cat(sprintf("SUCCESS: Results saved to %s\n", output_file))
  }
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})