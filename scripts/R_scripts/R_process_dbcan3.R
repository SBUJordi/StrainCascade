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
  make_option("--gff_file", type="character", default=NULL, help="Path to prodigal GFF file (provides positions for all predicted proteins)"),
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

# Function to import prodigal GFF file (provides positions for all predicted proteins)
import_gff_data <- function(gff_file) {
  gff_lines <- readLines(gff_file)
  gff_data_lines <- gff_lines[!grepl("^#", gff_lines)]

  if (length(gff_data_lines) == 0) {
    return(data.frame(
      locus_tag_dbcan3 = character(0),
      contig_tag_fasta = character(0),
      start_position = numeric(0),
      end_position = numeric(0),
      strand = character(0),
      stringsAsFactors = FALSE
    ))
  }

  gff_parsed <- read.delim(
    textConnection(gff_data_lines),
    header = FALSE, sep = "\t",
    col.names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"),
    stringsAsFactors = FALSE
  )

  # Extract protein ID from attributes (ID=contig_1_1;...)
  gff_parsed$locus_tag_dbcan3 <- sub("ID=([^;]+).*", "\\1", gff_parsed$attributes)

  gff_positions <- gff_parsed %>%
    select(
      locus_tag_dbcan3,
      contig_tag_fasta = seqname,
      start_position = start,
      end_position = end,
      strand
    ) %>%
    mutate(
      contig_tag_fasta = as.character(contig_tag_fasta),
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position),
      strand = as.character(strand)
    )

  return(gff_positions)
}

# Function to process dbcan3 data
# Uses GFF for positional data (all proteins), .txt for individual CAZyme annotations,
# and .out for CGC cluster info (when clusters exist)
process_dbcan3_data <- function(dbcan3_data, out_file, txt_file, gff_positions, fasta_data) {

  # --- Process .txt data (individual CAZyme annotations) ---
  # These are the primary results: every protein with a CAZyme hit
  dbcan3_txt_processed <- dbcan3_data$txt %>%
    dplyr::rename(
      locus_tag_dbcan3 = Gene.ID,
      EC_number_dbcan3 = EC.,
      HMMER_dbcan3 = dbCAN_hmm,
      dbCAN_sub = dbCAN_sub,
      DIAMOND_dbcan3 = DIAMOND,
      number_of_tools_dbcan3 = X.ofTools
    ) %>%
    mutate(source_file_dbcan3_txt = basename(txt_file))

  dbcan3_txt_processed[dbcan3_txt_processed == "-"] <- NA
  dbcan3_txt_processed[dbcan3_txt_processed == "null"] <- NA

  # --- Process .out data (CGC cluster info, may be empty) ---
  has_cgc <- nrow(dbcan3_data$out) > 0

  cgc_info <- NULL
  cgc_positions <- NULL
  if (has_cgc) {
    dbcan3_out_processed <- dbcan3_data$out %>%
      dplyr::rename(
        CAZyme_gene_clusters_dbcan3 = CGC.,
        gene_type_dbcan3 = Gene.Type,
        contig_tag_fasta = Contig.ID,
        locus_tag_dbcan3 = Protein.ID,
        start_position = Gene.Start,
        end_position = Gene.Stop,
        strand = Gene.Strand,
        protein_family_dbcan3 = Gene.Annotation
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

    # CGC-specific columns (cluster membership info)
    cgc_info <- dbcan3_out_processed %>%
      select(locus_tag_dbcan3, CAZyme_gene_clusters_dbcan3, gene_type_dbcan3,
             protein_family_dbcan3, source_file_dbcan3_out)

    # Positions from .out as fallback when GFF is not available
    cgc_positions <- dbcan3_out_processed %>%
      select(locus_tag_dbcan3, contig_tag_fasta, start_position, end_position, strand) %>%
      distinct(locus_tag_dbcan3, .keep_all = TRUE)
  }

  # --- Determine positions for all relevant proteins ---
  relevant_ids <- unique(dbcan3_txt_processed$locus_tag_dbcan3)
  if (has_cgc) {
    relevant_ids <- unique(c(relevant_ids, cgc_info$locus_tag_dbcan3))
  }

  has_gff <- !is.null(gff_positions) && nrow(gff_positions) > 0

  if (has_gff) {
    # GFF available: positional data for ALL proteins
    positions <- gff_positions %>%
      filter(locus_tag_dbcan3 %in% relevant_ids)
  } else if (has_cgc) {
    # No GFF: positions only from .out (CGC members only)
    positions <- cgc_positions
  } else {
    # No positions available
    positions <- data.frame(
      locus_tag_dbcan3 = character(0),
      contig_tag_fasta = character(0),
      start_position = numeric(0),
      end_position = numeric(0),
      strand = character(0),
      stringsAsFactors = FALSE
    )
  }

  # --- Build combined results ---
  # Start from all relevant protein IDs, add positions, annotations, and CGC info
  dbcan3_results <- data.frame(locus_tag_dbcan3 = relevant_ids, stringsAsFactors = FALSE) %>%
    left_join(positions, by = "locus_tag_dbcan3") %>%
    left_join(dbcan3_txt_processed, by = "locus_tag_dbcan3")

  if (has_cgc) {
    dbcan3_results <- dbcan3_results %>%
      left_join(cgc_info, by = "locus_tag_dbcan3")
  } else {
    dbcan3_results <- dbcan3_results %>%
      mutate(
        CAZyme_gene_clusters_dbcan3 = NA_character_,
        gene_type_dbcan3 = NA_character_,
        protein_family_dbcan3 = NA_character_,
        source_file_dbcan3_out = NA_character_
      )
  }

  # Individual CAZymes (from .txt) without CGC assignment get gene_type "CAZyme"
  dbcan3_results <- dbcan3_results %>%
    mutate(gene_type_dbcan3 = if_else(
      !is.na(source_file_dbcan3_txt) & is.na(gene_type_dbcan3),
      "CAZyme",
      gene_type_dbcan3
    ))

  # Ensure character types for join-critical columns
  dbcan3_results <- dbcan3_results %>%
    mutate(
      contig_tag_fasta = as.character(contig_tag_fasta),
      strand = as.character(strand)
    )

  # --- Extract nucleotide sequences ---
  extract_nucleotide_sequence <- function(contig, start, end, strand) {
    tryCatch({
      contig <- as.character(contig)
      if (is.na(contig) || is.na(start) || is.na(end) || is.na(strand)) {
        return(NA_character_)
      }
      names(fasta_data) <- as.character(names(fasta_data))
      if (contig %in% names(fasta_data)) {
        seq <- as.character(subseq(fasta_data[[contig]], start, end))
        if (strand == "-") {
          seq <- as.character(reverseComplement(DNAString(seq)))
        }
        return(seq)
      }
      return(NA_character_)
    }, error = function(e) {
      warning(paste("Error extracting nucleotide sequence for contig", contig, ":", e$message))
      return(NA_character_)
    })
  }

  dbcan3_results <- dbcan3_results %>%
    rowwise() %>%
    mutate(nucleotide_code = extract_nucleotide_sequence(contig_tag_fasta, start_position, end_position, strand)) %>%
    ungroup() %>%
    mutate(nucleotide_code = as.character(unlist(nucleotide_code)))

  return(dbcan3_results)
}

# Main function
dbcan3_main <- function(out_file, txt_file, gff_file, fasta_file, version) {
  dbcan3_data <- import_dbcan3_data(out_file, txt_file)
  fasta_data <- import_fasta_file(fasta_file)

  # Import prodigal GFF if available (provides positions for all predicted proteins)
  gff_positions <- NULL
  if (!is.null(gff_file) && file.exists(gff_file)) {
    gff_positions <- import_gff_data(gff_file)
    cat(sprintf("INFO: Loaded prodigal GFF with %d gene predictions\n", nrow(gff_positions)))
  } else {
    cat("WARNING: Prodigal GFF file not available. Positional data limited to CGC members only.\n")
  }

  processed_data <- process_dbcan3_data(dbcan3_data, out_file, txt_file, gff_positions, fasta_data)

  # Keep entries that have annotation data (from .txt) or positional data (from .out/.gff)
  processed_data <- processed_data[
    !is.na(processed_data$source_file_dbcan3_txt) | !is.na(processed_data$start_position),
  ]

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

  # Add any missing columns as NA (handles edge cases)
  for (col in setdiff(new_column_order, names(processed_data))) {
    processed_data[[col]] <- NA
  }

  processed_data <- processed_data[, new_column_order]

  return(processed_data)
}

# Main execution
tryCatch({
  dbcan3_results <- dbcan3_main(opt$out_file, opt$txt_file, opt$gff_file, opt$fasta, opt$version)

  if (!is.null(dbcan3_results) && nrow(dbcan3_results) > 0) {
    output_file <- file.path(opt$output_dir, "dbcan3_results.qs")
    qsave(dbcan3_results, output_file, preset = "archive")
    cat(sprintf("SUCCESS: Results saved to %s (%d entries)\n", output_file, nrow(dbcan3_results)))
  } else {
    cat("WARNING: No CAZyme results to save\n")
  }
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})