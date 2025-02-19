#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_bakta.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(Biostrings)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(optparse)
})


# Parse command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character", default=NULL, help="Output directory for results"),
  make_option(c("--tsv"), type="character", default=NULL, help="Path to the TSV file"),
  make_option(c("--ffn"), type="character", default=NULL, help="Path to the FFN file"),
  make_option(c("--faa"), type="character", default=NULL, help="Path to the FAA file"),
  make_option(c("--fasta"), type="character", default=NULL, help="Path to the FASTA file"),
  make_option(c("--version"), type="character", default=NULL, help="StrainCascade version")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to check if bakta files are present
check_bakta_files <- function(tsv_file, ffn_file, faa_file) {
  if (!file.exists(tsv_file) || !file.exists(ffn_file) || !file.exists(faa_file)) {
    stop("ERROR: One or more Bakta files are missing.")
  }
  return(TRUE)
}

# Function to import fasta file
import_fasta_file <- function(fasta_file) {
  fasta <- readDNAStringSet(fasta_file)
  return(as.data.frame(fasta))
}

# Function to import bakta data
import_bakta_data <- function(tsv_file, ffn_file, faa_file) {
  bakta_data <- read.delim(tsv_file, header = TRUE, sep = "\t", fill = TRUE, skip = 5)
  bakta_nucleotide_seq <- readDNAStringSet(ffn_file)
  bakta_aminoacid_seq <- readAAStringSet(faa_file)
  
  return(list(tsv = bakta_data, ffn = bakta_nucleotide_seq, faa = bakta_aminoacid_seq))
}


# Function to process bakta data
process_bakta_data <- function(bakta_data, fasta_file) {
  # Process .tsv data
  bakta_tsv_processed <- bakta_data$tsv %>%
    dplyr::rename(
      contig_tag_bakta = X.Sequence.Id,
      gene_bakta = Gene,
      gene_product_bakta = Product,
      locus_tag_bakta = Locus.Tag,
      database_references = DbXrefs,
      start_position = Start,
      end_position = Stop,
      strand = Strand,
      type_bakta = Type
    ) %>%
    mutate(
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position),
      length_bp = end_position - start_position + 1,
      source_file_bakta = "annotation_bakta"
    )
  
  # Split the database_references column into individual entries based on commas
  bakta_tsv_databases <- bakta_tsv_processed |>
    select(locus_tag_bakta, database_references) |>
    separate_rows(database_references, sep = ",")
  
  # Trim whitespace
  bakta_tsv_databases$database_references <- trimws(bakta_tsv_databases$database_references)
  
  # Extract everything up to but not including the colon into a new column
  bakta_tsv_databases$column_names <- sub("^(.*?):.*$", "\\1", bakta_tsv_databases$database_references)
  
  # Remove the extracted pattern and colon from the original column
  bakta_tsv_databases$database_references <- sub("^[^:]*?:", "", bakta_tsv_databases$database_references)
  
  # Check if col1 starts with "COG" and conditionally update col2 and col1
  bakta_tsv_databases$column_names <- ifelse(grepl("^COG", bakta_tsv_databases$database_references), paste(bakta_tsv_databases$column_names, "_number", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$database_references <- ifelse(grepl("^COG", bakta_tsv_databases$database_references), sub("^COG", "", bakta_tsv_databases$database_references), bakta_tsv_databases$database_references)
  
  bakta_tsv_databases$column_names <- ifelse(grepl("^UniRef50_", bakta_tsv_databases$database_references), paste(bakta_tsv_databases$column_names, "50", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$database_references <- ifelse(grepl("^UniRef50_", bakta_tsv_databases$database_references), sub("^UniRef50_", "", bakta_tsv_databases$database_references), bakta_tsv_databases$database_references)
  
  bakta_tsv_databases$column_names <- ifelse(grepl("^UniRef90_", bakta_tsv_databases$database_references), paste(bakta_tsv_databases$column_names, "90", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$database_references <- ifelse(grepl("^UniRef90_", bakta_tsv_databases$database_references), sub("^UniRef90_", "", bakta_tsv_databases$database_references), bakta_tsv_databases$database_references)
  
  bakta_tsv_databases$column_names <- ifelse(grepl("^UniRef100_", bakta_tsv_databases$database_references), paste(bakta_tsv_databases$column_names, "100", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$database_references <- ifelse(grepl("^UniRef100_", bakta_tsv_databases$database_references), sub("^UniRef100_", "", bakta_tsv_databases$database_references), bakta_tsv_databases$database_references)
  
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "COG", paste(bakta_tsv_databases$column_names, "_category", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "EC", paste(bakta_tsv_databases$column_names, "_number", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "GO", paste(bakta_tsv_databases$column_names, "_number", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "SO", paste(bakta_tsv_databases$column_names, "_number", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "UniParc", paste(bakta_tsv_databases$column_names, "_ID", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "VFDB", paste(bakta_tsv_databases$column_names, "_ID", sep = ""), bakta_tsv_databases$column_names)
  bakta_tsv_databases$column_names <- ifelse(bakta_tsv_databases$column_names == "KEGG", "K_number", bakta_tsv_databases$column_names)
  
  bakta_tsv_databases[bakta_tsv_databases == ""] <- NA
  bakta_tsv_databases <- bakta_tsv_databases[!is.na(bakta_tsv_databases$column_names),]
  
  bakta_tsv_databases <- bakta_tsv_databases %>%
    pivot_wider(names_from = column_names, 
                values_from = database_references,
                values_fn = list(database_references = ~str_c(.x, collapse = "; ")))
  
  # Add the appendix to all but the first column name
  colnames(bakta_tsv_databases)[-1] <- paste0(colnames(bakta_tsv_databases)[-1], "_bakta")
  
  # Merge the databases back
  bakta_tsv_processed <- bakta_tsv_processed |>
    select(-database_references) |>
    left_join(bakta_tsv_databases, by = "locus_tag_bakta")
  
  
  # Process .ffn data
  bakta_ffn_processed <- as.data.frame(bakta_data$ffn) %>%
    dplyr::rename(nucleotide_code = x) %>%
    mutate(
      locus_tag_bakta = sub(" .*", "", rownames(.)),
      gene_product_bakta = sub("^[^ ]* ", "", rownames(.)),
      rownames = NULL
    )
  
  # Process .ffn data
  bakta_faa_processed <- as.data.frame(bakta_data$faa) %>%
    dplyr::rename(amino_acid_code = x) %>%
    mutate(
      locus_tag_bakta = sub(" .*", "", rownames(.)),
      gene_product_bakta = sub("^[^ ]* ", "", rownames(.)),
      rownames = NULL
    )
  
  # Add contig_tag_fasta as a new column 
  fasta_file$contig_tag_fasta <- rownames(fasta_file)
  
  # Generate contig_tag_bakta and contig_tag_SC
  fasta_file <- fasta_file %>%
    rename(nucleotide_code = x) %>%
    mutate(
      contig_tag_bakta = paste("contig_", row_number(), sep = ""),
      contig_tag_SC = paste("SC_contig_", row_number(), sep = ""),
      rownames = NULL
    )
  
  contig_tag_diccionary <- fasta_file |>
    select(contig_tag_bakta, contig_tag_SC, contig_tag_fasta)
  
  fasta_file <- fasta_file %>%
    select(nucleotide_code, contig_tag_bakta)
  
  ## Index the fasta file
  #
  fasta_file_processed <- fasta_file
  
  # Determine the maximum length of the nucleotide sequences
  max_length <- max(nchar(fasta_file_processed$nucleotide_code))
  
  # Create a new data frame with a "position" column
  position <- 1:max_length
  result <- data.frame(position = position)
  
  # Loop through each contig and split the nucleotide sequence into individual bases
  for (i in 1:nrow(fasta_file_processed)) {
    contig_name <- fasta_file_processed$contig_tag_bakta[i]
    nucleotide_sequence <- strsplit(fasta_file_processed$nucleotide_code[i], split = "")[[1]]
    
    # Pad the sequence with NAs if it's shorter than the max length
    padded_sequence <- c(nucleotide_sequence, rep(NA, max_length - length(nucleotide_sequence)))
    
    # Add this sequence as a new column in the result data frame
    result[[as.character(contig_name)]] <- padded_sequence
  }
  
  squence_data <- merge(bakta_ffn_processed, bakta_faa_processed, by = c("locus_tag_bakta", "gene_product_bakta"), all = TRUE)
  bakta_results_processed <- merge(bakta_tsv_processed, squence_data, by = c("locus_tag_bakta", "gene_product_bakta"), all = TRUE)
  
  ## Retrieve missing nucleotide sequences from the indexed fasta file
  # Function to complement bases
  complement_bases <- function(sequence) {
    # Replace each base with its complement
    complemented_sequence <- chartr("ATCG", "TAGC", sequence)
    return(complemented_sequence)
  }
  
## ALTERNATIVE WAY TO RETRIEVE NUCLEOTIDE SEQUENCE DIRECTLY FROM FASTA FILE (fills nucleotide code where NA)
# Loop through each row in the bakta_results_processed dataframe
for (i in 1:nrow(bakta_results_processed)) {
  # Check if nucleotide_code is NA and start_position and end_position are not NA
  if (is.na(bakta_results_processed$nucleotide_code[i]) && !is.na(bakta_results_processed$start_position[i]) && !is.na(bakta_results_processed$end_position[i])) {

    # Extract start and end positions
    start_pos <- bakta_results_processed$start_position[i]
    end_pos <- bakta_results_processed$end_position[i]

    # Get the contig name (assuming bakta_results_processed has a contig_tag_bakta or similar column matching result)
    contig_name <- bakta_results_processed$contig_tag_bakta[i]

    # Ensure the contig_name exists in the result data frame
    if (contig_name %in% colnames(result)) {
      # Extract the bases corresponding to start_pos and end_pos
      contig_bases <- result[result$position >= start_pos & result$position <= end_pos, contig_name]

      # If contig_bases is a data frame, extract it as a vector
      if (is.data.frame(contig_bases)) {
        contig_bases <- unlist(contig_bases)
      }

      # Ensure valid bases are extracted
      if (length(contig_bases) > 0 && !all(is.na(contig_bases))) {

        # Collapse the valid bases into a single string
        concatenated_string <- paste(contig_bases, collapse = "")

        # Check the strand and modify the sequence if necessary
        if (bakta_results_processed$strand[i] == "-") {
          # Reverse the sequence
          concatenated_string <- rev(strsplit(concatenated_string, NULL)[[1]])
          concatenated_string <- paste(concatenated_string, collapse = "")

          # Complement the sequence
          concatenated_string <- complement_bases(concatenated_string)
        }

        # Assign the final string to bakta_results_processed$nucleotide_code[i]
        bakta_results_processed$nucleotide_code[i] <- concatenated_string
      }
    }
  }
}

  # Add contig_tag_fasta and contig_tag_SC
  bakta_results_processed <- merge(bakta_results_processed, contig_tag_diccionary, by = "contig_tag_bakta", all.x = TRUE)
  
  # Replace empty values with NA
  bakta_results_processed <- bakta_results_processed %>%
    mutate(across(where(is.character), ~ na_if(., "")))
  
  return(list(bakta_results_processed, fasta_file))
}


# Function to replace "NA" strings with actual NA values, excluding date-time columns
replace_na_strings <- function(df) {
  # Identify date-time columns
  date_cols <- sapply(df, inherits, what = c("POSIXt", "Date"))
  
  # Perform replacement on non-date-time columns
  df[!date_cols] <- df[!date_cols] %>%
    mutate(across(everything(), ~ replace(.x, .x == "NA", NA)))
  
  return(df)
}



## Main execution
# Check Bakta files
tryCatch({
  if (!check_bakta_files(opt$tsv, opt$ffn, opt$faa)) {
    stop("Bakta files check failed")
  }
}, error = function(e) {
  cat("ERROR: Failed to check Bakta files.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Import FASTA file
tryCatch({
  fasta_file <- import_fasta_file(opt$fasta)
}, error = function(e) {
  cat("ERROR: Failed to import FASTA file.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Import Bakta data
tryCatch({
  bakta_data <- import_bakta_data(opt$tsv, opt$ffn, opt$faa)
}, error = function(e) {
  cat("ERROR: Failed to import Bakta data.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Process Bakta data
tryCatch({
  processed_data <- process_bakta_data(bakta_data, fasta_file)
  bakta_results <- processed_data[[1]]
}, error = function(e) {
  cat("ERROR: Failed to process Bakta data.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Add metadata
tryCatch({
  bakta_results$software_version <- paste("StrainCascade v", opt$version, sep = "")
  bakta_results$source_file_bakta <- paste(basename(opt$tsv), basename(opt$ffn), basename(opt$faa), sep = "; ")
}, error = function(e) {
  cat("ERROR: Failed to add metadata to results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Select columns
tryCatch({
  columns_to_select <- c(
    "locus_tag_bakta", "contig_tag_bakta", "contig_tag_fasta", "contig_tag_SC", "gene_bakta", "gene_product_bakta",
    "start_position", "end_position", "length_bp", "strand", "type_bakta",
    "nucleotide_code", "amino_acid_code", "EC_number_bakta", "GO_number_bakta",
    "COG_number_bakta", "COG_category_bakta", "K_number_bakta", "SO_number_bakta",
    "UniRef100_bakta", "UniRef90_bakta", "UniRef50_bakta", "UniParc_ID_bakta",
    "RefSeq_bakta", "NCBIProtein_bakta", "VFDB_ID_bakta", "PFAM_bakta",
    "RFAM_bakta", "NCBIFam_bakta", "IS_bakta", "BlastRules_bakta",
    "source_file_bakta", "software_version"
  )
  
  bakta_results <- bakta_results %>%
    select(any_of(columns_to_select))
}, error = function(e) {
  cat("ERROR: Failed to select columns from results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Replace NA strings
tryCatch({
  bakta_results <- replace_na_strings(bakta_results)
}, error = function(e) {
  cat("ERROR: Failed to replace NA strings in results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Remove erroneous (duplicate) rows
bakta_results <- bakta_results[!is.na(bakta_results$contig_tag_bakta) | !is.na(bakta_results$contig_tag_fasta) | !is.na(bakta_results$start_position),]

# Save results
tryCatch({
  output_file <- file.path(opt$output_dir, "bakta_results.qs")
  qsave(bakta_results, output_file, preset = "archive")
  cat("SUCCESS: Results saved to", output_file, "\n")
}, error = function(e) {
  cat("ERROR: Failed to save results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})