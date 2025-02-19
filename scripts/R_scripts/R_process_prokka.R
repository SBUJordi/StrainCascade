#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_prokka.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
 library(Biostrings)
  library(stringr)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(qs)
  library(optparse)
})


# Parse command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character", default=NULL, help="Output directory for results"),
  make_option(c("--gff"), type="character", default=NULL, help="Path to the GFF file"),
  make_option(c("--ffn"), type="character", default=NULL, help="Path to the FFN file"),
  make_option(c("--faa"), type="character", default=NULL, help="Path to the FAA file"),
  make_option(c("--fasta"), type="character", default=NULL, help="Path to the FASTA file"),
  make_option(c("--version"), type="character", default=NULL, help="StrainCascade version")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to check if prokka files are present
check_prokka_files <- function(gff_file, ffn_file, faa_file) {
  if (!file.exists(gff_file) || !file.exists(ffn_file) || !file.exists(faa_file)) {
    stop("ERROR: One or more Prokka files are missing.")
  }
  return(TRUE)
}

# Function to import fasta file
import_fasta_file <- function(fasta_file) {
  fasta <- readDNAStringSet(fasta_file)
  return(as.data.frame(fasta))
}

# Function to import prokka data
import_prokka_data <- function(gff_file, ffn_file, faa_file) {
  prokka_gff <- readGFF(gff_file, version=3, columns=NULL, tags=NULL, filter=NULL, nrows=-1, raw_data=FALSE)
  prokka_nucleotide_seq <- readDNAStringSet(ffn_file)
  prokka_aminoacid_seq <- readAAStringSet(faa_file)
  
  return(list(gff = prokka_gff, ffn = prokka_nucleotide_seq, faa = prokka_aminoacid_seq))
}

# Function to process prokka data
process_prokka_data <- function(prokka_data, fasta_file) {
  
  # Process .gff data
  prokka_gff <- prokka_data$gff
  prokka_gff_processed <- data.frame(
    locus_tag_prokka = prokka_gff@listData$locus_tag,
    contig_tag_prokka = prokka_gff@listData$seqid,  
    gene_prokka = prokka_gff@listData$gene,
    gene_product_prokka = prokka_gff@listData$product,
    start_position = prokka_gff@listData$start,
    end_position = prokka_gff@listData$end,
    length_bp = prokka_gff@listData$end - prokka_gff@listData$start + 1,  # Calculate length in bp
    strand = prokka_gff@listData$strand,
    type_prokka = prokka_gff@listData$type,
    COG_number_prokka = prokka_gff@listData$db_xref,
    EC_number_prokka = prokka_gff@listData$eC_number,
    temp = prokka_gff@listData$ID  # Assuming 'ID' can represent the contig tag
  ) |>
  mutate(
    start_position = as.numeric(start_position),
    end_position = as.numeric(end_position)
  )
  
  # Remove all gene rows (keep CDS etc.) to avoid duplicates
  prokka_gff_processed <- prokka_gff_processed %>%
    filter(!grepl("_gene", temp)) |>
    select(-temp)
  
  # Remove everything before and including the last underscore
  prokka_gff_processed$contig_tag_prokka <- sub(".*_", "", prokka_gff_processed$contig_tag_prokka)
  
  # Process .ffn data
  prokka_ffn_processed <- as.data.frame(prokka_data$ffn) %>%
    rename(nucleotide_code = x) %>%
    mutate(
      locus_tag_prokka = sub(" .*", "", rownames(.)),
      gene_product_prokka = sub("^[^ ]* ", "", rownames(.)),
      rownames = NULL
    )
  
  
  # Process .faa data
  prokka_faa_processed <- as.data.frame(prokka_data$faa) %>%
    rename(amino_acid_code = x) %>%
    mutate(
      locus_tag_prokka = sub(" .*", "", rownames(.)),
      gene_product_prokka = sub("^[^ ]* ", "", rownames(.)),
      rownames = NULL
    )
  
  # Add contig_tag_fasta as a new column 
  fasta_file$contig_tag_fasta <- rownames(fasta_file)
  
  # Generate contig_tag_prokka and contig_tag_SC
  fasta_file <- fasta_file %>%
    rename(nucleotide_code = x) %>%
    mutate(
      contig_tag_prokka = row_number(),
      contig_tag_SC = paste("SC_contig_", row_number(), sep = ""),
      rownames = NULL
    )
  
  contig_tag_diccionary <- fasta_file |>
    select(contig_tag_prokka, contig_tag_SC, contig_tag_fasta)
  
  fasta_file <- fasta_file %>%
    select(nucleotide_code, contig_tag_prokka)
  
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
    contig_name <- fasta_file_processed$contig_tag_prokka[i]
    nucleotide_sequence <- strsplit(fasta_file_processed$nucleotide_code[i], split = "")[[1]]
    
    # Pad the sequence with NAs if it's shorter than the max length
    padded_sequence <- c(nucleotide_sequence, rep(NA, max_length - length(nucleotide_sequence)))
    
    # Add this sequence as a new column in the result data frame
    result[[as.character(contig_name)]] <- padded_sequence
  }
  
  squence_data <- merge(prokka_ffn_processed, prokka_faa_processed, by = c("locus_tag_prokka", "gene_product_prokka"), all = TRUE)
  prokka_results_processed <- merge(prokka_gff_processed, squence_data, by = c("locus_tag_prokka", "gene_product_prokka"), all = TRUE)
  
  ## Retrieve missing nucleotide sequences from the indexed fasta file
  # Function to complement bases
  complement_bases <- function(sequence) {
    # Replace each base with its complement
    complemented_sequence <- chartr("ATCG", "TAGC", sequence)
    return(complemented_sequence)
  }
  
 
## ALTERNATIVE WAY TO RETRIEVE NUCLEOTIDE SEQUENCE DIRECTLY FROM FASTA FILE (fills nucleotide code where NA)
# Loop through each row in the prokka_results_processed dataframe
for (i in 1:nrow(prokka_results_processed)) {
  # Check if nucleotide_code is NA and start_position and end_position are not NA
  if (is.na(prokka_results_processed$nucleotide_code[i]) && !is.na(prokka_results_processed$start_position[i]) && !is.na(prokka_results_processed$end_position[i])) {

    # Extract start and end positions
    start_pos <- prokka_results_processed$start_position[i]
    end_pos <- prokka_results_processed$end_position[i]

    # Get the contig name (assuming prokka_results_processed has a contig_tag_prokka or similar column matching result)
    contig_name <- prokka_results_processed$contig_tag_prokka[i]

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
        if (prokka_results_processed$strand[i] == "-") {
          # Reverse the sequence
          concatenated_string <- rev(strsplit(concatenated_string, NULL)[[1]])
          concatenated_string <- paste(concatenated_string, collapse = "")

          # Complement the sequence
          concatenated_string <- complement_bases(concatenated_string)
        }

        # Assign the final string to prokka_results_processed$nucleotide_code[i]
        prokka_results_processed$nucleotide_code[i] <- concatenated_string
      }
    }
  }
}
 
  # Add contig_tag_fasta and contig_tag_SC
  prokka_results_processed <- merge(prokka_results_processed, contig_tag_diccionary, by = "contig_tag_prokka", all.x = TRUE)
  
  # Replace empty values with NA
  prokka_results_processed <- prokka_results_processed %>%
    mutate(across(where(is.character), ~ na_if(., "")))
  
  return(list(prokka_results_processed, fasta_file))
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
# Check Prokka files
tryCatch({
  if (!check_prokka_files(opt$gff, opt$ffn, opt$faa)) {
    stop("Prokka files check failed")
  }
}, error = function(e) {
  cat("ERROR: Failed to check Prokka files.\n")
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

# Import Prokka data
tryCatch({
  prokka_data <- import_prokka_data(opt$gff, opt$ffn, opt$faa)
}, error = function(e) {
  cat("ERROR: Failed to import Prokka data.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Process Prokka data
tryCatch({
  processed_data <- process_prokka_data(prokka_data, fasta_file)
  prokka_results <- processed_data[[1]]
}, error = function(e) {
  cat("ERROR: Failed to process Prokka data.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Add metadata
tryCatch({
  prokka_results$software_version <- paste("StrainCascade v", opt$version, sep = "")
  prokka_results$source_file_prokka <- paste(basename(opt$gff), sep = "; ")
}, error = function(e) {
  cat("ERROR: Failed to add metadata to results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Select columns
tryCatch({
  columns_to_select <- c(
    "locus_tag_prokka", "contig_tag_prokka", "contig_tag_fasta", "contig_tag_SC", "gene_prokka", "gene_product_prokka",
    "start_position", "end_position", "length_bp", "strand", "type_prokka",
    "nucleotide_code", "amino_acid_code", "EC_number_prokka", "COG_number_prokka",
    "source_file_prokka", "software_version"
  )
  
  prokka_results <- prokka_results %>%
    select(any_of(columns_to_select))
}, error = function(e) {
  cat("ERROR: Failed to select columns from results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Replace NA strings
tryCatch({
  prokka_results <- replace_na_strings(prokka_results)
}, error = function(e) {
  cat("ERROR: Failed to replace NA strings in results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})

# Remove erroneous (duplicate) rows
prokka_results <- prokka_results[!is.na(prokka_results$contig_tag_prokka) | !is.na(prokka_results$contig_tag_fasta) | !is.na(prokka_results$start_position),]

# Save results
tryCatch({
  output_file <- file.path(opt$output_dir, "prokka_results.qs")
  qsave(prokka_results, output_file, preset = "archive")
  cat("SUCCESS: Results saved to", output_file, "\n")
}, error = function(e) {
  cat("ERROR: Failed to save results.\n")
  cat(sprintf("Error message: %s\n", e$message))
  quit(save = "no", status = 1)
})