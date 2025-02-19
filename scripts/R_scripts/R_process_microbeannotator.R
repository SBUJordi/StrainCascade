#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_microbeannotator.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(qs)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character", default=NULL, help="Output directory for results"),
  make_option(c("--annot"), type="character", default=NULL, help="Path to the MicrobeAnnotator annotation file"),
  make_option(c("--version"), type="character", default=NULL, help="StrainCascade version")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to check if MicrobeAnnotator file is present
check_microbeannotator_file <- function(annot_file) {
  if (!file.exists(annot_file)) {
    stop("ERROR: MicrobeAnnotator annotation file is missing.")
  }
  return(TRUE)
}

# Function to determine the annotation tool used
determine_annotation_tool <- function(file_path) {
  if (grepl("bakta", file_path)) {
    return("bakta")
  } else if (grepl("prokka", file_path)) {
    return("prokka")
  } else {
    return("non_StainCascade_tool")
  }
}

# Function to import MicrobeAnnotator data
import_microbeannotator_data <- function(file_path) {
  microbeannotator_data <- read.delim(file_path, header = TRUE, sep = "\t", fill = TRUE)
  return(microbeannotator_data)
}

# Function to process MicrobeAnnotator data
process_microbeannotator_data <- function(microbeannotator_data, annotation_tool, file_path) {
  locus_tag_col <- paste0("locus_tag_", annotation_tool)
  
  processed_data <- microbeannotator_data %>%
    select(query_id, ko_number, ko_product, database) %>%
    rename(
      !!locus_tag_col := query_id,
      gene_product_microbeannotator = ko_product,
      K_number_microbeannotator = ko_number,
      database_microbeannotator = database
    ) %>%
    mutate(
      source_file_microbeannotator = basename(file_path),
      software_version = paste("StrainCascade v", opt$version, sep = "")
    )
  
  return(processed_data)
}

# Function to replace "NA" strings with actual NA values, excluding date-time columns
replace_na_strings <- function(df) {
  date_cols <- sapply(df, inherits, what = c("POSIXt", "Date"))
  df[!date_cols] <- df[!date_cols] %>%
    mutate(across(everything(), ~ replace(.x, .x == "NA", NA)))
  return(df)
}

# Main execution
tryCatch({
  if (check_microbeannotator_file(opt$annot)) {
    annotation_tool <- determine_annotation_tool(opt$annot)
    microbeannotator_data <- import_microbeannotator_data(opt$annot)
    processed_data <- process_microbeannotator_data(microbeannotator_data, annotation_tool, opt$annot)
    
    microbeannotator_results <- processed_data
    
    # Extract the EC numbers from the gene_product_microbeannotator column
    microbeannotator_results$EC_number_microbeannotator <- ifelse(
      grepl("\\[EC:[0-9.-]+( [0-9.-]+)*\\]", microbeannotator_results$gene_product_microbeannotator),
      gsub(".*\\[EC:([0-9.-]+( [0-9.-]+)*)\\].*", "\\1", microbeannotator_results$gene_product_microbeannotator),
      NA
    )
    
    # Replace spaces with semicolons to match the desired format
    microbeannotator_results$EC_number_microbeannotator <- gsub(
      pattern = " ",
      replacement = "; ",
      x = microbeannotator_results$EC_number_microbeannotator
    )
    
    # Remove the EC numbers and brackets from the gene_product_microbeannotator column
    microbeannotator_results$gene_product_microbeannotator <- gsub(
      pattern = "\\s*\\[EC:[0-9.-]+( [0-9.-]+)*\\]",
      replacement = "",
      x = microbeannotator_results$gene_product_microbeannotator
    )
    
    microbeannotator_results <- microbeannotator_results %>%
      select(
        starts_with("locus_tag"),
        gene_product_microbeannotator,
        EC_number_microbeannotator,
        K_number_microbeannotator,
        database_microbeannotator,
        source_file_microbeannotator,
        software_version
      )
    
    microbeannotator_results <- replace_na_strings(microbeannotator_results)
    
    # Save the results
    output_file <- file.path(opt$output_dir, "microbeannotator_results.qs")
    qsave(microbeannotator_results, output_file, preset = "archive")
    cat("SUCCESS: Results saved to", output_file, "\n")
  }
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})