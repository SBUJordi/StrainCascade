#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_deepfri.R

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
  make_option(c("--output_dir"), type="character", default=NULL, help="Output directory for results"),
  make_option(c("--input_dir"), type="character", default=NULL, help="Input directory with DeepFRI CSV files"),
  make_option(c("--sample_name"), type="character", default=NULL, help="Sample name"),
  make_option(c("--faa"), type="character", default=NULL, help="Path to the FAA file"),
  make_option(c("--version"), type="character", default=NULL, help="StrainCascade version")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to determine annotation tool from filename
determine_annotation_tool <- function(file_path) {
  if (grepl("bakta", file_path)) {
    return("bakta")
  } else if (grepl("prokka", file_path)) {
    return("prokka")
  } else {
    return("unknown")
  }
}

# Function to import protein sequences and extract locus tags
import_protein_sequences <- function(faa_file, annotation_tool) {
  aa_sequences <- readAAStringSet(faa_file)
  
  faa_data <- data.frame(
    protein_id = names(aa_sequences),
    amino_acid_code = as.character(aa_sequences),
    stringsAsFactors = FALSE
  )
  
  # Extract locus tag from protein ID
  # Bakta format: "locus_tag info..."
  # Prokka format: "locus_tag info..."
  # DeepFRI uses only the locus tag as protein ID, so update protein_id to match
  faa_data$protein_id <- sub("\\s.*", "", faa_data$protein_id)
  
  locus_tag_col <- paste0("locus_tag_", annotation_tool)
  faa_data[[locus_tag_col]] <- faa_data$protein_id
  
  return(faa_data)
}

# Function to import DeepFRI predictions
import_deepfri_predictions <- function(csv_file, ontology) {
  if (!file.exists(csv_file)) {
    return(NULL)
  }
  
  tryCatch({
    # Skip comment lines starting with # (DeepFRI adds header comment)
    predictions <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
    
    if (nrow(predictions) == 0) {
      return(NULL)
    }
    
    # Standardize column names
    colnames(predictions) <- c("protein_id", "term_id", "score", "term_name")
    
    # Add ontology information
    predictions$ontology <- ontology
    
    return(predictions)
  }, error = function(e) {
    cat(sprintf("Warning: Failed to read %s: %s\n", csv_file, e$message))
    return(NULL)
  })
}

# Function to process EC numbers
process_ec_predictions <- function(ec_data) {
  if (is.null(ec_data) || nrow(ec_data) == 0) {
    return(data.frame())
  }
  
  # Get best EC prediction (highest score) for each protein
  # DeepFRI uses structure-based prediction (orthogonal to sequence homology tools)
  # so we use only the most confident prediction rather than matching vote count
  ec_best <- ec_data %>%
    dplyr::group_by(.data$protein_id) %>%
    dplyr::slice_max(order_by = .data$score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select("protein_id", "term_id", "score") %>%
    dplyr::rename(
      EC_number_deepfri = "term_id",
      EC_number_deepfri_score = "score"
    )
  
  # Also keep all EC predictions for reference (with scores)
  ec_all <- ec_data %>%
    group_by(.data$protein_id) %>%
    summarise(
      EC_number_deepfri_all = paste(.data$term_id, collapse = "; "),
      EC_number_deepfri_all_scores = paste(sprintf("%.5f", .data$score), collapse = "; "),
      .groups = "drop"
    )
  
  # Join best and all predictions
  ec_processed <- ec_best %>%
    left_join(ec_all, by = "protein_id")
  
  return(ec_processed)
}

# Function to process GO terms
process_go_predictions <- function(go_data, ontology_name) {
  if (is.null(go_data) || nrow(go_data) == 0) {
    return(data.frame())
  }
  
  # Use full ontology names instead of abbreviations
  full_name <- switch(ontology_name,
    "MF" = "molecular_function",
    "BP" = "biological_process",
    "CC" = "cellular_component",
    ontology_name
  )
  col_prefix <- paste0("GO_", full_name)
  
  # Group by protein and aggregate GO terms
  go_processed <- go_data %>%
    group_by(.data$protein_id) %>%
    summarise(
      !!paste0(col_prefix, "_number_deepfri") := paste(.data$term_id, collapse = "; "),
      !!paste0(col_prefix, "_score_deepfri") := paste(sprintf("%.5f", .data$score), collapse = "; "),
      !!paste0(col_prefix, "_names_deepfri") := paste(.data$term_name, collapse = "; "),
      .groups = "drop"
    )
  
  return(go_processed)
}

# Function to extract gene products from predictions
extract_gene_products <- function(all_predictions) {
  if (is.null(all_predictions) || nrow(all_predictions) == 0) {
    return(data.frame())
  }
  
  # Get highest scoring prediction per protein
  gene_products <- all_predictions %>%
    dplyr::group_by(.data$protein_id) %>%
    dplyr::slice_max(order_by = .data$score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select("protein_id", "term_name", "score") %>%
    dplyr::rename(
      gene_product_deepfri = "term_name",
      gene_product_deepfri_score = "score"
    )
  
  return(gene_products)
}

# Function to replace "NA" strings with actual NA values
replace_na_strings <- function(df) {
  date_cols <- sapply(df, inherits, what = c("POSIXt", "Date"))
  na_pattern <- "^(NA|N/A|N/a|n/A|n/a|)$"
  
  df[!date_cols] <- df[!date_cols] %>%
    mutate(across(everything(), ~ replace(.x, grepl(na_pattern, .x), NA)))
  
  return(df)
}

# Main execution
tryCatch({
  # Check if required arguments are provided
  if (is.null(opt$output_dir) || is.null(opt$input_dir) || 
      is.null(opt$sample_name) || is.null(opt$faa) || is.null(opt$version)) {
    stop("All arguments must be supplied.")
  }
  
  # Determine annotation tool
  annotation_tool <- determine_annotation_tool(opt$faa)
  if (annotation_tool == "unknown") {
    stop("Could not determine annotation tool from FAA filename.")
  }
  
  cat(sprintf("Using annotation tool: %s\n", annotation_tool))
  
  # Import protein sequences
  protein_data <- import_protein_sequences(opt$faa, annotation_tool)
  cat(sprintf("Imported %d protein sequences\n", nrow(protein_data)))
  
  # Import DeepFRI predictions for all ontologies
  ontologies <- list(
    mf = "MF",
    bp = "BP",
    cc = "CC",
    ec = "EC"
  )
  
  all_predictions <- list()
  
  for (ont in names(ontologies)) {
    ont_upper <- ontologies[[ont]]
    csv_file <- file.path(opt$input_dir, 
                          sprintf("%s_deepfri_%s_%s_predictions.csv", 
                                  opt$sample_name, ont, ont_upper))
    
    if (file.exists(csv_file)) {
      predictions <- import_deepfri_predictions(csv_file, ont)
      if (!is.null(predictions)) {
        all_predictions[[ont]] <- predictions
        cat(sprintf("Imported %d predictions for ontology %s\n", 
                    nrow(predictions), ont))
      }
    }
  }
  
  if (length(all_predictions) == 0) {
    stop("No DeepFRI predictions found.")
  }
  
  # Combine all predictions for gene product extraction
  combined_predictions <- bind_rows(all_predictions)
  
  # Extract gene products (highest scoring prediction)
  gene_products <- extract_gene_products(combined_predictions)
  
  # Process EC numbers
  ec_results <- if ("ec" %in% names(all_predictions)) {
    process_ec_predictions(all_predictions$ec)
  } else {
    data.frame()
  }
  
  # Process GO terms
  go_mf_results <- if ("mf" %in% names(all_predictions)) {
    process_go_predictions(all_predictions$mf, "MF")
  } else {
    data.frame()
  }
  
  go_bp_results <- if ("bp" %in% names(all_predictions)) {
    process_go_predictions(all_predictions$bp, "BP")
  } else {
    data.frame()
  }
  
  go_cc_results <- if ("cc" %in% names(all_predictions)) {
    process_go_predictions(all_predictions$cc, "CC")
  } else {
    data.frame()
  }
  
  # Merge all results
  deepfri_results <- protein_data
  
  if (nrow(gene_products) > 0) {
    deepfri_results <- deepfri_results %>%
      left_join(gene_products, by = "protein_id")
  }
  
  if (nrow(ec_results) > 0) {
    deepfri_results <- deepfri_results %>%
      left_join(ec_results, by = "protein_id")
  }
  
  if (nrow(go_mf_results) > 0) {
    deepfri_results <- deepfri_results %>%
      left_join(go_mf_results, by = "protein_id")
  }
  
  if (nrow(go_bp_results) > 0) {
    deepfri_results <- deepfri_results %>%
      left_join(go_bp_results, by = "protein_id")
  }
  
  if (nrow(go_cc_results) > 0) {
    deepfri_results <- deepfri_results %>%
      left_join(go_cc_results, by = "protein_id")
  }
  
  # Add metadata
  deepfri_results$source_file_deepfri <- sprintf("%s_deepfri_predictions", opt$sample_name)
  deepfri_results$software_version <- paste0("StrainCascade v", opt$version)
  
  # Remove protein_id column (keep locus_tag instead)
  deepfri_results <- deepfri_results %>%
    select(-protein_id)
  
  # Replace NA strings
  deepfri_results <- replace_na_strings(deepfri_results)
  
  # Save results
  output_file <- file.path(opt$output_dir, "deepfri_results.qs")
  qsave(deepfri_results, output_file, preset = "archive")
  
  cat(sprintf("SUCCESS: Results saved to %s\n", output_file))
  cat(sprintf("Processed %d proteins with DeepFRI annotations\n", nrow(deepfri_results)))
  
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 1)
})
