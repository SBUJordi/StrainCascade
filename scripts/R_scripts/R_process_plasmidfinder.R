#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_plasmidfinder.R

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(qs)
  library(optparse)
  library(purrr)
  library(tidyr)
  library(Biostrings)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="character"),
  make_option(c("--json"), type="character", default=NULL, 
              help="Path to the PlasmidFinder JSON file", metavar="character"),
  make_option(c("--fasta"), type="character", default=NULL,
              help="Path to FASTA file", metavar="character"),
  make_option(c("--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$json) || is.null(opt$fasta) || is.null(opt$version)) {
  print_help(opt_parser)
  stop("All arguments must be supplied.", call.=FALSE)
}

# Function to check if file exists
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
}

# Import and process plasmid data from JSON file
import_and_process_plasmid_data <- function(json_file_path) {
  check_file_exists(json_file_path)
  data <- fromJSON(json_file_path)
  
  # Process user input
  user_input <- data$user_input %||% data$plasmidfinder$user_input
  
  # Process run info
  run_info <- data$run_info %||% data$plasmidfinder$run_info
  
  # Process results
  results <- data$plasmid_data %||% data$plasmidfinder$results
  
  # Extract plasmid data
  plasmid_data <- if (is.null(data$plasmid_data)) {
    lapply(names(results), function(category) {
      lapply(names(results[[category]]), function(subcategory) {
        list(
          category = category,
          subcategory = subcategory,
          result = results[[category]][[subcategory]]
        )
      })
    }) %>% unlist(recursive = FALSE)
  } else {
    results
  }
  
  return(list(
    user_input = user_input,
    run_info = run_info,
    plasmid_data = plasmid_data
  ))
}

# Function to process plasmid result
process_plasmid_result <- function(result) {
  if (is.character(result) && result == "No hit found") {
    return(data.frame(hit_found = FALSE, stringsAsFactors = FALSE))
  } else if (is.list(result)) {
    # Extract the first (and only) element of the result list if it exists
    hit <- if (length(result) > 0) result[[1]] else result
    if (is.list(hit)) {
      # Convert the hit to a data frame
      df <- as.data.frame(hit, stringsAsFactors = FALSE)
      df$hit_found <- TRUE
      return(df)
    }
  }
  # Return NA for any other case
  return(data.frame(hit_found = NA, stringsAsFactors = FALSE))
}

# Function to process a single plasmid entry
process_plasmid_entry <- function(entry) {
  result_df <- process_plasmid_result(entry$result)
  cbind(
    data.frame(
      category = entry$category,
      subcategory = entry$subcategory,
      stringsAsFactors = FALSE
    ),
    result_df
  )
}

# Function to import FASTA file
import_fasta_file <- function(fasta_file) {
  fasta_data <- readDNAStringSet(fasta_file)
  return(fasta_data)
}

# Function to extract nucleotide sequence
extract_nucleotide_sequence <- function(fasta_data, contig, start, end, strand) {
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


# Create a dataframe from processed plasmid data
create_plasmid_dataframe <- function(processed_data, json_file_path, fasta_data) {
  plasmid_df <- map_df(processed_data$plasmid_data, process_plasmid_entry)
  
  plasmid_df$date <- processed_data$run_info$date
  plasmid_df$method <- processed_data$user_input$method
  plasmid_df$file_format <- processed_data$user_input$`file_format` %||% processed_data$user_input$file_format
  
  # Rename columns flexibly
  new_names <- sapply(names(plasmid_df), function(col) {
    paste0(col, "_plasmidfinder")
  })
  
  required_columns <- c("contig_name", "positions_in_contig", "position_in_ref")
  
  if (all(required_columns %in% names(plasmid_df)) && !all(plasmid_df == FALSE)) {
    plasmid_df <- plasmid_df %>%
      rename_with(~ new_names[.x], .cols = everything()) %>%
      dplyr::rename(contig_tag_fasta = contig_name_plasmidfinder) %>%
      mutate(
        start_position = as.integer(sub("(\\d+)\\.\\.(\\d+)", "\\1", positions_in_contig_plasmidfinder)),
        end_position = as.integer(sub("(\\d+)\\.\\.(\\d+)", "\\2", positions_in_contig_plasmidfinder)),
        start_position_ref = as.integer(sub("(\\d+)\\.\\.(\\d+)", "\\1", position_in_ref_plasmidfinder)),
        end_position_ref = as.integer(sub("(\\d+)\\.\\.(\\d+)", "\\2", position_in_ref_plasmidfinder))
      ) %>%
      rowwise() %>%
      mutate(
        nucleotide_code = extract_nucleotide_sequence(
          fasta_data, 
          contig_tag_fasta, 
          start_position, 
          end_position, 
          "+"  # Assuming forward strand, adjust if needed
        ),
        start_position = as.numeric(start_position),
        end_position = as.numeric(end_position),
        start_position_ref = as.numeric(start_position_ref),
        end_position_ref = as.numeric(end_position_ref)
      ) %>%
      ungroup() %>%
      select(-positions_in_contig_plasmidfinder, -position_in_ref_plasmidfinder, -date_plasmidfinder)
  }
  
  plasmid_df$source_file_plasmidfinder = basename(json_file_path)
  
  return(plasmid_df)
}

# Main execution
main <- function() {
  # Import and process plasmid data
  processed_plasmid_data <- import_and_process_plasmid_data(opt$json)
  
  # Import FASTA data
  fasta_data <- import_fasta_file(opt$fasta)
  
  # Create plasmid dataframe
  plasmid_profile <- create_plasmid_dataframe(processed_plasmid_data, opt$json, fasta_data)
  
  # Add comments to each column
# comments_list <- list(
#   category_plasmidfinder = "Category of the plasmid",
#   subcategory_plasmidfinder = "Subcategory of the plasmid",
#   hit_found_plasmidfinder = "Whether a hit was found",
#   method_plasmidfinder = "Method used for analysis",
#   file_format_plasmidfinder = "Format of the input file",
#   source_file_plasmidfinder = "Source JSON file for PlasmidFinder data",
#   contig_tag_fasta = "Contig tag from input assembly (.fasta)",
#   start_position = "Start position of the plasmid on the contig",
#   end_position = "End position of the plasmid on the contig",
#   nucleotide_code = "Nucleotide sequence of the plasmid"
# )
  
  # Assign comments to plasmid_profile columns
# for (col in names(comments_list)) {
#   if (col %in% names(plasmid_profile)) {
#     comment(plasmid_profile[[col]]) <- comments_list[[col]]
#   } else {
#     warning(paste("Column '", col, "' not found in plasmid_profile. Skipping comment assignment.", sep = ""))
#   }
# }
  
  # Save the results
  output_file <- file.path(opt$output_dir, "plasmidfinder_results.qs")
  qsave(plasmid_profile, output_file, preset = "archive")
  cat(paste("SUCCESS: Results saved to", output_file, "\n"))
}

# Run the main function
tryCatch({
  main()
}, error = function(e) {
  cat(paste("ERROR:", e$message, "\n"))
  quit(status = 1)
})