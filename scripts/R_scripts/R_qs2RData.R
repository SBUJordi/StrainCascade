#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

set.seed(42)

# Load required libraries
suppressPackageStartupMessages({
  library(qs)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = ".", 
              help = "Directory containing .qs files to import"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".", 
              help = "Directory to save RData object and RProject"),
  make_option(c("-s", "--sample_name"), type = "character", default = "default",
              help = "Sample name to be included in output file names")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list = option_list))

# Function to import .qs file if it exists and process dataframe if applicable
import_and_process_qs <- function(file_path) {
  if (file.exists(file_path)) {
    cat(paste("Importing:", file_path, "\n"))
    data <- qs::qread(file_path)
    
    # Check if data is a dataframe
    if (is.data.frame(data)) {
      columns_to_check <- c("contig_tag_fasta", "start_position", "end_position", "length_bp", "strand", "nucleotide_code", "amino_acid_code")
      char_columns <- c("contig_tag_fasta", "strand", "nucleotide_code", "amino_acid_code")
      numeric_columns <- c("start_position", "end_position", "length_bp")
      
      for (col in columns_to_check) {
        if (col %in% names(data)) {
          if (col %in% char_columns) {
            data[[col]] <- as.character(data[[col]])
          } else if (col %in% numeric_columns) {
            data[[col]] <- as.numeric(data[[col]])
          }
        }
      }
      cat(paste("Processed dataframe:", file_path, "\n"))
    }
    
    return(data)
  } else {
    cat(paste("File not found:", file_path, "\n"))
    return(NULL)
  }
}

# List of possible .qs files
qs_files <- c(
  "selected_assembly_results.qs",
  "circlator_results.qs",
  "annotation_results_integrated.qs",
  "annotation_results_aggregated.qs",
  "microbeannotator_KEGG_modules_results.qs",
  "checkm2_results.qs",
  "quast_results.qs",
  "gtdbtk_classification_results.qs",
  "gtdbtk_classification_newick_tree.qs",
  "gtdbtk_25_tree_plot.qs",
  "dbcan3_results.qs",
  "plasmidfinder_results.qs",
  "amrfinderplus_results.qs",
  "resfinder_results.qs",
  "islandpath_results.qs",
  "prokka_results.qs",
  "bakta_results.qs",
  "microbeannotator_results.qs",
  "resfinder_AMR_profile.qs",
  "virsorter2_results.qs",
  "isescan_results.qs"
)

# Import available .qs files as individual objects
for (file in qs_files) {
  file_path <- file.path(opt$input_dir, file)
  
  # Skip aggregated if integrated is available
  if (file == "annotation_results_aggregated.qs" && exists("annotation_results_integrated")) {
    next
  }

  data <- import_and_process_qs(file_path)
  if (!is.null(data)) {
    # Remove the file extension and create a variable name dynamically
    var_name <- gsub(".qs", "", file)
    
    # Assign the data to an object in the environment with the constructed name
    assign(var_name, data, envir = .GlobalEnv)
  }
}

# Check if any files were imported
if (!any(sapply(qs_files, function(f) exists(gsub(".qs", "", f), envir = .GlobalEnv)))) {
  cat("No .qs files were found or imported.\n")
  quit(save = "no", status = 0)
}

# Save imported data objects as RData object
save_name <- paste0("StrainCascade_Results_", opt$sample_name, "_", format(Sys.time(), "%y%m%d"), ".RData")
save(list = ls(pattern = ".*_results|resfinder_AMR_profile|gtdbtk_25_tree_plot"), file = file.path(opt$output_dir, save_name))
cat(paste("Saved imported data as:", save_name, "\n"))

# Create RProject
project_name <- paste0("StrainCascade_Analysis_", opt$sample_name)
project_dir <- file.path(opt$output_dir)
setwd(project_dir)

# Create .Rproj file
rproj_content <- c(
  "Version: 1.0",
  "",
  "RestoreWorkspace: Default",
  "SaveWorkspace: Default",
  "AlwaysSaveHistory: Default",
  "",
  "EnableCodeIndexing: Yes",
  "UseSpacesForTab: Yes",
  "NumSpacesForTab: 2",
  "Encoding: UTF-8",
  "",
  "RnwWeave: Sweave",
  "LaTeX: pdfLaTeX"
)
writeLines(rproj_content, paste0(project_name, ".Rproj"))

cat(paste("Created RProject:", project_dir, "\n"))
cat("RData file saved in the project directory.\n")