# BIOMES_result_processing_for_R.r - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# List of required packages
packages <- c("Biostrings", "genbankr", "stringr", "dplyr", "tidyr", "stringdist")

# Function to check and install packages if necessary in the library associated with the current R session (in this case in the miniconda env.)
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    # Check if the package is a Bioconductor package
    if (pkg %in% c("Biostrings", "genbankr")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, lib = .libPaths()[1])
    } else {
      install.packages(pkg, lib = .libPaths()[1], dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Check and install required packages
sapply(packages, check_and_install)

# Load the required packages
library(Biostrings)
library(genbankr)
library(stringr)
library(dplyr)
library(tidyr)
library(stringdist)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define the base directory
main_results_dir_abs <- args[1]  # assuming main_results_dir_abs is the 1. argument
R_analysis_dir <- args[2]  # assuming R_analysis_dir is the 2. argument
sample_name <- args[3]

# Save time of execution
date_of_execution <- Sys.Date()

# Find and extract the chosen assembly name
chosen_assembly_paths <- list.files(path = main_results_dir_abs, pattern = "_chosen_assembly.fasta", recursive = TRUE, full.names = FALSE)
chosen_assembly <- unique(sub("_chosen_assembly.fasta.*", "", chosen_assembly_paths))
chosen_assembly <- sub(".*?\\Qassembly_\\E", "", chosen_assembly, perl = TRUE)
chosen_assembly <- paste(chosen_assembly, "_assembly", sep = "")

# Define the file patterns for files that should be imported
file_patterns <- c("annotation_prokka.tsv", "annotation_bakta.tsv", "bakta.hypotheticals.tsv", 
                   "annotation_prokka.ffn", "annotation_bakta.ffn", 
                   "annotation_prokka.faa", "annotation_bakta.faa",
                   "annotation_prokka.gbk", "annotation_bakta.gbff",
                   "annotation_bakta_microbeannotator.annot",
                   "module_completeness.tab", "assembly_evaluation.tsv", 
                   "gtdbtk_organism_name.txt", "taxonomy_gtdbtk_summary.tsv", "assembly_evaluation.tsv")

# Initialize an empty vector to store the file paths
file_paths <- c()

# -----
## Import annotation data
# Data import and basic processing
# Loop over the file patterns
for (pattern in file_patterns) {
  # Get a list of file paths that match the pattern
  paths <- list.files(path = main_results_dir_abs, pattern = pattern, recursive = TRUE, full.names = TRUE)
  
  # Add the file paths to the vector
  file_paths <- c(file_paths, paths)
}

# Exclude files that start with "~$" (-> temporary files)
file_paths <- file_paths[grep("~\\$", file_paths, invert = TRUE)]

# Loop over the file paths
for (file_path in file_paths) {
  # Determine the number of lines to skip based on the file name
  if (grepl("annotation_bakta.tsv", file_path)) {
    skip <- 5
  } else if (grepl("bakta.hypotheticals.tsv", file_path)) {
    skip <- 2
  } else {
    skip <- 0
  }
  
  # Read the TSV file
  if (grepl("\\.(tsv|annot|ko|tab)$", file_path)) {
    data <- read.delim(file_path, header = TRUE, sep = "\t", fill = TRUE, skip = skip)
    
    # Get the base name of the file path
    base_name <- basename(file_path)
    
    # Remove everything before the pattern "annotation"
    base_name <- sub(".*annotation", "annotation", base_name)
    
    # Remove everything before the pattern "annotation"
    base_name <- sub(".*assembly_evaluation", "assembly_evaluation", base_name)
    
    # Remove everything before the pattern "annotation"
    base_name <- sub(".*taxonomy_gtdbtk_summary.tsv", "taxonomy_summary", base_name)
    
    # Remove everything before the pattern "annotation"
    base_name <- sub(".*organism_name", "organism_name", base_name)
    
    # Remove the pattern ".tsv", ".annot", ".ko", ".tab"
    base_name <- sub("\\.(tsv|annot|ko|tab|txt)$", "", base_name)
    
    # Add a new column to the data frame that contains the name of the file
    data$source_file <- base_name
    
    # Rename the columns
    colnames(data)[colnames(data) == "Gene"] <- "gene"
    colnames(data)[colnames(data) == "Product"] <- "product"
    colnames(data)[colnames(data) == "Locus.Tag"] <- "locus_tag"
    colnames(data)[colnames(data) == "DbXrefs"] <- "database_references"
    colnames(data)[colnames(data) == "Start"] <- "start"
    colnames(data)[colnames(data) == "Stop"] <- "stop"
    colnames(data)[colnames(data) == "Strand"] <- "strand"
    colnames(data)[colnames(data) == "Type"] <- "type"
    colnames(data)[colnames(data) == "ftype"] <- "type"

    if ("annotation_bakta_microbeannotator" %in% data$source_file) {
      
      colnames(data)[colnames(data) == "query_id"] <- "locus_tag_bakta"
      colnames(data)[colnames(data) == "ko_product"] <- "KO_name_microbeannotator"
      colnames(data)[colnames(data) == "ko_number"] <- "K_number_microbeannotator"
      
      # Extracting EC numbers and splitting if multiple EC numbers exist
      data$EC_number_microbeannotator_1 <- gsub(".*\\[EC:|\\].*", "", data$KO_name)
      data$EC_number_microbeannotator_1 <- ifelse(data$EC_number_microbeannotator_1 == data$KO_name, NA, data$EC_number_microbeannotator_1)
      
      max_splits <- max(str_count(data$EC_number_microbeannotator_1, " ")[!is.na(str_count(data$EC_number_microbeannotator_1, " "))]) + 1
      
      for (i in 2:max_splits) {
        new_col <- paste0("EC_number_microbeannotator_", i)
        prev_col <- paste0("EC_number_microbeannotator_", i - 1)
        data[[new_col]] <- ifelse(grepl(" ", data[[prev_col]]), sub("^[^ ]+ ", "", data[[prev_col]]), NA)
        data[[prev_col]] <- ifelse(grepl(" ", data[[prev_col]]), sub(" .*$", "", data[[prev_col]]), data[[prev_col]])
      }
      
      # Removing everything after " [EC:"
      data$KO_name_microbeannotator_1 <- gsub(" \\[EC:.*", "", data$KO_name)
      
      max_splits <- max(str_count(data$KO_name_microbeannotator_1, " / ")[!is.na(str_count(data$KO_name_microbeannotator_1, " / "))]) + 1
      
      for (i in 2:max_splits) {
        new_col <- paste0("KO_name_microbeannotator_", i)
        prev_col <- paste0("KO_name_microbeannotator_", i - 1)
        
        # Splitting if multiple values exist separated by " / "
        data[[new_col]] <- ifelse(grepl(" / ", data[[prev_col]]), sub(".*? / ", "", data[[prev_col]]), NA)
        data[[prev_col]] <- ifelse(grepl(" / ", data[[prev_col]]), sub(" / .*$", "", data[[prev_col]]), data[[prev_col]])
      }
      
      #remove all columns only containing NA
      data$KO_name_microbeannotator <- NULL
      data <- data[, colSums(is.na(data)) != nrow(data)]
      
    }
    
    # Add a new column that is the difference between the "stop" and "start" columns if both of these columns exist and the "length_bp" column does not exist
    if ("start" %in% colnames(data) && "stop" %in% colnames(data) && !"length_bp" %in% colnames(data)) {
      data$length_bp <- data$stop - data$start
    }
    
    # Generate a unique name for the data frame
    df_name <- make.names(base_name, unique = TRUE)
  }
  
  
  # Read the txt file
  if (grepl(".txt$", file_path)) {
    data <- readLines(file_path, warn = FALSE)
    base_name <- basename(file_path)
    base_name <- sub(".*gtdbtk_organism_name.txt", "gtdbtk_organism_name", base_name)
    df_name <- make.names(base_name, unique = TRUE)
  }  
  
  # Read the ffn file
  if (grepl(".ffn$", file_path)) {
    data <- readDNAStringSet(file_path)
    base_name <- basename(file_path)
    base_name <- sub(".*annotation", "annotation", base_name)
    df_name <- make.names(base_name, unique = TRUE)
  }
  
  # Read the faa file
  if (grepl(".faa$", file_path)) {
    data <- readAAStringSet(file_path)
    base_name <- basename(file_path)
    base_name <- sub(".*annotation", "annotation", base_name)
    df_name <- make.names(base_name, unique = TRUE)
  }
  
  # Read the gbk file
  if (grepl("prokka.gbk$", file_path)) {
    data <- readGenBank(file_path, partial=TRUE)
    base_name <- basename(file_path)
    base_name <- sub(".*annotation", "annotation", base_name)
    df_name <- make.names(base_name, unique = TRUE)
    
    # Check if the data is a GenBankRecord object
    if (class(data) == "GenBankRecord") {
      
      # Extract the features
      features <- data
      
      # Extract the information
      gene_prokka <- features@cds$gene
      gene_product_prokka <- features@cds$product
      start_position <- features@cds@ranges@start
      end_position <- ((features@cds@ranges@start + features@cds@ranges@width) - 1)
      length_bp <- features@cds@ranges@width
      strand <- rep(features@cds@strand@values, times = features@cds@strand@lengths)
      type_prokka <- features@cds$type
      locus_tag_prokka <- features@cds$locus_tag
      #gene_id <- features@cds$gene_id
      EC_number_prokka <- as.character(features@cds$EC_number)
      COG_number_prokka <- as.character(features@cds$db_xref)
      COG_number_prokka <- gsub("^COG:", "", COG_number_prokka)
      K_number_prokka <- as.character(features@cds$db_xref)
      UniRef100_prokka <- as.character(features@cds$db_xref)
      UniRef90_prokka <- as.character(features@cds$db_xref)
      deviation_prokka <- as.character(features@cds$db_xref)
      amino_acid_code <- as.character(features@cds$translation)
      
      # Add the information to the summary data frame
      summary_df <- data.frame(gene_prokka, gene_product_prokka, start_position, end_position, length_bp, strand, type_prokka, locus_tag_prokka, EC_number_prokka, COG_number_prokka, K_number_prokka, UniRef100_prokka, UniRef90_prokka, deviation_prokka, amino_acid_code)
      summary_df$source_file_prokka <- basename(file_path)
      summary_df$K_number_prokka <- NA
      summary_df$UniRef100_prokka <- NA
      summary_df$UniRef90_prokka <- NA
      summary_df$deviation_prokka <- NA
      
      # Remove objects again
      gene_prokka <- NULL
      gene_product_prokka <- NULL
      start_position <- NULL
      end_position <- NULL
      length_bp <- NULL
      strand <- NULL
      type_prokka <- NULL
      # gene_id <- NULL
      EC_number_prokka <- NULL
      COG_number_prokka <- NULL
      amino_acid_code_prokka <- NULL
      
      # Remove columns that contain the pattern "group" in their names
      summary_df <- summary_df[ , !grepl("group", names(summary_df))]
      
      # Assign the summary data frame to a variable in the global environment
      assign(paste0("summary_", df_name), summary_df, envir = .GlobalEnv)
      
      # Set summary_df to NULL
      summary_df <- NULL
    }
  }
  
  # Read the gbff file
  if (grepl("bakta.gbff$", file_path)) {
    data <- readGenBank(file_path, partial=TRUE)
    base_name <- basename(file_path)
    base_name <- sub(".*annotation", "annotation", base_name)
    df_name <- make.names(base_name, unique = TRUE)
    
    # Check if the data is a GenBankRecord object
    if (class(data) == "GenBankRecord") {
      
      # Extract the features
      features <- data
      
      # Extract the information
      gene_bakta <- features@cds$gene
      gene_product_bakta <- features@cds$product
      start_position <- features@cds@ranges@start
      end_position <- ((features@cds@ranges@start + features@cds@ranges@width) - 1)
      length_bp <- features@cds@ranges@width
      strand <- rep(features@cds@strand@values, times = features@cds@strand@lengths)
      type_bakta <- features@cds$type
      locus_tag_bakta <- features@cds$locus_tag
      #gene_id <- features@cds$gene_id
      EC_number_bakta <- as.character(features@cds$EC_number)
      COG_number_bakta <- sapply(features@cds$note, paste, collapse = "; ") # Bakta has COG at different location than prokka
      K_number_bakta <- sapply(features@cds$note, paste, collapse = "; ") # Bakta has KEGG 
      UniRef100_bakta <- sapply(features@cds$note, paste, collapse = "; ") 
      UniRef90_bakta <- sapply(features@cds$note, paste, collapse = "; ")
      deviation_bakta <- sapply(features@cds$note, paste, collapse = "; ") 
      amino_acid_code <- as.character(features@cds$translation)
      
      # Add the information to the summary data frame
      summary_df <- data.frame(gene_bakta, gene_product_bakta, start_position, end_position, length_bp, strand, type_bakta, locus_tag_bakta, EC_number_bakta, COG_number_bakta, K_number_bakta, UniRef100_bakta, UniRef90_bakta, deviation_bakta, amino_acid_code)
      summary_df$source_file <- basename(file_path)
      
      # Append a ";" at the end of all the values in the COG_number column
      summary_df$COG_number_bakta <- ifelse(!is.na(summary_df$COG_number_bakta), 
                                            paste0(summary_df$COG_number_bakta, ";"), 
                                            NA)
      
      # Append a ";" at the end of all the values in the K_number column
      summary_df$K_number_bakta <- ifelse(!is.na(summary_df$K_number_bakta), 
                                             paste0(summary_df$K_number_bakta, ";"), 
                                             NA)
      
      # Append a ";" at the end of all the values in the UniRef90 column
      summary_df$UniRef90_bakta <- ifelse(!is.na(summary_df$UniRef90_bakta), 
                                          paste0(summary_df$UniRef90_bakta, ";"), 
                                          NA)
      
      # Append a ";" at the end of all the values in the UniRef100 column
      summary_df$UniRef100_bakta <- ifelse(!is.na(summary_df$UniRef100_bakta), 
                                           paste0(summary_df$UniRef100_bakta, ";"), 
                                           NA)
      
      # Remove the unnecessary characters from the COG_number column
      summary_df$COG_number_bakta <- ifelse(grepl("COG:", summary_df$COG_number_bakta), 
                                            sub(".*?COG:(.*?);.*", "\\1", summary_df$COG_number_bakta), 
                                            NA) # Replace the strings that do not contain "COG:" with NA
      
      # Remove the unnecessary characters from the K_number column
      summary_df$K_number_bakta <- ifelse(grepl("KEGG:", summary_df$K_number_bakta), 
                                             sub(".*?KEGG:(.*?);.*", "\\1", summary_df$K_number_bakta), 
                                             NA) # Replace the strings that do not contain "KEGG:" with NA
      
      # Remove the unnecessary characters from the UniRef90 column
      summary_df$UniRef90_bakta <- ifelse(grepl("UniRef:UniRef90", summary_df$UniRef90_bakta), 
                                          sub(".*?UniRef:UniRef90(.*?);.*", "\\1", summary_df$UniRef90_bakta), 
                                          NA) # Replace the strings that do not contain "UniRef:" with NA
      summary_df$UniRef90_bakta <- ifelse(!is.na(summary_df$UniRef90_bakta), 
                                          paste0("UniRef90", summary_df$UniRef90_bakta),
                                          NA)
      
      # Remove the unnecessary characters from the UniRef100 column
      summary_df$UniRef100_bakta <- ifelse(grepl("UniRef:UniRef100", summary_df$UniRef100_bakta), 
                                           sub(".*?UniRef:UniRef100(.*?);.*", "\\1", summary_df$UniRef100_bakta), 
                                           NA) # Replace the strings that do not contain "UniRef:" with NA
      summary_df$UniRef100_bakta <- ifelse(!is.na(summary_df$UniRef100_bakta), 
                                           paste0("UniRef100", summary_df$UniRef100_bakta),
                                           NA)
      
      # Append a ";" at the end of all the values in the deviation column
      summary_df$deviation_bakta <- ifelse(!is.na(summary_df$deviation_bakta), 
                                           paste0(summary_df$deviation_bakta, ";"), 
                                           NA)
      
      # Remove the unnecessary characters from the deviation column
      
      # Define the patterns
      patterns <- c("Nonsense", "Frameshift", "Internal")
      
      # Append a ";" at the end of all the values in the deviation column
      summary_df$deviation_bakta <- ifelse(!is.na(summary_df$deviation_bakta), 
                                           paste0(summary_df$deviation_bakta, ";"), 
                                           NA)
      
      # Check if any of the patterns is present in each string
      pattern_present <- apply(sapply(patterns, grepl, summary_df$deviation_bakta), 1, any)
      
      # Remove the unnecessary characters from the deviation column
      summary_df$deviation_bakta <- ifelse(pattern_present, 
                                           str_extract(summary_df$deviation_bakta, paste0("(" , paste(patterns, collapse = "|"), ")[^;]*;")), 
                                           NA) # Replace the strings that do not contain any of the patterns with NA
      
      # Remove the trailing semicolon
      summary_df$deviation_bakta <- sub(";$", "", summary_df$deviation_bakta)
      
      
      # Remove objects again
      gene_bakta <- NULL
      gene_product_bakta <- NULL
      start_position <- NULL
      end_position <- NULL
      length_bp <- NULL
      strand <- NULL
      type_bakta <- NULL
      # gene_id <- NULL
      EC_number_bakta <- NULL
      COG_number_bakta <- NULL
      amino_acid_code <- NULL
      
      # Remove columns that contain the pattern "group" in their names
      summary_df <- summary_df[ , !grepl("group", names(summary_df))]
      
      # Assign the summary data frame to a variable in the global environment
      assign(paste0("summary_", df_name), summary_df, envir = .GlobalEnv)
      
      # Set summary_df to NULL
      summary_df <- NULL
    }
  }
  
  # Assign the data frame to a variable in the global environment
  assign(df_name, data, envir = .GlobalEnv)
}

save.image(file = paste0(R_analysis_dir, "/", sample_name, "_BIOMES_WGseq.RData")) # just for development mode to have intermediate files


# -----
## Add the nucleotide sequences to the summary files from the .ffn files

# Get a list of all objects in the global environment
objects <- ls()

## For Prokka results
# Find the object that contains the patterns "prokka" and ".ffn" in its name
object_name <- objects[grepl("prokka", objects) & grepl(".ffn", objects)]

# Check if there is exactly one object that matches the patterns
if (length(object_name) == 1) {
  # Save the object as prokka_nucleotide_seq
  assign("prokka_nucleotide_seq", get(object_name), envir = .GlobalEnv)
} else {
  # Print a warning message if there is not exactly one object that matches the patterns
  print("Warning: There is not exactly one object that contains the patterns 'prokka' and '.ffn' in its name.")
}

# transform the DNAStringSet into a dataframe
prokka_nucleotide_seq <- as.data.frame(prokka_nucleotide_seq)

# Make row names the first column
prokka_nucleotide_seq <- data.frame(RowNames = rownames(prokka_nucleotide_seq), prokka_nucleotide_seq)

# Remove the row names
rownames(prokka_nucleotide_seq) <- NULL

# Rename the column names
names(prokka_nucleotide_seq) <- c("locus_tag_prokka", "nucleotide_code")

# Remove everything after the locus tag
prokka_nucleotide_seq$locus_tag_prokka <- sub(" .*", "", prokka_nucleotide_seq$locus_tag_prokka)

# Merge the nucleotide sequence to the summary file of prokka
summary_annotation_prokka.gbk <- merge(summary_annotation_prokka.gbk, prokka_nucleotide_seq, by = "locus_tag_prokka" , all = TRUE)

## for Bakta results
# Find the object that contains the patterns "bakta" and ".ffn" in its name
object_name <- objects[grepl("bakta", objects) & grepl(".ffn", objects)]

# Check if there is exactly one object that matches the patterns
if (length(object_name) == 1) {
  # Save the object as bakta_nucleotide_seq
  assign("bakta_nucleotide_seq", get(object_name), envir = .GlobalEnv)
} else {
  # Print a warning message if there is not exactly one object that matches the patterns
  print("Warning: There is not exactly one object that contains the patterns 'bakta' and '.ffn' in its name.")
}

# transform the DNAStringSet into a dataframe
bakta_nucleotide_seq <- as.data.frame(bakta_nucleotide_seq)

# Make row names the first column
bakta_nucleotide_seq <- data.frame(RowNames = rownames(bakta_nucleotide_seq), bakta_nucleotide_seq)

# Remove the row names
rownames(bakta_nucleotide_seq) <- NULL

# Rename the column names
names(bakta_nucleotide_seq) <- c("locus_tag_bakta", "nucleotide_code")

# Remove everything after the locus tag
bakta_nucleotide_seq$locus_tag_bakta <- sub(" .*", "", bakta_nucleotide_seq$locus_tag_bakta)

# Merge the nucleotide sequence to the summary file of bakta
summary_annotation_bakta.gbff <- merge(summary_annotation_bakta.gbff, bakta_nucleotide_seq, by = "locus_tag_bakta" , all = TRUE)


# -----
## Condense the data to summary files

# Define the columns to merge by to create a high fidelity data frame (overlap between bakta and prokka)
columns <- c("start_position", "end_position", "length_bp", "strand", "amino_acid_code", "nucleotide_code")

# Select only the columns of interest from each data frame
df1_selected <- summary_annotation_prokka.gbk[ , columns]
df2_selected <- summary_annotation_bakta.gbff[ , columns]

# Merge the selected columns from the data frames
high_fidelity_annotation <- merge(df1_selected, df2_selected, by = columns, all = FALSE)
high_fidelity_annotation <- unique(high_fidelity_annotation)

df1_selected <- summary_annotation_prokka.gbk
df2_selected <- summary_annotation_bakta.gbff

# Create the extended summary dataframe that is limited to high-fidelity annotations
annotation_summary_extended <- high_fidelity_annotation |>
  dplyr::left_join(df1_selected, by = c("start_position", "end_position", "length_bp", "strand", "amino_acid_code", "nucleotide_code"), relationship = "many-to-many") |>
  dplyr::left_join(df2_selected, by = c("start_position", "end_position", "length_bp", "strand", "amino_acid_code", "nucleotide_code"), relationship = "many-to-many") |>
  dplyr::left_join(
    dplyr::select(annotation_bakta_microbeannotator, locus_tag_bakta, dplyr::starts_with("KO_"), dplyr::starts_with("K_"), dplyr::starts_with("EC_")), 
    by = c("locus_tag_bakta"), relationship = "many-to-many") |>
  unique() |>
  dplyr::mutate(
    gene_name_dist = stringdist(tolower(gene_bakta), tolower(gene_prokka), method = "jw"),
    gene_name_consensus = ifelse(is.na(gene_bakta) & !is.na(gene_prokka), gene_prokka,
                            ifelse(!is.na(gene_bakta) & is.na(gene_prokka), gene_bakta,
                                   ifelse(!is.na(gene_bakta) & !is.na(gene_prokka) & (gene_bakta == gene_prokka), gene_bakta,
                                          ifelse(!is.na(gene_bakta) & !is.na(gene_prokka) & (gene_bakta != gene_prokka) & gene_name_dist < 0.15, gene_bakta,
                                                 ifelse(!is.na(gene_bakta) & !is.na(gene_prokka) & (gene_bakta != gene_prokka) & gene_name_dist > 0.15, gene_bakta,
                                                        NA
                                                 ))))),
    
    gene_name_consensus_confidence = ifelse(!is.na(gene_bakta) & is.na(gene_prokka), "only_one_result",
                                            ifelse(is.na(gene_bakta) & !is.na(gene_prokka), "only_one_result",
                                                   ifelse(gene_name_dist == 0, "complete_consensus",
                                                          ifelse(gene_name_dist > 0 & gene_name_dist <= 0.15, "high_probability_consensus",
                                                                 ifelse(gene_name_dist > 0.15 & gene_name_dist <= 0.2, "medium_probability_consensus",
                                                                        ifelse(gene_name_dist > 0.2 & gene_name_dist <= 0.45, "low_probability_consensus", 
                                                                               ifelse(gene_name_dist > 0.45, "no_consensus",
                                                                                      NA
                                                     )))))))
  ) |>
  dplyr::mutate(
    gene_product_dist = stringdist(tolower(gene_product_bakta), tolower(gene_product_prokka), method = "osa"),
    
    gene_product_consensus = ifelse(is.na(gene_product_bakta) & !is.na(gene_product_prokka), gene_product_prokka,
                                    ifelse(!is.na(gene_product_bakta) & is.na(gene_product_prokka), gene_product_bakta,
                                           ifelse(!is.na(gene_product_bakta) & !is.na(gene_product_prokka) & (gene_product_bakta == gene_product_prokka), gene_product_bakta,
                                                  ifelse(gene_product_bakta != "hypothetical protein" & !is.na(gene_product_bakta) & gene_product_prokka == "hypothetical protein", gene_product_bakta,
                                                         ifelse(gene_product_prokka != "hypothetical protein" & !is.na(gene_product_prokka) & gene_product_bakta == "hypothetical protein", gene_product_bakta,
                                                                ifelse(!is.na(gene_product_bakta) & !is.na(gene_product_prokka) & (gene_product_bakta != gene_product_prokka) & gene_product_dist <= 10, gene_product_bakta,
                                                                       ifelse(!is.na(gene_product_bakta) & !is.na(gene_product_prokka) & (gene_product_bakta != gene_product_prokka) & gene_product_dist > 10, gene_product_bakta,
                                                                              NA
                                                 ))))))),
    
    gene_product_consensus_confidence = ifelse(!is.na(gene_product_bakta) & is.na(gene_product_prokka), "only_one_result",
                                               ifelse(is.na(gene_product_bakta) & !is.na(gene_product_prokka), "only_one_result",
                                                      ifelse(gene_product_dist == 0, "complete_consensus",
                                                             ifelse(gene_product_dist > 0 & gene_product_dist <= 10, "high_probability_consensus",
                                                                    ifelse(gene_product_dist > 10 & gene_product_dist < 20, "medium_probability_consensus",
                                                                           ifelse(gene_product_dist >=20 & gene_product_dist <25, "low_probability_consensus",
                                                                                  ifelse(gene_product_dist >=25, "no_consensus",
                                                                                         NA
                                                          )))))))
    
  ) |>
  dplyr::select("gene_name_consensus", "gene_name_consensus_confidence", "gene_bakta", "gene_prokka", "gene_product_consensus", "gene_product_consensus_confidence", "gene_product_bakta", "gene_product_prokka", "start_position", "end_position", "length_bp", "strand", "amino_acid_code", "nucleotide_code", dplyr::starts_with("EC"), dplyr::starts_with("COG"), dplyr::starts_with("K_"), dplyr::starts_with("KO_"), dplyr::starts_with("UniRef100"), dplyr::starts_with("UniRef90"), dplyr::starts_with("deviation"), dplyr::starts_with("locus_tag"))

## Create the compact summary list
# Define a function to convert column values to a data frame
column_to_list <- function(df, columns, new_col_name) {
  # Create an empty list to store data frames
  df_values <- vector("list", nrow(df))
  
  # Iterate over each row
  for (i in 1:nrow(df)) {
    # Extract values from specified columns and store in a data frame
    row_values <- df[i, columns]
    # Extract values from additional columns
    additional_values <- df[i, c("locus_tag_prokka", "locus_tag_bakta", "nucleotide_code")]
    # Combine the two into a single data frame
    df_values[[i]] <- setNames(data.frame(unlist(row_values), additional_values), c(new_col_name, names(additional_values)))
  }
  
  return(df_values)
}


# Create the annotation_summary frame as copy of annotation_summary_extended
annotation_summary <- annotation_summary_extended

# Create a columns where genes have an unique appendix if they appear multiple times (important for the shiny app)
# Assuming your data frame is called df and the column with character entries is called gene_names
duplicates <- duplicated(annotation_summary$gene_name_consensus) | duplicated(annotation_summary$gene_name_consensus, fromLast = TRUE)
annotation_summary$gene_name_consensus_version_appendix[duplicates] <- ave(annotation_summary$gene_name_consensus[duplicates], annotation_summary$gene_name_consensus[duplicates], FUN = function(x) paste0(x, "_v", seq_along(x)))
annotation_summary$gene_name_consensus_version_appendix <- ifelse(is.na(annotation_summary$gene_name_consensus_version_appendix), annotation_summary$gene_name_consensus, annotation_summary$gene_name_consensus_version_appendix)

# Columns to convert to lists
EC_columns <- grep("^EC_number", names(annotation_summary), value = TRUE)
COG_columns <- grep("^COG_number", names(annotation_summary), value = TRUE)
K_columns <- grep("^K_number", names(annotation_summary), value = TRUE)
gene_name_columns <- grep("^(gene_name|gene_bakta|gene_prokka)", names(annotation_summary), value = TRUE)
gene_product_columns <- grep("^gene_product", names(annotation_summary), value = TRUE)

# Create new columns with lists of values
annotation_summary$EC_list <- column_to_list(annotation_summary, EC_columns, "EC_number")
annotation_summary$COG_list <- column_to_list(annotation_summary, COG_columns, "COG_number")
annotation_summary$K_list <- column_to_list(annotation_summary, K_columns, "K_number")
annotation_summary$gene_name_list <- column_to_list(annotation_summary, gene_name_columns, "gene_name")
annotation_summary$gene_product_list <- column_to_list(annotation_summary, gene_product_columns, "gene_product")

annotation_summary <- dplyr::select(annotation_summary, "gene_name_consensus", "gene_name_consensus_confidence", "gene_product_consensus", "gene_product_consensus_confidence", "start_position", "end_position", "length_bp", "strand", "amino_acid_code", "nucleotide_code", dplyr::ends_with("_list"), dplyr::starts_with("locus_tag"), gene_name_consensus_version_appendix)

# Rename the metabolic summary created by microbeannotator
metabolic_summary_microbeannotator <- metabolic_summary__module_completeness
metabolic_summary_microbeannotator$comment <- "KEGG_analysis_result_by_MicrobeAnnotator"
colnames(metabolic_summary_microbeannotator)[colnames(metabolic_summary_microbeannotator) == "module"] <- "M_number"
colnames(metabolic_summary_microbeannotator)[colnames(metabolic_summary_microbeannotator) == "name"] <- "description"
colnames(metabolic_summary_microbeannotator)[colnames(metabolic_summary_microbeannotator) == "pathway.group"] <- "pathway_group"
colnames(metabolic_summary_microbeannotator)[4] <- "module_completeness"


# -----
## Create the taxonomy overviews
taxonomy_summary <- as.data.frame(lapply(taxonomy_summary, as.character))
taxonomy_summary <- tidyr::pivot_longer(taxonomy_summary, colnames(taxonomy_summary[, 1:ncol(taxonomy_summary)]))

# Add the info of the chosen file
taxonomy_summary <- rbind(c("chosen_assembly", chosen_assembly), taxonomy_summary)
taxonomy_summary$name[taxonomy_summary$name == "user_genome"] <- "input_file_name"
taxonomy_summary[taxonomy_summary == "N/A"] <- NA

# Explanatory labels for each row
explanatory_labels <- c(
  "The assembly that was chosen for downstream analysis",
  "The name of the file that was used",
  "Classification result",
  "Reference used for FastANI calculation",
  "Radius used for FastANI calculation",
  "Taxonomy based on FastANI",
  "Average Nucleotide Identity (ANI) obtained from FastANI",
  "Average Fraction of Orthologous Genes (AF) obtained from FastANI",
  "Reference with closest placement",
  "Radius of closest placement",
  "Taxonomy of closest placement",
  "Average Nucleotide Identity (ANI) compared to closest placement",
  "Average Fraction of Orthologous Genes (AF) compared to closest placement",
  "Taxonomy obtained from pplacer",
  "Classification method used",
  "Additional notes",
  "Other related references: genome_id.species_name.radius.ANI.AF.",
  "MSA Percent",
  "Translation Table",
  "Red result",
  "Warnings",
  "Source file"
)

# Add explanatory column to the taxonomy_summary data frame
taxonomy_summary <- cbind(taxonomy_summary, Explanation = explanatory_labels)

# Remove unnecessary rows
taxonomy_summary <- taxonomy_summary[taxonomy_summary$name != "source_file",]
taxonomy_summary <- na.omit(taxonomy_summary)

# Create other related taxa table
other_related_organisms <- taxonomy_summary[grepl("other_related", taxonomy_summary$name),]

# Split the 'value' column by semicolons to create new rows
other_related_organisms <- strsplit(other_related_organisms$value, ";\\s*")[[1]]

# Split each row by commas to create separate columns
other_related_organisms <- strsplit(other_related_organisms, ",\\s*")

# Convert the result into a data frame
other_related_organisms <- data.frame(do.call(rbind, other_related_organisms))

# Rename columns
colnames(other_related_organisms) <- c("genome_id", "species_name", "radius", "ANI", "AF")

# Finalise taxonomy summaries
colnames(taxonomy_summary) <- c("data_type", "result", "explanation")
taxonomy_summary_short <- taxonomy_summary |>
  dplyr::filter(data_type == "classification" |
                  data_type == "fastani_ani" |
                  data_type == "fastani_af" |
                  data_type == "closest_placement_reference" |
                  data_type == "closest_placement_taxonomy" |
                  data_type == "classification_method"
                ) 

taxonomy_summary <- taxonomy_summary 


# -----
## Create the assembly overviews
assembly_evaluation <- as.data.frame(lapply(assembly_evaluation, as.character))

# Extract the assembler info
assembly_evaluation$Assembly <- sub(".*?\\Qassembly_\\E", "", assembly_evaluation$Assembly, perl = TRUE)
assembly_evaluation$Assembly <- paste(assembly_evaluation$Assembly, "_assembler", sep = "")

# Assuming your data frame is called assembly_evaluation
assembly_evaluation <- assembly_evaluation %>%
  pivot_longer(cols = -Assembly, names_to = "data_type", values_to = "result") %>%
  pivot_wider(names_from = Assembly, values_from = result)

# Define a mapping of original column names to new, more descriptive names
column_name_mapping <- c(
  "X..contigs.....0.bp."       = "Num contigs (>= 0 bp)",
  "X..contigs.....1000.bp."    = "Num contigs (>= 1000 bp)",
  "X..contigs.....5000.bp."    = "Num contigs (>= 5000 bp)",
  "X..contigs.....10000.bp."   = "Num contigs (>= 10000 bp)",
  "X..contigs.....25000.bp."   = "Num contigs (>= 25000 bp)",
  "X..contigs.....50000.bp."   = "Num contigs (>= 50000 bp)",
  "Total.length.....0.bp."     = "Total length (>= 0 bp)",
  "Total.length.....1000.bp."  = "Total length (>= 1000 bp)",
  "Total.length.....5000.bp."  = "Total length (>= 5000 bp)",
  "Total.length.....10000.bp." = "Total length (>= 10000 bp)",
  "Total.length.....25000.bp." = "Total length (>= 25000 bp)",
  "Total.length.....50000.bp." = "Total length (>= 50000 bp)",
  "X..contigs"                 = "Total number of contigs",
  "Largest.contig"             = "Largest contig",
  "Total.length"               = "Total assembly length",
  "GC...."                     = "GC percentage",
  "N50"                        = "N50",
  "N75"                        = "N75",
  "L50"                        = "L50",
  "L75"                        = "L75",
  "X..N.s.per.100.kbp"         = "N's per 100 kbp",
  "source_file"                = "Source file"
)

# Replace the original values with the new descriptive values
assembly_evaluation$data_type <- column_name_mapping[assembly_evaluation$data_type]

# Remove unnecessary row
assembly_evaluation <- assembly_evaluation[assembly_evaluation$data_type != "Source file",]

# Add chosen file row
pattern <- gsub("_assembly", "", chosen_assembly)

# Create a new row
new_row <- c("Selected assembly for downstream analysis",
             ifelse(grepl(pattern, names(assembly_evaluation)[-1]), "Yes", "No"))

# Bind the new row to the dataframe
assembly_evaluation <- rbind(assembly_evaluation, new_row)

# Explanatory labels for each row
explanatory_labels <- c(
  "Number of contigs with length greater than or equal to 0 base pairs",
  "Number of contigs with length greater than or equal to 1000 base pairs",
  "Number of contigs with length greater than or equal to 5000 base pairs",
  "Number of contigs with length greater than or equal to 10000 base pairs",
  "Number of contigs with length greater than or equal to 25000 base pairs",
  "Number of contigs with length greater than or equal to 50000 base pairs",
  "Total length of contigs with length greater than or equal to 0 base pairs",
  "Total length of contigs with length greater than or equal to 1000 base pairs",
  "Total length of contigs with length greater than or equal to 5000 base pairs",
  "Total length of contigs with length greater than or equal to 10000 base pairs",
  "Total length of contigs with length greater than or equal to 25000 base pairs",
  "Total length of contigs with length greater than or equal to 50000 base pairs",
  "Total number of contigs",
  "Length of the largest contig",
  "Total assembly length",
  "GC percentage",
  "N50: Contig length such that 50% of the assembly length is contained in contigs of this size or larger",
  "N75: Contig length such that 75% of the assembly length is contained in contigs of this size or larger",
  "L50: Number of contigs whose lengths sum to reach or exceed 50% of the total assembly length",
  "L75: Number of contigs whose lengths sum to reach or exceed 75% of the total assembly length",
  "N's per 100 kbp",
  "Selected assembly for analyses: taxonomy classification, genome annotation and metabolic function"
)

# Add explanatory column to the taxonomy_summary data frame
assembly_evaluation <- cbind(assembly_evaluation, Explanation = explanatory_labels)


# Finalise the assembly summaries 
assembly_evaluation <- dplyr::select(assembly_evaluation, data_type, lja_assembler, spades_assembler, canu_assembler, cisa_assembler, Explanation)

assembly_summary <- assembly_evaluation

assembly_summary_short <- assembly_evaluation |>
  dplyr::filter(
    data_type == "Total assembly length" |
      data_type == "Total number of contigs" |
      data_type == "Largest contig" |
      data_type == "N50" |
      data_type == "N75" |
      data_type == "L50" |
      data_type == "L75" |
      data_type == "Selected assembly for downstream analysis"
  )


# -----
## Remove unecessary objects
rm(all_dir,
   assembly_evaluation,
   annotation_bakta,
   annotation_bakta.hypotheticals,
   annotation_prokka,
   annotation_bakta_microbeannotator,
   bakta_nucleotide_seq,
   prokka_nucleotide_seq,
   data,
   df1_selected,
   df2_selected,
   features,
   high_fidelity_annotation,
   summary_annotation_bakta.gbff,
   summary_annotation_prokka.gbk,
   amino_acid_code,
   amino_acid_code_prokka,
   base_name,
   check_and_install,
   COG_columns,
   COG_number_bakta,
   COG_number_prokka,
   chosen_assembly_paths,
   columns,
   column_name_mapping,
   deviation_bakta,
   deviation_prokka,
   df_name,
   duplicates,
   EC_columns,
   EC_number_bakta,
   EC_number_prokka,
   end_position,
   explanatory_labels,
   file_path,
   file_paths,
   file_patterns,
   gene_bakta,
   gene_name_columns,
   gene_product_bakta,
   gene_product_columns,
   gene_product_prokka,
   gene_prokka,
   i,
   K_columns,
   K_number_bakta,
   K_number_prokka,
   length_bp,
   locus_tag_bakta,
   locus_tag_prokka,
   main_result_dir,
   max_splits,
   new_col,
   new_row,
   object_name,
   objects,
   main_results_dir_abs,
   paths,
   pattern,
   patterns,
   pattern_present,
   prev_col,
   skip,
   start_position,
   strand,
   summary_df,
   type_bakta,
   type_prokka,
   UniRef100_bakta,
   UniRef100_prokka,
   UniRef90_bakta,
   UniRef90_prokka,
   column_to_list,
   metabolic_summary__module_completeness,
   R_analysis_dir
)

# Save the R environment at a specified location
save.image(file = paste0(R_analysis_dir, "/", sample_name, "_BIOMES_WGseq.RData"))