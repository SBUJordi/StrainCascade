#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# genome_annotation_integration.R

set.seed(42)

# Load necessary libraries
suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(stringdist)
  library(optparse)
})

# Parse command-line arguments
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, 
              help="Path to the input qs file (annotation_results_aggregated.qs)", metavar="FILE"),
  make_option(c("--output_dir"), type="character", default=".", 
              help="Output directory [default= %default]", metavar="DIRECTORY"),
  make_option(c("--version"), type="character", default=NULL,
              help="StrainCascade version", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$input_file) || is.null(opt$version)) {
  print_help(opt_parser)
  stop("Input file and version must be supplied.", call.=FALSE)
}

# Set the random seed
set.seed(77)

# Import the aggregated annotation results
annotation_results_integrated <- qread(opt$input_file)

# Exit if no files were imported
if (is.null(annotation_results_integrated)) {
  stop("ERROR: No files available for import. Exiting script.")
}

# Process data to assign a new locus tag
random_number <- sprintf("%03d", sample(0:999, 1))           # Random 3-digit number
random_letter <- sample(c(letters, LETTERS), 1)              # Random letter (upper or lower case)
current_date <- format(Sys.Date(), "%y%m%d")                 # Current date in 'yymmdd' format
count <- sprintf("%05d", seq(1, nrow(annotation_results_integrated))) # Sequential count

# Create the 'locus_tag_SC' column
annotation_results_integrated <- annotation_results_integrated %>%
  mutate(locus_tag_SC = paste0("SC", random_number, random_letter, current_date, "_", count))

rm(
  count,
  current_date,
  random_letter,
  random_number
)

# Functions to handle EC number processing

# Function to check if one EC number is less complete than another
is_less_complete <- function(ec1, ec2) {
  parts1 <- strsplit(ec1, "\\.")[[1]]
  parts2 <- strsplit(ec2, "\\.")[[1]]
  
  for (i in seq_along(parts1)) {
    if (parts1[i] == "-" && parts2[i] != "-") return(TRUE)
    if (parts1[i] != parts2[i]) return(FALSE)
  }
  
  return(FALSE)
}

# Function to filter out less complete EC numbers within each group
filter_less_complete <- function(group) {
  EC_numbers <- group$EC_numbers
  to_remove <- logical(length(EC_numbers))
  
  for (i in seq_along(EC_numbers)) {
    for (j in seq_along(EC_numbers)) {
      if (i != j && is_less_complete(EC_numbers[i], EC_numbers[j])) {
        to_remove[i] <- TRUE
      }
    }
  }
  
  group[!to_remove, ]
}

# Function to calculate the score based on the completeness of the EC number
calculate_EC_score <- function(EC_number) {
  parts <- strsplit(EC_number, "\\.")[[1]]
  score <- sum(parts != "-")
  return(score)
}

# Function to check if EC numbers match or are partial matches
is_partial_match <- function(ec1, ec2) {
  parts1 <- strsplit(ec1, "\\.")[[1]]
  parts2 <- strsplit(ec2, "\\.")[[1]]
  
  for (i in seq_along(parts1)) {
    if (parts1[i] != parts2[i] && parts1[i] != "-" && parts2[i] != "-") {
      return(FALSE)
    }
  }
  
  return(TRUE)
}

# Main processing of EC numbers
# Select relevant columns
EC_number_raw_data <- annotation_results_integrated %>%
  select(locus_tag_SC, starts_with("EC"))

# Reshape and clean the EC numbers
cleaned_EC_numbers <- EC_number_raw_data %>%
  pivot_longer(cols = starts_with("EC"), names_to = "EC_column", values_to = "EC_numbers") %>%
  separate_rows(EC_numbers, sep = "; ") %>%
  filter(!is.na(EC_numbers))

# Group by locus_tag and EC_column, then filter less complete EC numbers
filtered_EC_data <- cleaned_EC_numbers %>%
  group_by(locus_tag_SC, EC_column) %>%
  group_modify(~ filter_less_complete(.x)) %>%
  distinct() %>%
  ungroup() %>%
  mutate(EC_score = sapply(EC_numbers, calculate_EC_score))

# Group by locus_tag_SC and process EC numbers
grouped_EC_results <- filtered_EC_data %>%
  group_by(locus_tag_SC) %>%
  summarise(
    EC_numbers = list(EC_numbers),
    EC_columns = list(EC_column),
    EC_scores = list(EC_score),
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(
    grouped_EC = list({
      EC_numbers <- unlist(EC_numbers)
      EC_columns <- unlist(EC_columns)
      scores <- unlist(EC_scores)
      
      groups <- list()
      assigned <- logical(length(EC_numbers))
      
      for (i in seq_along(EC_numbers)) {
        if (!assigned[i]) {
          group <- c(i)
          for (j in seq_along(EC_numbers)) {
            if (i != j && !assigned[j] && is_partial_match(EC_numbers[i], EC_numbers[j])) {
              group <- c(group, j)
            }
          }
          groups[[length(groups) + 1]] <- group
          assigned[group] <- TRUE
        }
      }
      
      lapply(groups, function(group) {
        most_complete_EC <- EC_numbers[group][which.max(scores[group])]
        data.frame(
          EC_number_SC = most_complete_EC,
          EC_number_SC_score = sum(scores[group]),
          EC_columns = paste(unique(EC_columns[group]), collapse = ", ")
        )
      })
    })
  ) %>%
  unnest(grouped_EC)

# Add locus_tag_SC to the final grouped data by looping through each element in grouped_EC_results$grouped_EC
for (i in seq_along(grouped_EC_results$grouped_EC)) {
  # Add the new column locus_tag_SC with the value of result$locus_tag_SC[1]
  grouped_EC_results$grouped_EC[[i]] <- grouped_EC_results$grouped_EC[[i]] %>%
    mutate(locus_tag_SC = grouped_EC_results$locus_tag_SC[i])
}

# Final result combining the grouped EC data
EC_number_scores <- bind_rows(grouped_EC_results$grouped_EC) %>%
  select(locus_tag_SC, EC_number_SC, EC_number_SC_score) %>%
  distinct() %>%
  arrange(locus_tag_SC, desc(EC_number_SC_score))

# Summarize the final results
processed_EC_numbers <- EC_number_scores %>%
  group_by(locus_tag_SC) %>%
  summarise(
    EC_number_SC_best = paste0(EC_number_SC[EC_number_SC_score == max(EC_number_SC_score)], collapse = "; "),
    EC_number_SC_all = paste0(EC_number_SC[order(-EC_number_SC_score)], collapse = "; ")
  ) %>%
  ungroup()

# Remove the unprocessed EC numbers and add the processed EC numbers
annotation_results_integrated <- annotation_results_integrated |>
  select(-starts_with("EC_number")) |>
  left_join(processed_EC_numbers, by = "locus_tag_SC")

# Add documentation to the variables
base::comment(annotation_results_integrated$EC_number_SC_best) <- 
  "The best (highest scoring) EC number(s) for this locus tag. Multiple numbers are separated by semicolons if they have the same highest score."

base::comment(annotation_results_integrated$EC_number_SC_all) <- 
  "All EC numbers (from all available tools) for this locus tag, ordered by decreasing score. Numbers are separated by semicolons."

base::comment(EC_number_scores$EC_number_SC_score) <- 
  "Cumulative score for this EC number, calculated by summing individual scores from all contributing tools. Each EC number is scored based on its completeness (e.g., 1.1.1.1 scores higher than 1.1.1.-) and the number of tools that predicted it. The score for each EC number prediction is calculated using the calculate_EC_score() function, which assigns higher scores to more complete EC numbers. When multiple tools predict the same or partially matching EC numbers, their scores are combined. This cumulative score represents the confidence and consensus across different prediction tools."

# Make the process more robust
if (!"EC_number_SC_best" %in% names(annotation_results_integrated)) {
  annotation_results_integrated$EC_number_SC_best <- NA_character_
  warning("EC_number_SC_best not found in annotation_results_integrated. Added as NA.")
}

if (!"EC_number_SC_all" %in% names(annotation_results_integrated)) {
  annotation_results_integrated$EC_number_SC_all <- NA_character_
  warning("EC_number_SC_all not found in annotation_results_integrated. Added as NA.")
}

if (!"EC_number_SC_score" %in% names(EC_number_scores)) {
  EC_number_scores$EC_number_SC_score <- NA_real_
  warning("EC_number_SC_score not found in EC_number_scores. Added as NA.")
}

rm(
  cleaned_EC_numbers,
  processed_EC_numbers,
  filtered_EC_data,
  grouped_EC_results,
  calculate_EC_score,
  filter_less_complete,
  is_less_complete,
  is_partial_match,
  i
)

# processed_EC_numbers now contains the cleaned and processed EC numbers


# Main processing of COG numbers
available_columns <- c("locus_tag_SC", "COG_number_bakta", "COG_category_bakta", "COG_number_prokka")
available_columns <- available_columns[available_columns %in% names(annotation_results_integrated)]

if (length(available_columns) == 1 && available_columns == "locus_tag_SC") {
  message("Only 'locus_tag_SC' is available. No COG information available.")
} else {
  message("Selected columns: ", paste(available_columns, collapse = ", "))
}

COG_number_raw_data <- annotation_results_integrated %>% select(all_of(available_columns))

if ("COG_number_bakta" %in% available_columns && !"COG_number_prokka" %in% available_columns) {
  COG_number_raw_data <- COG_number_raw_data %>%
    dplyr::rename(COG_number_SC = COG_number_bakta) %>%
    mutate(COG_number_SC_score = ifelse(!is.na(COG_number_SC), 1, NA_real_))
  
  comment(COG_number_raw_data$COG_number_SC) <- "COG_number_SC: Derived from COG_number_bakta. Score indicates presence (1)."
  comment(COG_number_raw_data$COG_number_SC_score) <- "COG_number_SC_score: 1 if COG_number_SC is not NA."
  
} else if (!"COG_number_bakta" %in% available_columns && "COG_number_prokka" %in% available_columns) {
  COG_number_raw_data <- COG_number_raw_data %>%
    dplyr::rename(COG_number_SC = COG_number_prokka) %>%
    mutate(COG_number_SC_score = ifelse(!is.na(COG_number_SC), 1, NA_real_))
  
  comment(COG_number_raw_data$COG_number_SC) <- "COG_number_SC: Derived from COG_number_prokka. Score indicates presence (1)."
  comment(COG_number_raw_data$COG_number_SC_score) <- "COG_number_SC_score: 1 if COG_number_SC is not NA."
  
} else if ("COG_number_bakta" %in% available_columns && "COG_number_prokka" %in% available_columns) {
  COG_number_raw_data <- COG_number_raw_data %>%
    mutate(
      COG_number_SC = ifelse(!is.na(COG_number_bakta), COG_number_bakta, COG_number_prokka),
      COG_number_SC_score = case_when(
        !is.na(COG_number_bakta) & !is.na(COG_number_prokka) & (COG_number_bakta == COG_number_prokka) ~ 2,
        !is.na(COG_number_bakta) & !is.na(COG_number_prokka) & (COG_number_bakta != COG_number_prokka) ~ 0,
        !is.na(COG_number_SC) ~ 1,
        TRUE ~ NA_real_
      )
    )
  
  comment(COG_number_raw_data$COG_number_SC) <- "COG_number_SC: Derived from COG_number_bakta if available; otherwise, from COG_number_prokka."
  comment(COG_number_raw_data$COG_number_SC_score) <- "COG_number_SC_score: 2 if both bakta and prokka agree, 0 if they disagree, 1 if only one is available."
}

if ("COG_category_bakta" %in% available_columns) {
  COG_number_raw_data <- COG_number_raw_data %>%
    mutate(COG_category_bakta = ifelse(is.na(COG_number_bakta), NA, COG_category_bakta)) %>%
    dplyr::rename(COG_category_SC = COG_category_bakta)
  
  comment(COG_number_raw_data$COG_category_SC) <- "COG_category_SC: Derived from COG_category_bakta, set to NA if COG_number_bakta is NA."
}

if (exists("COG_number_raw_data") && nrow(COG_number_raw_data) > 0) {
  annotation_results_integrated <- annotation_results_integrated %>%
    select(-all_of(setdiff(available_columns, "locus_tag_SC"))) %>%
    left_join(COG_number_raw_data %>% select(-all_of(ends_with(c("bakta", "prokka", "SC_score")))), by = "locus_tag_SC")
  
  COG_number_scores <- COG_number_raw_data %>%
    select(contains("SC"), ends_with("_score"))
  
  # Attach comments to the final columns in annotation_results_integrated
  comment(annotation_results_integrated$COG_number_SC) <- comment(COG_number_raw_data$COG_number_SC)
  comment(annotation_results_integrated$COG_category_SC) <- comment(COG_number_raw_data$COG_category_SC)
  
  # Attach comments to the COG_number_scores dataframe
  comment(COG_number_scores$COG_number_SC) <- comment(COG_number_raw_data$COG_number_SC)
  comment(COG_number_scores$COG_number_SC_score) <- comment(COG_number_raw_data$COG_number_SC_score)
} else {
  message("No COG_number_raw_data generated. Skipping the join and score creation steps.")
}

rm(
  available_columns
)

# Final Result: annotation_results_integrated now contains the new COG columns and scores


# Main processing of K numbers

available_columns <- c("locus_tag_SC", "K_number_bakta", "K_number_microbeannotator")
available_columns <- available_columns[available_columns %in% names(annotation_results_integrated)]

if (length(available_columns) == 1 && available_columns == "locus_tag_SC") {
  message("Only 'locus_tag_SC' is available. No K number information available.")
} else {
  message("Selected columns: ", paste(available_columns, collapse = ", "))
}

K_number_raw_data <- annotation_results_integrated %>% select(all_of(available_columns))

if ("K_number_bakta" %in% available_columns && !"K_number_microbeannotator" %in% available_columns) {
  K_number_raw_data <- K_number_raw_data %>%
    dplyr::rename(K_number_SC = K_number_bakta) %>%
    mutate(K_number_SC_score = ifelse(!is.na(K_number_SC), 1, NA_real_))
  
  comment(K_number_raw_data$K_number_SC) <- "K_number_SC: Derived from K_number_bakta. Score indicates presence (1)."
  comment(K_number_raw_data$K_number_SC_score) <- "K_number_SC_score: 1 if K_number_SC is not NA."
  
} else if (!"K_number_bakta" %in% available_columns && "K_number_microbeannotator" %in% available_columns) {
  K_number_raw_data <- K_number_raw_data %>%
    dplyr::rename(K_number_SC = K_number_microbeannotator) %>%
    mutate(K_number_SC_score = ifelse(!is.na(K_number_SC), 1, NA_real_))
  
  comment(K_number_raw_data$K_number_SC) <- "K_number_SC: Derived from K_number_microbeannotator. Score indicates presence (1)."
  comment(K_number_raw_data$K_number_SC_score) <- "K_number_SC_score: 1 if K_number_SC is not NA."
  
} else if ("K_number_bakta" %in% available_columns && "K_number_microbeannotator" %in% available_columns) {
  K_number_raw_data <- K_number_raw_data %>%
    mutate(
      K_number_SC = ifelse(!is.na(K_number_bakta), K_number_bakta, K_number_microbeannotator),
      K_number_SC_score = case_when(
        !is.na(K_number_bakta) & !is.na(K_number_microbeannotator) & (K_number_bakta == K_number_microbeannotator) ~ 2,
        !is.na(K_number_bakta) & !is.na(K_number_microbeannotator) & (K_number_bakta != K_number_microbeannotator) ~ 0,
        !is.na(K_number_SC) ~ 1,
        TRUE ~ NA_real_
      )
    )
  
  comment(K_number_raw_data$K_number_SC) <- "K_number_SC: Derived from K_number_bakta if available; otherwise, from K_number_microbeannotator."
  comment(K_number_raw_data$K_number_SC_score) <- "K_number_SC_score: 2 if both bakta and microbeannotator agree, 0 if they disagree, 1 if only one is available."
}

if (exists("K_number_raw_data") && nrow(K_number_raw_data) > 0) {
  annotation_results_integrated <- annotation_results_integrated %>%
    select(-all_of(setdiff(available_columns, "locus_tag_SC"))) %>%
    left_join(K_number_raw_data %>% select(-all_of(ends_with(c("bakta", "microbeannotator", "SC_score")))), by = "locus_tag_SC")
  
  K_number_scores <- K_number_raw_data %>%
    select(contains("SC"), ends_with("_score"))
  
  # Attach comments to the final columns in annotation_results_integrated
  comment(annotation_results_integrated$K_number_SC) <- comment(K_number_raw_data$K_number_SC)
  
  # Attach comments to the K_number_scores dataframe
  comment(K_number_scores$K_number_SC) <- comment(K_number_raw_data$K_number_SC)
  comment(K_number_scores$K_number_SC_score) <- comment(K_number_raw_data$K_number_SC_score)
} else {
  message("No K_number_raw_data generated. Skipping the join and score creation steps.")
}

rm(available_columns, K_number_scores)

# Final Result: annotation_results_integrated now contains the new K number columns and scores


# Main processing of gene names and products
available_columns <- c("locus_tag_SC", "gene_bakta", "gene_product_bakta", "gene_prokka", "gene_product_prokka", "gene_product_microbeannotator")
available_columns <- available_columns[available_columns %in% names(annotation_results_integrated)]
message("Selected columns: ", paste(available_columns, collapse = ", "))
gene_data_raw <- annotation_results_integrated %>% select(all_of(available_columns))

# Block 1: Process gene names
gene_data_raw <- gene_data_raw %>%
  mutate(
    # Use Jaro-Winkler (jw) distance for gene names as they are typically short
    # JW is better for short strings as it gives more weight to characters that match at the beginning
    gene_name_dist = stringdist(tolower(gene_bakta), tolower(gene_prokka), method = "jw"),
    
    # Determine consensus gene name
    # Algorithm:
    # 1. If only one gene name is available, use that
    # 2. If both are available and identical, use that
    # 3. If both are available but different, prefer bakta's gene name
    # 4. If neither is available, use NA
    gene_SC = case_when(
      is.na(gene_bakta) & !is.na(gene_prokka) ~ gene_prokka,
      !is.na(gene_bakta) & is.na(gene_prokka) ~ gene_bakta,
      !is.na(gene_bakta) & !is.na(gene_prokka) & (gene_bakta == gene_prokka) ~ gene_bakta,
      !is.na(gene_bakta) & !is.na(gene_prokka) & (gene_bakta != gene_prokka) ~ gene_bakta, # Preference given to bakta
      TRUE ~ NA_character_
    ),
    
    # Determine confidence level of consensus
    # Algorithm:
    # 1. If only one result is available, mark as "only_one_result"
    # 2. If both are available and identical, mark as "complete_consensus"
    # 3. If both are available but different, use Jaro-Winkler distance to determine confidence level:
    #    - 0.00 to 0.15: high probability
    #    - 0.15 to 0.20: medium probability
    #    - 0.20 to 0.45: low probability
    #    - > 0.45: no consensus
    gene_SC_confidence = case_when(
      (!is.na(gene_bakta) & is.na(gene_prokka)) | (is.na(gene_bakta) & !is.na(gene_prokka)) ~ "only_one_result",
      gene_name_dist == 0 ~ "complete_consensus",
      gene_name_dist > 0 & gene_name_dist <= 0.15 ~ "high_probability_consensus",
      gene_name_dist > 0.15 & gene_name_dist <= 0.2 ~ "medium_probability_consensus",
      gene_name_dist > 0.2 & gene_name_dist <= 0.45 ~ "low_probability_consensus",
      gene_name_dist > 0.45 ~ "no_consensus",
      TRUE ~ NA_character_
    )
  )

# Block 2: Process gene products
gene_data_raw <- gene_data_raw %>%
  mutate(
    # Use Optimal String Alignment (osa) distance for gene products as they are typically longer
    # OSA is better for longer strings as it allows for transpositions which are common in longer descriptions
    gene_product_dist_bakta_prokka = stringdist(tolower(gene_product_bakta), tolower(gene_product_prokka), method = "osa"),
    gene_product_dist_bakta_micro = stringdist(tolower(gene_product_bakta), tolower(gene_product_microbeannotator), method = "osa"),
    gene_product_dist_prokka_micro = stringdist(tolower(gene_product_prokka), tolower(gene_product_microbeannotator), method = "osa"),
    
    # Determine consensus gene product
    # Algorithm:
    # 1. Prefer non-hypothetical proteins
    # 2. If all are hypothetical or NA, prefer in order: bakta, prokka, microbeannotator
    gene_product_SC = case_when(
      !is.na(gene_product_bakta) & gene_product_bakta != "hypothetical protein" ~ gene_product_bakta,
      !is.na(gene_product_prokka) & gene_product_prokka != "hypothetical protein" ~ gene_product_prokka,
      !is.na(gene_product_microbeannotator) & gene_product_microbeannotator != "hypothetical protein" ~ gene_product_microbeannotator,
      !is.na(gene_product_bakta) ~ gene_product_bakta,
      !is.na(gene_product_prokka) ~ gene_product_prokka,
      !is.na(gene_product_microbeannotator) ~ gene_product_microbeannotator,
      TRUE ~ NA_character_
    ),
    
    # Determine confidence level of consensus
    # Algorithm:
    # 1. If all three sources agree, mark as "complete_consensus"
    # 2. If Bakta is very similar (distance <= 10) to at least one other source, mark as "high_probability_consensus"
    # 3. If Bakta is somewhat similar (10 < distance <= 20) to at least one other source, mark as "medium_probability_consensus"
    # 4. If Bakta is less similar to other sources (20 < distance <= 25), mark as "low_probability_consensus"
    # 5. If Bakta is very different from other sources (distance > 25), mark as "no_consensus"
    # 6. If only Bakta is available, mark as "only_bakta_result"
    # 7. If only one source (not Bakta) is available, mark as "only_one_result"
    
    gene_product_SC_confidence = case_when(
      # Complete consensus among all three
      !is.na(gene_product_bakta) & !is.na(gene_product_prokka) & !is.na(gene_product_microbeannotator) &
        (gene_product_bakta == gene_product_prokka & gene_product_bakta == gene_product_microbeannotator) ~ "complete_consensus",
      
      # High probability consensus
      # Bakta is very similar (distance <= 10) to at least one other source
      !is.na(gene_product_bakta) & 
        (gene_product_dist_bakta_prokka <= 10 | gene_product_dist_bakta_micro <= 10) ~ "high_probability_consensus",
      
      # Medium probability consensus
      # Bakta is somewhat similar (15 < distance <= 30) to at least one other source
      !is.na(gene_product_bakta) & 
        ((gene_product_dist_bakta_prokka > 10 & gene_product_dist_bakta_prokka <= 20) |
           (gene_product_dist_bakta_micro > 10 & gene_product_dist_bakta_micro <= 20)) ~ "medium_probability_consensus",
      
      # Low probability consensus
      # Bakta is present but different from other sources (30 < distance <= 45)
      !is.na(gene_product_bakta) & 
        ((gene_product_dist_bakta_prokka > 20 & gene_product_dist_bakta_prokka <= 25) |
           (gene_product_dist_bakta_micro > 20 & gene_product_dist_bakta_micro <= 25)) ~ "low_probability_consensus",
      
      # No consensus
      !is.na(gene_product_bakta) & 
        (gene_product_dist_bakta_prokka > 25 | gene_product_dist_bakta_micro > 25) ~ "no_consensus",
      
      # Only Bakta result
      (!is.na(gene_product_bakta) & is.na(gene_product_prokka) & is.na(gene_product_microbeannotator)) ~ "only_bakta_result",
      
      # Only one result (not Bakta)
      (is.na(gene_product_bakta) & (!is.na(gene_product_prokka) | !is.na(gene_product_microbeannotator))) ~ "only_one_result",
      
      TRUE ~ NA_character_
    )
  )

# Create overall confidence variable
gene_data_raw <- gene_data_raw %>%
  mutate(
    # Convert confidence levels to numeric scores
    gene_confidence_score = case_when(
      gene_SC_confidence == "complete_consensus" ~ 5,
      gene_SC_confidence == "high_probability_consensus" ~ 4,
      gene_SC_confidence == "medium_probability_consensus" ~ 3,
      gene_SC_confidence == "low_probability_consensus" ~ 2,
      gene_SC_confidence == "no_consensus" ~ 1,
      gene_SC_confidence == "only_one_result" ~ 3,  # Moderate confidence for single source
      TRUE ~ 0  # For NA or unexpected values
    ),
    product_confidence_score = case_when(
      gene_product_SC_confidence == "complete_consensus" ~ 5,
      gene_product_SC_confidence == "high_probability_consensus" ~ 4,
      gene_product_SC_confidence == "medium_probability_consensus" ~ 3,
      gene_product_SC_confidence == "low_probability_consensus" ~ 2,
      gene_product_SC_confidence == "no_consensus" ~ 1,
      gene_product_SC_confidence == "only_one_result" ~ 3,  # Moderate confidence for single source
      TRUE ~ 0  # For NA or unexpected values
    ),
    
    # Calculate weighted average, giving  more weight to gene name confidence
    overall_confidence_score = (gene_confidence_score * 0.8 + product_confidence_score * 0.2),
    
    # Convert score back to categorical variable
    overall_gene_SC_confidence = case_when(
      overall_confidence_score >= 4.5 ~ "very_high",
      overall_confidence_score >= 3.5 ~ "high",
      overall_confidence_score >= 2.5 ~ "medium",
      overall_confidence_score >= 1.5 ~ "low",
      overall_confidence_score > 0 ~ "very_low",
      TRUE ~ "unknown"
    )
  )

# Complete consensus between tools for either gene name or gene product should give very high confidence 
gene_data_raw$overall_gene_SC_confidence[gene_data_raw$gene_SC_confidence == "complete_consensus"] <- "very_high"
gene_data_raw$overall_gene_SC_confidence[gene_data_raw$gene_product_SC_confidence == "complete_consensus"] <- "very_high"

# Update annotation_results_integrated with new columns
if (nrow(gene_data_raw) > 0) {
  annotation_results_integrated <- annotation_results_integrated %>%
    select(-all_of(setdiff(available_columns, "locus_tag_SC"))) %>%
    left_join(gene_data_raw %>% select(locus_tag_SC, ends_with(c("_SC", "all_gene_SC_confidence"))),
              by = "locus_tag_SC")
  
  # Attach detailed comments to the final columns in annotation_results_integrated
  comment(annotation_results_integrated$gene_SC) <- "gene_SC: Consensus gene name derived from gene_bakta and gene_prokka. 
  Algorithm: 
  1. If only one gene name is available, use that. 
  2. If both are available and identical, use that. 
  3. If both are available but different, prefer bakta's gene name. 
  4. If neither is available, use NA."
  
  comment(gene_data_raw$gene_SC_confidence) <- "gene_SC_confidence: Confidence level of the gene name consensus, based on Jaro-Winkler distance between bakta and prokka gene names. 
  Levels: 
  1. 'only_one_result': Only one source provided a gene name. 
  2. 'complete_consensus': Both sources agree completely (distance = 0). 
  3. 'high_probability_consensus': Very similar names (0 < distance <= 0.10). 
  4. 'medium_probability_consensus': Somewhat similar names (0.10 < distance <= 0.20). 
  5. 'low_probability_consensus': Less similar names (0.20 < distance <= 0.25). 
  6. 'no_consensus': Very different names (distance > 0.25)."
  
  comment(annotation_results_integrated$gene_product_SC) <- "gene_product_SC: Consensus gene product derived from gene_product_bakta, gene_product_prokka, and gene_product_microbeannotator. 
  Algorithm: 
  1. Prefer non-hypothetical proteins in the order: bakta, prokka, microbeannotator. 
  2. If all are hypothetical or NA, prefer in order: bakta, prokka, microbeannotator. 
  3. If all sources are NA, use NA."
  
  comment(gene_data_raw$gene_product_SC_confidence) <- "gene_product_SC_confidence: Confidence level of the gene product consensus, based on pairwise Optimal String Alignment (OSA) distances between gene products from all three sources. 
  Levels: 
  1. 'complete_consensus': All three sources agree completely. 
  2. 'high_probability_consensus': Two sources agree exactly, or all sources are very similar (all pairwise distances <= 10). 
  3. 'medium_probability_consensus': Sources are somewhat similar (10 < at least one distance < 20). 
  4. 'low_probability_consensus': Sources are less similar (20 <= at least one distance < 25). 
  5. 'no_consensus': Sources are very different (at least one distance >= 25). 
  6. 'only_one_result': Only one source provided a gene product."
  
  comment(annotation_results_integrated$overall_gene_SC_confidence) <- "overall_confidence: Composite measure of annotation confidence based on both gene name and gene product confidences. 
  Algorithm:
  1. Convert gene_SC_confidence and gene_product_SC_confidence to numeric scores (1-5).
  2. Calculate weighted average: 40% gene confidence, 60% product confidence.
  3. Convert weighted average to categorical confidence level.
  Levels: 'very_high' (≥4.5), 'high' (≥3.5), 'medium' (≥2.5), 'low' (≥1.5), 'very_low' (>0), 'unknown' (0)."
  
} else {
  message("No gene_data_raw generated. Skipping the join and confidence creation steps.")
}

gene_data_raw_scores <- gene_data_raw %>% select(all_of(available_columns), ends_with(c("_SC", "SC_confidence")))

# Clean up intermediate variables
rm(available_columns, gene_data_raw)

# Final Result: annotation_results_integrated now contains the new gene and gene product columns and confidence levels


# Integrage gene types of Bakta and Prokka
if ("type_bakta" %in% colnames(annotation_results_integrated)) {
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "cds"] <- "CDS"
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "sorf"] <- "sORF"
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "crispr"] <- "CRISPR"
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "crispr-repeat"] <- "CRISPR_repeat"
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "crispr-spacer"] <- "CRISPR_spacer"
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "ncRNA"] <- "misc_RNA"
  annotation_results_integrated$type_bakta[annotation_results_integrated$type_bakta == "ncRNA-region"] <- "misc_RNA_region"
}

annotation_results_integrated <- annotation_results_integrated %>%
  mutate(
    type_prokka = if ("type_prokka" %in% colnames(.)) as.character(type_prokka) else NA_character_,
    type_bakta = if ("type_bakta" %in% colnames(.)) as.character(type_bakta) else NA_character_,
    type_SC = case_when(
      !is.na(type_bakta) ~ type_bakta,
      !is.na(type_prokka) ~ type_prokka,
      TRUE ~ NA_character_
    )
  )

# Last details
# Add prefix "COG" to each value in COG_number_SC, keeping NA as NA
annotation_results_integrated$COG_number_SC <- ifelse(
  is.na(annotation_results_integrated$COG_number_SC),
  NA,
  paste0("COG", annotation_results_integrated$COG_number_SC)
)

# Replace patterns COGCOG:COG, COG:COG, and COG: with COG
annotation_results_integrated$COG_number_SC <- gsub(
  pattern = "COGCOG:COG|COG:COG|COG:",
  replacement = "COG",
  x = annotation_results_integrated$COG_number_SC
)

# Add prefix "GO:" to each number in GO_number_bakta, keeping NA as NA
annotation_results_integrated$GO_number_bakta <- sapply(annotation_results_integrated$GO_number_bakta, function(x) {
  if (is.na(x)) {
    return(NA)
  } else {
    return(paste0("GO:", unlist(strsplit(x, "; ")), collapse = "; "))
  }
})

# Add prefix "SO:" to each number in SO_number_bakta, keeping NA as NA
annotation_results_integrated$SO_number_bakta <- sapply(annotation_results_integrated$SO_number_bakta, function(x) {
  if (is.na(x)) {
    return(NA)
  } else {
    return(paste0("SO:", unlist(strsplit(x, "; ")), collapse = "; "))
  }
})

# Define the desired order of columns
desired_column_order <- c(
  "contig_tag_fasta",
  "contig_tag_prokka",
  "contig_tag_bakta",
  "contig_tag_SC",
  "locus_tag_prokka", 
  "locus_tag_bakta", 
  "locus_tag_SC", 
  "gene_SC", 
  "gene_product_SC", 
  "overall_gene_SC_confidence",
  "start_position", 
  "end_position", 
  "length_bp", 
  "strand",
  "EC_number_SC_best", 
  "EC_number_SC_all", 
  "COG_category_SC", 
  "COG_number_SC", 
  "K_number_SC",
  "type_prokka",
  "type_bakta",
  "type_SC",
  "GO_number_bakta", 
  "SO_number_bakta",
  "nucleotide_code", 
  "amino_acid_code",
  "UniRef100_bakta", 
  "UniRef90_bakta", 
  "UniRef50_bakta", 
  "UniParc_ID_bakta", 
  "RefSeq_bakta", 
  "NCBIProtein_bakta", 
  "VFDB_ID_bakta", 
  "PFAM_bakta", 
  "RFAM_bakta", 
  "NCBIFam_bakta",
  "IS_bakta",
  "BlastRules_bakta", 
  "database_microbeannotator", 
  "source_file_prokka",
  "source_file_bakta",
  "source_file_microbeannotator",
  "software_version"
)

# Get the current column names of the dataframe
current_columns <- colnames(annotation_results_integrated)

# Find which of the desired columns are actually present in the dataframe
available_columns <- intersect(desired_column_order, current_columns)

# Find any columns in the dataframe that weren't in our desired order
remaining_columns <- setdiff(current_columns, desired_column_order)

# Combine the available columns (in our desired order) with any remaining columns
new_column_order <- c(available_columns, remaining_columns)

# Reorder the columns in the dataframe
annotation_results_integrated <- annotation_results_integrated[, new_column_order]

# Print information about the reordering process
cat("Columns reordered successfully.\n")
cat("Number of columns in desired order:", length(available_columns), "\n")
cat("Number of additional columns:", length(remaining_columns), "\n")

if (length(remaining_columns) > 0) {
  cat("Additional columns:\n")
  print(remaining_columns)
}

rm(new_column_order)

# At the end of your script, save the results:
output_file <- file.path(opt$output_dir, "annotation_results_integrated.qs")
qsave(annotation_results_integrated, output_file, preset = "archive")
cat(sprintf("SUCCESS: Results saved to %s\n", output_file))