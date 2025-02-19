#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_gtdbtk_taxonomy.R

set.seed(42)

# Load necessary libraries
suppressPackageStartupMessages({
  library(ggtree)
  library(ape)
  library(dplyr)
  library(tidytree)
  library(optparse)
  library(ggplot2)
  library(qs)
})


# Parse command-line arguments
option_list <- list(
  make_option("--output_dir", type="character", default=NULL, help="Output directory for data (QS file)"),
  make_option("--output_dir2", type="character", default=NULL, help="Output directory for plots (PDF file)"),
  make_option("--tree", type="character", default=NULL, help="Path to GTDB-Tk tree file"),
  make_option("--tsv", type="character", default=NULL, help="Path to GTDB-Tk summary TSV file"),
  make_option("--sample_name", type="character", default=NULL, help="Name of the input sample"),
  make_option("--version", type="character", default=NULL, help="StrainCascade version")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to replace "NA" strings with actual NA values, excluding date-time columns
replace_na_strings <- function(df) {
  date_cols <- sapply(df, inherits, what = c("POSIXt", "Date"))
  na_pattern <- "^(NA|N/A|N/a|n/A|n/a|)$"
  df[!date_cols] <- df[!date_cols] %>%
    mutate(across(everything(), ~ replace(.x, grepl(na_pattern, .x), NA)))
  return(df)
}

# Function to process TSV file
process_tsv_file <- function(tsv_file, sample_name, version) {
  summary_tsv <- read.delim(tsv_file, header = TRUE, sep = "\t", fill = TRUE)
  summary_tsv <- replace_na_strings(summary_tsv)
  
  colnames(summary_tsv) <- paste0(colnames(summary_tsv), "_gtdbtk_classification")
  
  summary_tsv$sample_name_gtdbtk_classification <- sample_name
  summary_tsv$source_file_gtdbtk_classification <- basename(tsv_file)
  summary_tsv$software_version <- paste0("StrainCascade v", version)
  
  return(summary_tsv)
}

# Function to create phylogenetic tree plot
create_tree_plot <- function(tree_file, sample_name, taxonomic_classification) {
  tree <- read.tree(tree_file)
  tree <- as.phylo(tree)
  
  filtered_tips <- tree$tip.label[grepl(sample_name, tree$tip.label)]
  tip_of_interest <- filtered_tips[1]
  
  n_neighbors <- 25
  tip_idx <- match(tip_of_interest, tree$tip.label)
  
  dist_matrix <- cophenetic(tree)
  neighbor_idxs <- order(dist_matrix[tip_idx, ])[1:(n_neighbors + 1)]
  neighbor_tips <- tree$tip.label[neighbor_idxs]
  
  subtree <- keep.tip(tree, neighbor_tips, trim.internal = TRUE)
  nudge <- max(subtree$edge.length) * 0.06
  subtree$tip.label[subtree$tip.label == filtered_tips] <- sample_name
  
  p <- ggtree(subtree) +
    geom_tiplab(aes(label = label), align = TRUE, linetype = "dotted") +
    geom_nodelab(size = 3, nudge_x = -1*nudge, nudge_y = 0.5) +
    theme_tree2() +
    ggtitle("Phylogenetic tree of closest taxa", 
            subtitle = paste("Analysed genome: ", sample_name, 
                             "\nResult: ", taxonomic_classification,
                             "\nTool: GTDB-Tk classification (classify_wf)", 
                             "\nInternal node values: non-parametric bootstrap support values", sep = "")) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 10, face = "italic"),
      axis.title.x = element_text(size = 10, face = "bold")
    ) +
    xlab("Expected substitutions per site")
  
  max_x <- max(p$data$x)
  p <- p + xlim(0, (max_x*1.4))
  
  p <- p + geom_tippoint(aes(color = ifelse(label == sample_name, "#A94322", "black"),
                             size = ifelse(label == sample_name, 3, 1)))
  
  p <- p + scale_color_manual(values = c("black" = "black", "#A94322" = "#A94322"),
                              labels = c("#A94322" = "Analysed genome", "black" = "25 closest reference genomes"),
                              name = "Taxa") +
    scale_size_continuous(range = c(1, 3), guide = "none")
  
  return(list(plot = p, tree_file = tree))
}

# Main function
process_gtdbtk_data <- function(tree_file, tsv_file, output_dir, output_dir2, sample_name, version) {
  processed_tsv <- process_tsv_file(tsv_file, sample_name, version)
  taxonomic_classification <- processed_tsv$classification_gtdbtk_classification[1]
  
  tree_results <- create_tree_plot(tree_file, sample_name, taxonomic_classification)
  
  tree_plot <- tree_results$plot
  newick_tree <- tree_results$tree_file

  return(list(data = processed_tsv, newick_tree = newick_tree, plot = tree_plot))
}

# Main execution
tryCatch({
  if (is.null(opt$output_dir) || is.null(opt$output_dir2) || is.null(opt$tree) || is.null(opt$tsv) || is.null(opt$sample_name) || is.null(opt$version)) {
    stop("All arguments (output_dir, output_dir2, tree, tsv, sample_name, and version) must be provided.")
  }
  
  result <- process_gtdbtk_data(opt$tree, opt$tsv, opt$output_dir, opt$output_dir2, opt$sample_name, opt$version)
  
  # Save processed TSV as QS file
  qs_file <- file.path(opt$output_dir, "gtdbtk_classification_results.qs")
  qsave(result$data, qs_file, preset = "archive")
  cat(sprintf("SUCCESS: QS (dataframe) file saved to %s\n", qs_file))

    qs_file2 <- file.path(opt$output_dir, "gtdbtk_classification_newick_tree.qs")
  qsave(result$newick_tree, qs_file2, preset = "archive")
  cat(sprintf("SUCCESS: QS (Newick) file saved to %s\n", qs_file2))

    qs_file3 <- file.path(opt$output_dir, "gtdbtk_25_tree_plot.qs")
  qsave(result$plot, qs_file3, preset = "archive")
  cat(sprintf("SUCCESS: QS (ggplot) file saved to %s\n", qs_file3))
  
  # Save the plot as PDF
  pdf_file <- file.path(opt$output_dir2, "gtdbtk_25_tree_plot.pdf")
  ggsave(pdf_file, result$plot, width = 12, height = 8)
  cat(sprintf("SUCCESS: PDF plot saved to %s\n", pdf_file))
  
}, error = function(e) {
  cat(sprintf("ERROR: An unexpected error occurred: %s\n", e$message))
  quit(save = "no", status = 0)
})