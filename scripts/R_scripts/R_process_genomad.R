#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# R_process_genomad.R

set.seed(42)

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(qs)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="Output directory for qs files", metavar="CHARACTER"),
  make_option(c("-s", "--virus_summary"), type="character", default=NULL, 
              help="Input virus_summary.tsv file path", metavar="CHARACTER"),
  make_option(c("-g", "--virus_genes"), type="character", default=NULL, 
              help="Input virus_genes.tsv file path", metavar="CHARACTER"),
  make_option(c("-p", "--plasmid_summary"), type="character", default=NULL, 
              help="Input plasmid_summary.tsv file path (optional)", metavar="CHARACTER"),
  make_option(c("-l", "--plasmid_genes"), type="character", default=NULL, 
              help="Input plasmid_genes.tsv file path (optional)", metavar="CHARACTER"),
  make_option(c("-v", "--version"), type="character", default=NULL, 
              help="StrainCascade version", metavar="CHARACTER")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$output_dir) || is.null(opt$virus_summary) || is.null(opt$virus_genes) || is.null(opt$version)) {
  stop("Required arguments (output_dir, virus_summary, virus_genes, and version) must be supplied.")
}

# Function to safely read files
safe_read <- function(file_path, read_func, ...) {
  if (is.null(file_path) || !file.exists(file_path)) {
    warning(paste("File does not exist or was not provided:", file_path))
    return(data.frame())
  }
  tryCatch({
    df <- read_func(file_path, ...)
    if (nrow(df) == 0) {
      warning(paste("File is empty or contains only headers:", file_path))
    }
    return(df)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(data.frame())
  })
}

# Function to create empty data frame with specified columns
create_empty_df <- function(columns) {
  df <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(df) <- columns
  return(df)
}

# Read and process virus summary file
virus_summary <- safe_read(opt$virus_summary, read.delim, header = TRUE, sep = "\t", fill = TRUE)
if (nrow(virus_summary) > 0) {
  virus_summary <- virus_summary %>%
    dplyr::rename(
      contig_tag_genomad = seq_name,
      length_bp = length,
      topology_genomad = topology,
      coordinates_genomad = coordinates,
      n_genes_genomad = n_genes,
      genetic_code_genomad = genetic_code,
      virus_score_genomad = virus_score,
      fdr_genomad = fdr,
      n_hallmarks_genomad = n_hallmarks,
      marker_enrichment_genomad = marker_enrichment,
      taxonomy_genomad = taxonomy
    ) %>%
    mutate(
      source_file_genomad_virus_summary = basename(opt$virus_summary),
      software_version = paste0("StrainCascade v", opt$version)
    )
} else {
  virus_summary <- create_empty_df(c("contig_tag_genomad", "length_bp", "topology_genomad", 
                                      "coordinates_genomad", "n_genes_genomad", "genetic_code_genomad", 
                                      "virus_score_genomad", "fdr_genomad", "n_hallmarks_genomad", 
                                      "marker_enrichment_genomad", "taxonomy_genomad", 
                                      "source_file_genomad_virus_summary", "software_version"))
}

# Read and process virus genes file
virus_genes <- safe_read(opt$virus_genes, read.delim, header = TRUE, sep = "\t", fill = TRUE)
if (nrow(virus_genes) > 0) {
  virus_genes <- virus_genes %>%
    dplyr::rename(
      gene_id_genomad = gene,
      start_position = start,
      end_position = end,
      length_bp_gene = length,
      strand_genomad = strand,
      gc_content_genomad = gc_content,
      genetic_code_genomad = genetic_code,
      rbs_motif_genomad = rbs_motif,
      marker_genomad = marker,
      evalue_genomad = evalue,
      bitscore_genomad = bitscore,
      uscg_genomad = uscg,
      plasmid_hallmark_genomad = plasmid_hallmark,
      virus_hallmark_genomad = virus_hallmark,
      taxid_genomad = taxid,
      taxname_genomad = taxname,
      annotation_conjscan_genomad = annotation_conjscan,
      annotation_amr_genomad = annotation_amr,
      annotation_accessions_genomad = annotation_accessions,
      annotation_description_genomad = annotation_description
    ) %>%
    mutate(
      source_file_genomad_virus_genes = basename(opt$virus_genes),
      software_version = paste0("StrainCascade v", opt$version)
    )
} else {
  virus_genes <- create_empty_df(c("gene_id_genomad", "start_position", "end_position", 
                                    "length_bp_gene", "strand_genomad", "gc_content_genomad", 
                                    "genetic_code_genomad", "rbs_motif_genomad", "marker_genomad", 
                                    "evalue_genomad", "bitscore_genomad", "uscg_genomad", 
                                    "plasmid_hallmark_genomad", "virus_hallmark_genomad", 
                                    "taxid_genomad", "taxname_genomad", "annotation_conjscan_genomad", 
                                    "annotation_amr_genomad", "annotation_accessions_genomad", 
                                    "annotation_description_genomad", "source_file_genomad_virus_genes", 
                                    "software_version"))
}

# Extract contig_tag_fasta from gene_id_genomad
if (nrow(virus_genes) > 0) {
  virus_genes$contig_tag_genomad <- sub("_\\d+$", "", virus_genes$gene_id_genomad)
}

# Merge virus data
genomad_virus_results <- create_empty_df(c("contig_tag_fasta", "contig_tag_genomad", "topology_genomad", 
                                           "coordinates_genomad", "n_genes_genomad", "length_bp", 
                                           "genetic_code_genomad", "virus_score_genomad", "fdr_genomad", 
                                           "n_hallmarks_genomad", "marker_enrichment_genomad", 
                                           "taxonomy_genomad", "gene_id_genomad", "start_position", 
                                           "end_position", "length_bp_gene", "strand_genomad", 
                                           "gc_content_genomad", "rbs_motif_genomad", "marker_genomad", 
                                           "evalue_genomad", "bitscore_genomad", "uscg_genomad", 
                                           "plasmid_hallmark_genomad", "virus_hallmark_genomad", 
                                           "taxid_genomad", "taxname_genomad", "annotation_conjscan_genomad", 
                                           "annotation_amr_genomad", "annotation_accessions_genomad", 
                                           "annotation_description_genomad", "software_version", 
                                           "source_file_genomad_virus_summary", "source_file_genomad_virus_genes"))

if (nrow(virus_summary) > 0 || nrow(virus_genes) > 0) {
  genomad_virus_results <- full_join(virus_summary, virus_genes, 
                                     by = c("contig_tag_genomad", "genetic_code_genomad", "software_version"))
  
  # Extract contig_tag_fasta (provirus format: contig|provirus_start_end)
  genomad_virus_results$contig_tag_fasta <- sub("\\|provirus.*", "", genomad_virus_results$contig_tag_genomad)
  
  # Ensure all desired columns exist
  for (col in colnames(genomad_virus_results)) {
    if (!(col %in% names(genomad_virus_results))) {
      genomad_virus_results[[col]] <- NA
    }
  }
}

# Read and process plasmid summary file (optional)
plasmid_summary <- safe_read(opt$plasmid_summary, read.delim, header = TRUE, sep = "\t", fill = TRUE)
if (nrow(plasmid_summary) > 0) {
  plasmid_summary <- plasmid_summary %>%
    dplyr::rename(
      contig_tag_genomad = seq_name,
      length_bp = length,
      topology_genomad = topology,
      n_genes_genomad = n_genes,
      genetic_code_genomad = genetic_code,
      plasmid_score_genomad = plasmid_score,
      fdr_genomad = fdr,
      n_hallmarks_genomad = n_hallmarks,
      marker_enrichment_genomad = marker_enrichment,
      conjugation_genes_genomad = conjugation_genes,
      amr_genes_genomad = amr_genes
    ) %>%
    mutate(
      source_file_genomad_plasmid_summary = basename(opt$plasmid_summary),
      software_version = paste0("StrainCascade v", opt$version)
    )
} else {
  plasmid_summary <- create_empty_df(c("contig_tag_genomad", "length_bp", "topology_genomad", 
                                        "n_genes_genomad", "genetic_code_genomad", "plasmid_score_genomad", 
                                        "fdr_genomad", "n_hallmarks_genomad", "marker_enrichment_genomad", 
                                        "conjugation_genes_genomad", "amr_genes_genomad", 
                                        "source_file_genomad_plasmid_summary", "software_version"))
}

# Read and process plasmid genes file (optional)
plasmid_genes <- safe_read(opt$plasmid_genes, read.delim, header = TRUE, sep = "\t", fill = TRUE)
if (nrow(plasmid_genes) > 0) {
  plasmid_genes <- plasmid_genes %>%
    dplyr::rename(
      gene_id_genomad = gene,
      start_position = start,
      end_position = end,
      length_bp_gene = length,
      strand_genomad = strand,
      gc_content_genomad = gc_content,
      genetic_code_genomad = genetic_code,
      rbs_motif_genomad = rbs_motif,
      marker_genomad = marker,
      evalue_genomad = evalue,
      bitscore_genomad = bitscore,
      uscg_genomad = uscg,
      plasmid_hallmark_genomad = plasmid_hallmark,
      virus_hallmark_genomad = virus_hallmark,
      taxid_genomad = taxid,
      taxname_genomad = taxname,
      annotation_conjscan_genomad = annotation_conjscan,
      annotation_amr_genomad = annotation_amr,
      annotation_accessions_genomad = annotation_accessions,
      annotation_description_genomad = annotation_description
    ) %>%
    mutate(
      source_file_genomad_plasmid_genes = basename(opt$plasmid_genes),
      software_version = paste0("StrainCascade v", opt$version)
    )
  
  # Extract contig_tag_genomad from gene_id_genomad
  plasmid_genes$contig_tag_genomad <- sub("_\\d+$", "", plasmid_genes$gene_id_genomad)
} else {
  plasmid_genes <- create_empty_df(c("gene_id_genomad", "start_position", "end_position", 
                                      "length_bp_gene", "strand_genomad", "gc_content_genomad", 
                                      "genetic_code_genomad", "rbs_motif_genomad", "marker_genomad", 
                                      "evalue_genomad", "bitscore_genomad", "uscg_genomad", 
                                      "plasmid_hallmark_genomad", "virus_hallmark_genomad", 
                                      "taxid_genomad", "taxname_genomad", "annotation_conjscan_genomad", 
                                      "annotation_amr_genomad", "annotation_accessions_genomad", 
                                      "annotation_description_genomad", "source_file_genomad_plasmid_genes", 
                                      "software_version"))
}

# Merge plasmid data
genomad_plasmid_results <- create_empty_df(c("contig_tag_fasta", "contig_tag_genomad", "topology_genomad", 
                                             "n_genes_genomad", "length_bp", "genetic_code_genomad", 
                                             "plasmid_score_genomad", "fdr_genomad", "n_hallmarks_genomad", 
                                             "marker_enrichment_genomad", "conjugation_genes_genomad", 
                                             "amr_genes_genomad", "gene_id_genomad", "start_position", 
                                             "end_position", "length_bp_gene", "strand_genomad", 
                                             "gc_content_genomad", "rbs_motif_genomad", "marker_genomad", 
                                             "evalue_genomad", "bitscore_genomad", "uscg_genomad", 
                                             "plasmid_hallmark_genomad", "virus_hallmark_genomad", 
                                             "taxid_genomad", "taxname_genomad", "annotation_conjscan_genomad", 
                                             "annotation_amr_genomad", "annotation_accessions_genomad", 
                                             "annotation_description_genomad", "software_version", 
                                             "source_file_genomad_plasmid_summary", "source_file_genomad_plasmid_genes"))

if (nrow(plasmid_summary) > 0 || nrow(plasmid_genes) > 0) {
  genomad_plasmid_results <- full_join(plasmid_summary, plasmid_genes, 
                                       by = c("contig_tag_genomad", "genetic_code_genomad", "software_version"))
  
  # Extract contig_tag_fasta
  genomad_plasmid_results$contig_tag_fasta <- genomad_plasmid_results$contig_tag_genomad
  
  # Ensure all desired columns exist
  for (col in colnames(genomad_plasmid_results)) {
    if (!(col %in% names(genomad_plasmid_results))) {
      genomad_plasmid_results[[col]] <- NA
    }
  }
}

# Save virus results as qs
virus_output_file <- file.path(opt$output_dir, "genomad_virus_results.qs")
qsave(genomad_virus_results, virus_output_file, preset = "archive")
cat("geNomad virus results processed and saved to:", virus_output_file, "\n")
cat("Number of rows in genomad_virus_results:", nrow(genomad_virus_results), "\n")

# Save plasmid results as qs
plasmid_output_file <- file.path(opt$output_dir, "genomad_plasmid_results.qs")
qsave(genomad_plasmid_results, plasmid_output_file, preset = "archive")
cat("geNomad plasmid results processed and saved to:", plasmid_output_file, "\n")
cat("Number of rows in genomad_plasmid_results:", nrow(genomad_plasmid_results), "\n")
