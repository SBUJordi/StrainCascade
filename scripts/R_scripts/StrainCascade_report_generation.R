#!/usr/bin/env Rscript

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_report_generation.R

set.seed(42)

library(optparse)

# Define command line options
option_list <- list(
  make_option(c("--sample_name"), type="character", default=NULL, 
              help="Sample name", metavar="character"),
  make_option(c("--RData"), type="character", default=NULL,
              help="Path to RData file", metavar="character"),
  make_option(c("--output_dir"), type="character", default=NULL, 
              help="Output directory path", metavar="character")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$sample_name) || is.null(opt$RData) || is.null(opt$output_dir)) {
  stop("All arguments (sample_name, RData, output_dir) must be supplied.")
}

# Assign parsed arguments to variables
sample_name <- opt$sample_name
rdata_file <- opt$RData
output_dir <- opt$output_dir

generate_rmarkdown <- function() {
  rmd_content <- "
---
title: 'StrainCascade'
output: 
  flexdashboard::flex_dashboard:
    orientation: columns
    vertical_layout: fill
    theme: readable
params:
  rdata_file: !r rdata_file
---
  
```{css, echo=FALSE}
:root {
  --main-navbar-color: #336392; /* Dark gray with 90% opacity */
}

/* Set border-top-color to cutsom color for active tabs */
.nav-tabs-custom > .nav-tabs > li.active {
  border-top-color: #336392 !important;
}

/* Keep the default font for navigation elements */
.navbar, .nav {
  font-family: 'Source Sans Pro', sans-serif !important;
}

/* Apply custom font to table contents and titles */
.dataTables_wrapper, .dataTable, .table, h3 {
  font-family: 'Source Sans Pro', sans-serif !important;
}

/* Apply custom font to all headings */
h1, h2, h3, h4, h5, h6 {
  font-family: 'Source Sans Pro', sans-serif !important;
}

/* Specifically target the titles we want to change */
.chart-title, .chart-stage h3 {
  font-family: 'Source Sans Pro', sans-serif !important;
}

/* Customize main navbar */
.navbar-inverse {
  background-color: var(--main-navbar-color) !important;
  border-color: var(--main-navbar-color) !important;
}

/* Ensure text is readable on custom backgrounds */
.navbar-inverse .navbar-brand,
.navbar-inverse .navbar-nav > li > a {
  color: white !important;
```

```{r setup, include=FALSE}
suppressPackageStartupMessages({
  library(flexdashboard)
  library(DT)
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(kableExtra)
  library(knitr)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

knitr::opts_chunk$set(echo = FALSE)
load(params$rdata_file)

# Function to safely get data
safe_get <- function(data_name, default = NULL) {
  if (exists(data_name)) {
    get(data_name)
  } else {
    default
  }
}

# Function to set custom font for plots
set_font <- function(p) {
  p %>% layout(font = list(family = 'Source Sans Pro, sans-serif'))
}

```

Overview
=====================================

Column {data-width=600}
-------------------------------------

### Overview table

```{r overview_table}
# Initialize a list to hold the available results
results <- list()

# Retrieve data from each source and add available data to the results list
gtdbtk_data <- safe_get('gtdbtk_classification_results')
if (!is.null(gtdbtk_data)) {
  results[['Analysed genome']] <- gtdbtk_data$user_genome_gtdbtk_classification %||% 'N/A'
}

quast_data <- safe_get('quast_results')
if (!is.null(quast_data)) {
  results[['Number of contigs']] <- quast_data$number_of_contigs_quast %||% 'N/A'
  results[['Genome size (bp)']] <- quast_data$total_length_quast %||% 'N/A'
}

checkm2_data <- safe_get('checkm2_results')
if (!is.null(checkm2_data)) {
  results[['Completeness (%)']] <- round(checkm2_data$completeness_checkm2, 2) %||% 'N/A'
  results[['Contamination (%)']] <- round(checkm2_data$contamination_checkm2, 2) %||% 'N/A'
}

gtdbtk_data <- safe_get('gtdbtk_classification_results')
if (!is.null(gtdbtk_data)) {
  results[['Taxonomy']] <- gtdbtk_data$classification_gtdbtk_classification %||% 'N/A'
}

amrfinderplus_data <- safe_get('amrfinderplus_results')
if (!is.null(amrfinderplus_data)) {
  results[['Number of AMR genes (AMRFinderPlus)']] <- nrow(amrfinderplus_data)
}

resfinder_data <- safe_get('resfinder_results')
if (!is.null(resfinder_data)) {
  results[['Number of AMR genes (ResFinder)']] <- nrow(resfinder_data)
}

islandpath_data <- safe_get('islandpath_results')
if (!is.null(islandpath_data)) {
  results[['Number of Genomic Islands']] <- nrow(islandpath_data)
}

quast_data <- safe_get('quast_results')
if (!is.null(quast_data)) {
  results[['Software']] <- quast_data$software_version %||% 'N/A'
}

# Add the report date
results[['Report date']] <- as.character(Sys.Date())

# Check if we have any results to display
if (length(results) > 1) {
  # Convert the results to a data frame for display
  summary_table <- data.frame(Value = unlist(results))
  
  # Print the summary table using kable
  kable(summary_table, class = 'table') %>%
    kable_styling(bootstrap_options = c('hover', 'condensed', 'responsive'))
} else {
  # If no data is available, print the warning
  cat('No data available for a summary table.')
}
```

Assembly
=====================================

Column {data-width=300}
-------------------------------------

### Assembly information by QUAST

```{r assembly_quast}
quast_data <- safe_get('quast_results')
if (!is.null(quast_data)) {
  
  # Define the nicer formatted column names
  formatted_col_names <- c(
    'Number of contigs',
    'Largest contig',
    'Total length',
    'GC (%)',
    'N50',
    'N90',
    'Number of Ns per 100kbp'
    )

  selected_quast_data <- quast_data %>%
    mutate(
      number_of_contigs_quast = round(number_of_contigs_quast, 0),
      largest_contig_quast = round(largest_contig_quast, 0),
      total_length_quast = round(total_length_quast, 0),
      GC_percentage_quast = round(GC_percentage_quast, 2),
      N50_quast = round(N50_quast, 0),
      N90_quast = round(N90_quast, 0),
      number_of_N_per100kbp_quast = round(number_of_N_per100kbp_quast, 0)
    )

  selected_quast_data <- t(selected_quast_data[, c('number_of_contigs_quast', 'largest_contig_quast', 'total_length_quast', 'GC_percentage_quast', 'N50_quast', 'N90_quast', 'number_of_N_per100kbp_quast')])  
  rownames(selected_quast_data) <- formatted_col_names

  kable(selected_quast_data,
        col.names = c('Value'))
} else {
  cat('QUAST assembly information not available.')
}
```

### Quality assessment by CheckM2

```{r checkm2}
checkm2_data <- safe_get('checkm2_results')
if (!is.null(checkm2_data)) {
  # Define the nicer formatted column names
  formatted_col_names <- c(
    'Genome completeness (%)',
    'Genome contamination (%)',
    'Coding density (%)'
    )
  
  selected_checkm2_data <- checkm2_data %>%
    mutate(
      completeness_checkm2 = round(completeness_checkm2, 2),
      contamination_checkm2 = round(contamination_checkm2, 2),
      coding_density_checkm2 = round((coding_density_checkm2*100), 2)
    )

  selected_checkm2_data <- t(checkm2_data[, c('completeness_checkm2', 'contamination_checkm2', 'coding_density_checkm2')])
  rownames(selected_checkm2_data) <- formatted_col_names

  kable(selected_checkm2_data,
        col.names = c('Value'))
} else {
  cat('CheckM2 results not available.')
}
```

Column {data-width=500}
-------------------------------------

### Assembly information by StrainCascade and Circlator

```{r assembly_sc_circlator}
assembly_data <- safe_get('selected_assembly_results')
if (!is.null(assembly_data)) {
  
  SC_circlator_assembly_data <- assembly_data |>
    select(
      contig_tag_fasta,
      contig_length_bp,
      A_percentage,
      C_percentage,
      G_percentage,
      T_percentage,
      N_percentage,
      source_file_fasta
      ) |>
    mutate(
      A_percentage = round(A_percentage, 2),
      C_percentage = round(C_percentage, 2),
      G_percentage = round(G_percentage, 2),
      T_percentage = round(T_percentage, 2),
      N_percentage = round(N_percentage, 2)
    )

  circlator_data <- safe_get('circlator_results')
  if (!is.null(circlator_data)) {  
    selected_circlator_data <- circlator_data |>
    select(
      contig_tag_fasta,
      circularisation_status,
      source_file_circlator
      ) |>
    mutate(
    circularisation_status = ifelse(circularisation_status == 1, 'Yes', ifelse(
      circularisation_status == 0, 'No', 'N/A'))
    )
  }

  if (!is.null(selected_circlator_data)) { 
    SC_circlator_assembly_data <- merge(SC_circlator_assembly_data, selected_circlator_data, by = 'contig_tag_fasta', all = TRUE)
    SC_circlator_assembly_data <- SC_circlator_assembly_data |>
    select(
      contig_tag_fasta,
      circularisation_status,
      contig_length_bp,
      A_percentage,
      C_percentage,
      G_percentage,
      T_percentage,
      N_percentage,
      source_file_fasta,
      source_file_circlator
      )|> 
    dplyr::rename(
      `Contig` = contig_tag_fasta,
      `Circularised` = circularisation_status,
      `Contig length (bp)` = contig_length_bp,
      `A (%)` = A_percentage,
      `T (%)` = T_percentage,
      `G (%)` = G_percentage,
      `C (%)` = C_percentage,
      `N (%)` = N_percentage,
      `Source file (assembly)` = source_file_fasta,
      `Source file (Circlator)` = source_file_circlator,
      )
  } else {
        SC_circlator_assembly_data <- SC_circlator_assembly_data |>
    select(
      contig_tag_fasta,
      contig_length_bp,
      A_percentage,
      C_percentage,
      G_percentage,
      T_percentage,
      N_percentage,
      source_file_fasta
      )|> 
    dplyr::rename(
      `Contig` = contig_tag_fasta,
      `Contig length (bp)` = contig_length_bp,
      `A (%)` = A_percentage,
      `T (%)` = T_percentage,
      `G (%)` = G_percentage,
      `C (%)` = C_percentage,
      `N (%)` = N_percentage,
      `Source file (assembly)` = source_file_fasta
      )
  }

  kable(SC_circlator_assembly_data, align = 'c') %>%
    kable_styling(position = 'center', full_width = FALSE)

} else {
  cat('StrainCascade assembly information not available.')
}
```

### Contig length distribution by QUAST

```{r contig_length_distribution}
if (!is.null(quast_data)) {
  # Define length categories and calculate counts
  contig_data <- data.frame(
    category = factor(c('0 to 999 bp', '1000 to 4999 bp', '5000 to 9999 bp', 
                        '10000 to 24999 bp', '25000 to 49999 bp', '>=50000 bp'),
                      levels = c('0 to 999 bp', '1000 to 4999 bp', '5000 to 9999 bp', 
                                 '10000 to 24999 bp', '25000 to 49999 bp', '>=50000 bp'),
                      ordered = TRUE),
    count = c(
      quast_data$number_of_contigs_quast - quast_data$number_of_contigs_greater_or_equal_1000bp_quast,
      quast_data$number_of_contigs_greater_or_equal_1000bp_quast - quast_data$number_of_contigs_greater_or_equal_5000bp_quast,
      quast_data$number_of_contigs_greater_or_equal_5000bp_quast - quast_data$number_of_contigs_greater_or_equal_10000bp_quast,
      quast_data$number_of_contigs_greater_or_equal_10000bp_quast - quast_data$number_of_contigs_greater_or_equal_25000bp_quast,
      quast_data$number_of_contigs_greater_or_equal_25000bp_quast - quast_data$number_of_contigs_greater_or_equal_50000bp_quast,
      quast_data$number_of_contigs_greater_or_equal_50000bp_quast
    )
  )
  
  # Ensure all counts are non-negative
  contig_data$count <- pmax(0, contig_data$count)
  
  # Define a single color for all bars
  single_color <- '#336392'

  # Calculate y-axis tick values
  max_count <- max(contig_data$count)
  tick_step <- max(1, round(max_count / 10))  # Ensure step is at least 1
  y_ticks <- seq(0, max_count, by = tick_step)

  # Create the bar chart with all bars in the same color
  plot_ly(contig_data, x = ~category, y = ~count, type = 'bar', marker = list(color = single_color)) %>%
    layout(
      xaxis = list(title = 'Contig length range'),
      yaxis = list(title = 'Number of contigs', tickvals = y_ticks, ticktext = y_ticks)
    )
} else {
  cat('Contig length distribution data not available.')
}
```

Column {.tabset data-width=400}
-------------------------------------

### Identified plasmids by PlasmidFinder

```{r plasmidfinder_table}
plasmidfinder_data <- safe_get('plasmidfinder_results')

if (!is.null(plasmidfinder_data)) {
  if (all(plasmidfinder_data$hit_found_plasmidfinder == FALSE)) {
    cat('No plasmid identified.')
  } else {
    plasmidfinder_data_adapted <- plasmidfinder_data |>
      filter(hit_found_plasmidfinder != FALSE) |>
      select(
        plasmid_plasmidfinder,
        category_plasmidfinder,
        subcategory_plasmidfinder,
        contig_tag_fasta,
        start_position,
        end_position,
        identity_plasmidfinder,
        coverage_plasmidfinder,
        HSP_length_plasmidfinder,
        template_length_plasmidfinder,
        note_plasmidfinder,
        accession_plasmidfinder,
        hit_id_plasmidfinder,
        source_file_plasmidfinder
      ) |>
      rename(
        `Plasmid` = plasmid_plasmidfinder,
        `Category` = category_plasmidfinder,
        `Subcategory` = subcategory_plasmidfinder,
        `Contig tag` = contig_tag_fasta,
        `Start position` = start_position,
        `End position` = end_position,
        `Identity %` = identity_plasmidfinder,
        `Coverage %` = coverage_plasmidfinder,
        `HSP length` = HSP_length_plasmidfinder,
        `Template length` = template_length_plasmidfinder,
        `Note` = note_plasmidfinder,
        `Accession` = accession_plasmidfinder,
        `Hit ID` = hit_id_plasmidfinder,
        `Source file` = source_file_plasmidfinder
      )
    
    datatable(plasmidfinder_data_adapted,
              options = list(pageLength = 10))
  }
} else {
  cat('PlasmidFinder results not available.')
}
```

### Nucleotide code

```{r plasmidfinder_table_nc}
plasmidfinder_data <- safe_get('plasmidfinder_results')

if (!is.null(plasmidfinder_data)) {
  if (all(plasmidfinder_data$hit_found_plasmidfinder == FALSE)) {
    cat('No plasmid sequences available.')
  } else {
    plasmidfinder_sequences <- plasmidfinder_data |>
      filter(hit_found_plasmidfinder != FALSE) |>
      arrange(start_position) |>
      select(
        plasmid_plasmidfinder,
        category_plasmidfinder,
        subcategory_plasmidfinder,
        nucleotide_code
      ) |>
      rename(
        `Plasmid` = plasmid_plasmidfinder,
        `Category` = category_plasmidfinder,
        `Subcategory` = subcategory_plasmidfinder,
        `Nucleotide code` = nucleotide_code
      )
    
    datatable(plasmidfinder_sequences,
              options = list(pageLength = 10))
  }
} else {
  cat('PlasmidFinder sequence data not available.')
}
```

Taxonomy
=====================================

Column {data-width=800}
-------------------------------------

### Taxonomic classification by GTDB-Tk classification (classify_wf)

```{r taxonomy_table}
gtdbtk_data <- safe_get('gtdbtk_classification_results')
if (!is.null(gtdbtk_data)) {

    gtdbtk_data_adapted <- gtdbtk_data |>
    select(user_genome_gtdbtk_classification, classification_gtdbtk_classification, closest_genome_reference_gtdbtk_classification, closest_genome_ani_gtdbtk_classification, closest_genome_af_gtdbtk_classification, closest_genome_taxonomy_gtdbtk_classification, classification_method_gtdbtk_classification, note_gtdbtk_classification, msa_percent_gtdbtk_classification, source_file_gtdbtk_classification) |>
    t()
  formatted_col_names <- c(
    'Analysed genome',
    'Assigned taxonomy',
    'Closest reference genome',
    'ANI with closest reference genome',
    'Alignment fraction with closest reference genome',
    'Taxonomy of closest reference genome',
    'Classification method',
    'Note regarding classification',            
    'MSA (%)',
    'Source file'
  )
  
  rownames(gtdbtk_data_adapted) <- formatted_col_names
  
  kable(gtdbtk_data_adapted,
        col.names = c('Value'))


} else {
  cat('GTDB-Tk classification results not available.')
}
```

Column {fig-width=400}
-------------------------------------

### Taxonomic tree by GTDB-Tk classification (classify_wf)

```{r taxonomy_plot, fig.width=15, fig.height=10}
tree_plot <- safe_get('gtdbtk_25_tree_plot')
if (!is.null(tree_plot)) {
  tree_plot_adapt <- tree_plot + 
    labs(title = NULL, subtitle = NULL) + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'))
  print(tree_plot_adapt)
} else {
  cat('Taxonomic tree not available.')
}
```

Gene Annotation
=====================================

Column {data-width=600}
-------------------------------------

### Genome annotation summary
```{r annotation_short_table}
annotation_data <- safe_get('annotation_results_integrated')
if (!is.null(annotation_data)) {
  # Check if the type_SC column exists
  if ('type_SC' %in% colnames(annotation_data)) {
    summary_table <- data.frame(
      Category = c('Total number of genes', 'Coding sequences (CDS)', 'Non-coding sequences', 'Uncategorised sequences', 'Hypothetical protein percentage of CDS'),
      Count = c(
        nrow(annotation_data),
        sum(annotation_data$type_SC == 'CDS', na.rm = TRUE),
        sum(annotation_data$type_SC != 'CDS', na.rm = TRUE),
        sum(is.na(annotation_data$type_SC)),
        (sum(annotation_data$gene_product_SC == 'hypothetical protein', na.rm = TRUE) / sum(annotation_data$type_SC == 'CDS', na.rm = TRUE))*100
      )
    )
    summary_table$Count <- format(summary_table$Count, scientific = FALSE)
    summary_table$Count <- ifelse(summary_table$Category == 'Hypothetical protein percentage', round(as.numeric(summary_table$Count), 2), round(as.numeric(summary_table$Count)))
    kable(summary_table)

  } else {
  cat('Type results not available.')
  }
} else {
  cat('Integrated annotation results not available.')
}
```

### Gene (product) confidence

```{r annotation_confidence_plot}
if (!is.null(annotation_data)) {
  confidence_levels <- c('very high', 'high', 'medium', 'low', 'very low', 'unknown')
  confidence_abbr <- c('vH', 'H', 'M', 'L', 'vL', 'U')
  
  color_mapping <- c(
    'very high' = '#336392',
    'high' = '#72A6CE',
    'medium' = '#AAC0CF',
    'low' = '#E8C09C',
    'very low' = '#F3A360',
    'unknown' = '#D0D1CA'
  )
  
  confidence_counts <- annotation_data %>%
    mutate(overall_gene_SC_confidence = case_when(
      overall_gene_SC_confidence == 'very_high' ~ 'very high',
      overall_gene_SC_confidence == 'very_low' ~ 'very low',
      TRUE ~ overall_gene_SC_confidence
    )) %>%
    count(overall_gene_SC_confidence)
  
  confidence_counts$overall_gene_SC_confidence <- factor(
    confidence_counts$overall_gene_SC_confidence,
    levels = confidence_levels,
    ordered = TRUE
  )
  
  confidence_counts <- confidence_counts %>%
    arrange(overall_gene_SC_confidence)
  
  confidence_labels <- paste0(confidence_counts$overall_gene_SC_confidence, ' (', 
                              confidence_abbr[match(confidence_counts$overall_gene_SC_confidence, confidence_levels)], ')')
  
  plot_ly(confidence_counts, 
          labels = ~confidence_labels, 
          values = ~n, 
          type = 'pie',
          marker = list(colors = color_mapping[as.character(confidence_counts$overall_gene_SC_confidence)]),
          sort = FALSE,
          textinfo = 'label+percent',
          insidetextorientation = 'radial') %>%
    layout(title = 'Gene confidence') |>
    set_font()
} else {
  cat('Gene product confidence data not available.')
}
```

Column {.tabset data-width=1000}
-------------------------------------

### Gene table with color-coded confidence

```{r annotation_table}
if (!is.null(annotation_data)) {
  # Function to map factor levels to colors and abbreviations with a smaller font size for abbreviations
  get_confidence_color_and_abbr <- function(confidence_level) {
    color <- switch(confidence_level,
        'very high' = '#336392',
        'high' = '#72A6CE',
        'medium' = '#AAC0CF',
        'low' = '#E8C09C',
        'very low' = '#F3A360',
        'unknown' = '#D0D1CA',
        '#D0D1CA')
    abbr <- switch(confidence_level,
        'very high' = 'vH',
        'high' = 'H',
        'medium' = 'M',
        'low' = 'L',
        'very low' = 'vL',
        'unknown' = 'U',
        'U')
    return(list(color = color, abbr = abbr))
  }
  
  # Modify the table to include circles with abbreviations in the first column
  gene_table <- annotation_data %>%
    arrange(start_position) %>%
    select(gene_SC, type_SC, gene_product_SC, start_position, end_position, length_bp, strand, contig_tag_fasta, EC_number_SC_best, COG_number_SC, overall_gene_SC_confidence) %>%
    mutate(
      overall_gene_SC_confidence = case_when(
        overall_gene_SC_confidence == 'very_high' ~ 'very high',
        overall_gene_SC_confidence == 'very_low' ~ 'very low',
        TRUE ~ overall_gene_SC_confidence
      ),
      # Color and abbreviation circles in place of row numbers
      row_color_circle = mapply(function(conf) {
        color_abbr <- get_confidence_color_and_abbr(conf)
        sprintf('<span style=\"display:inline-block; width:20px; height:20px; border-radius:50%%; background-color:%s; text-align:center; line-height:20px; color:white; font-size:10px;\">%s</span>', 
                color_abbr$color, color_abbr$abbr)
      }, overall_gene_SC_confidence)
    ) %>%
    rename(
      `Confidence` = row_color_circle,
      `Gene` = gene_SC,
      `Type` = type_SC,
      `Gene product` = gene_product_SC,
      `Start position` = start_position,
      `End position` = end_position,
      `Length (bp)` = length_bp,
      `Strand` = strand,
      `Contig` = contig_tag_fasta,
      `EC number (best hit)` = EC_number_SC_best,
      `COG number` = COG_number_SC
    ) %>%
    select(Confidence, everything(), -overall_gene_SC_confidence)
  
  # Display the table, inserting the color-coded circles and removing default row numbers
  datatable(gene_table,
            escape = FALSE,  # Allow rendering of HTML in the table
            rownames = FALSE,  # Disable row numbers
            options = list(pageLength = 20, 
                           columnDefs = list(list(className = 'dt-left', targets = '_all'))))
} else {
  cat('Gene annotation data not available.')
}

```

### Nucleotide and amino acid code

```{r annotation_table_nc}
if (!is.null(annotation_data)) {
  # Function to map factor levels to colors and abbreviations with a smaller font size for abbreviations
  get_confidence_color_and_abbr <- function(confidence_level) {
    color <- switch(confidence_level,
        'very high' = '#336392',
        'high' = '#72A6CE',
        'medium' = '#AAC0CF',
        'low' = '#E8C09C',
        'very low' = '#F3A360',
        'unknown' = '#D0D1CA',
        '#D0D1CA')
    abbr <- switch(confidence_level,
        'very high' = 'vH',
        'high' = 'H',
        'medium' = 'M',
        'low' = 'L',
        'very low' = 'vL',
        'unknown' = 'U',
        'U')
    return(list(color = color, abbr = abbr))
  }
  
  # Modify the table to include circles with abbreviations in the first column
  gene_table2 <- annotation_data %>%
    arrange(start_position) %>%
    select(gene_SC, gene_product_SC, EC_number_SC_best, COG_number_SC, nucleotide_code, amino_acid_code, overall_gene_SC_confidence) %>%
    mutate(
      overall_gene_SC_confidence = case_when(
        overall_gene_SC_confidence == 'very_high' ~ 'very high',
        overall_gene_SC_confidence == 'very_low' ~ 'very low',
        TRUE ~ overall_gene_SC_confidence
      ),
      # Color and abbreviation circles in place of row numbers
      row_color_circle = mapply(function(conf) {
        color_abbr <- get_confidence_color_and_abbr(conf)
        sprintf('<span style=\"display:inline-block; width:20px; height:20px; border-radius:50%%; background-color:%s; text-align:center; line-height:20px; color:white; font-size:10px;\">%s</span>', 
                color_abbr$color, color_abbr$abbr)
      }, overall_gene_SC_confidence)
    ) %>%
    rename(
      `Confidence` = row_color_circle,
      `Gene` = gene_SC,
      `Nucleotide code` = nucleotide_code,
      `Amino acid code` = amino_acid_code
    )
  
  # Display the table, inserting the color-coded circles and removing default row numbers
  datatable(gene_table2 %>% select(Confidence, everything(), -overall_gene_SC_confidence),
            escape = FALSE,  # Allow rendering of HTML in the table
            rownames = FALSE,  # Disable row numbers
            options = list(pageLength = 20, 
                           columnDefs = list(list(className = 'dt-left', targets = '_all'))))
} else {
  cat('Gene annotation data not available.')
}
```

Column {data-width=1000}
-------------------------------------

### KEGG modules (MicrobeAnnotator)

```{r KEEG_plot}
KEGG_modules <- safe_get('microbeannotator_KEGG_modules_results')
if (!is.null(KEGG_modules)) {

# Filter out rows with completeness = 0
df_filtered <- microbeannotator_KEGG_modules_results %>%
  filter(KEGG_module_completeness_microbeannotator > 0)

# Count occurrences of each KEGG module group and order x-axis by group size
group_counts <- df_filtered %>%
  count(KEGG_module_group_microbeannotator, sort = TRUE)

# Maintain y-axis ordering based on `ordering` column, but order x-axis by group size
df_ordered <- df_filtered %>%
  group_by(KEGG_module_group_microbeannotator) %>%
  mutate(ordering = sum(KEGG_module_completeness_microbeannotator)) %>%
  ungroup() %>%
  arrange(KEGG_module_group_microbeannotator, M_number_microbeannotator) %>%  
  mutate(KEGG_module_group_microbeannotator = factor(
    KEGG_module_group_microbeannotator, 
    levels = group_counts$KEGG_module_group_microbeannotator  
  ))

y_ordered <- df_ordered %>%
  arrange(ordering) %>%
  pull(KEGG_module_group_microbeannotator) %>%
  unique()

x_ordered <- df_ordered %>%
  arrange(KEGG_module_group_microbeannotator) %>%
  pull(M_number_microbeannotator)

# Create hover text
hover_text <- paste(
  'M number: ', df_ordered$M_number_microbeannotator, '<br>Module: ', 
  df_ordered$KEGG_module_microbeannotator, '<br>Completeness: ', 
  df_ordered$KEGG_module_completeness_microbeannotator, '<br>Group: ', 
  df_ordered$KEGG_module_group_microbeannotator, sep = ''
)

# Define custom color scale
custom_colorscale <- list(
  c(0, '#336392'),  # Start (low values)
  c(0.5, '#D0D1CA'),  # Middle
  c(1, '#F3A360')  # End (high values)
)

# Generate the heatmap with reordered x and y axis
plot_ly(
  z = as.numeric(df_ordered$KEGG_module_completeness_microbeannotator),  
  x = factor(df_ordered$M_number_microbeannotator, levels = x_ordered), 
  y = factor(df_ordered$KEGG_module_group_microbeannotator, levels = y_ordered), 
  type = 'heatmap', 
  colorscale = custom_colorscale,
  colorbar = list(title = list(text = 'Completeness')),
  hoverinfo = 'text',
  text = hover_text
) %>%
  layout(
    xaxis = list(title = list(text = 'M number')),
    yaxis = list(title = list(text = 'KEGG module group'), categoryorder = 'array', categoryarray = y_ordered),
    showlegend = TRUE
    )

} else {
  cat('KEGG module data (MicrobeAnnotator) not available.')
}
```

### COG categories (Bakta)

```{r COG_plot}
if (!is.null(annotation_data)) {
# Define functional groups and color mapping
cog_labels <- data.frame(
  COG_category_SC = c('J', 'A', 'K', 'L', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 'X', 
                      'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S', 'NA'),  # Add 'NA'
  COG_description = c('Translation, ribosomal structure, and biogenesis', 
                      'RNA processing and modification',
                      'Transcription', 
                      'Replication, recombination, and repair', 
                      'Cell cycle control, cell division, and chromosome partitioning', 
                      'Nuclear structure', 
                      'Defense mechanisms', 
                      'Signal transduction mechanisms', 
                      'Cell wall/membrane/envelope biogenesis', 
                      'Cell motility', 
                      'Cytoskeleton', 
                      'Extracellular structures', 
                      'Intracellular trafficking, secretion, and vesicular transport', 
                      'Posttranslational modification, protein turnover, chaperones',
                      'Mobilome: prophages, transposons',
                      'Energy production and conversion', 
                      'Carbohydrate transport and metabolism', 
                      'Amino acid transport and metabolism', 
                      'Nucleotide transport and metabolism', 
                      'Coenzyme transport and metabolism', 
                      'Lipid transport and metabolism', 
                      'Inorganic ion transport and metabolism', 
                      'Secondary metabolites biosynthesis, transport, and catabolism', 
                      'General function prediction only', 
                      'Function unknown', 
                      'Genes without assigned COG category'),
  Functional_Group = c(rep('Information storage and processing', 4),
                       rep('Cellular processes and signaling', 11),
                       rep('Metabolism', 8),
                       rep('Poorly characterized', 3))  # NA under Poorly Characterized
)

# Define colors for each functional group
color_mapping <- c(
  'Information storage and processing' = '#336392',
  'Cellular processes and signaling' = '#934A55',
  'Metabolism' = '#31837F',
  'Poorly characterized' = 'grey'
)

# Define the desired order of categories
category_order <- c('J', 'A', 'K', 'L', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 'X', 
                     'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S', 'NA')

# Split multi-letter COG categories into individual categories
annotation_data_adapted <- annotation_data %>%
  separate_rows(COG_category_SC, sep = '')

# Replace NA values with 'NA' to represent unassigned COG category
annotation_data_adapted <- annotation_data_adapted %>%
  mutate(COG_category_SC = ifelse(is.na(COG_category_SC), 'NA', COG_category_SC))

# Check if there are any COG categories in the data not defined in cog_labels
missing_categories <- setdiff(unique(annotation_data_adapted$COG_category_SC), cog_labels$COG_category_SC)

# Count occurrences of each COG category
cog_counts <- annotation_data_adapted %>%
  count(COG_category_SC) %>%
  filter(COG_category_SC %in% cog_labels$COG_category_SC) %>%
  arrange(desc(n))

# Merge with the full description and functional group for hover text and coloring
cog_counts <- cog_counts %>%
  left_join(cog_labels, by = 'COG_category_SC')

# Convert COG_category_SC to a factor with levels in the desired order
cog_counts$COG_category_SC <- factor(cog_counts$COG_category_SC, levels = category_order)

# Create a bar plot with the COG category letter as the x-axis, full description on hover, and color by functional group
plot_ly(cog_counts, 
        x = ~COG_category_SC, 
        y = ~n, 
        type = 'bar', 
        hoverinfo = 'text', 
        text = ~paste('Functional group:', Functional_Group, '<br>COG category:', COG_description, '<br>Gene count:', n),
        marker = list(color = ~color_mapping[Functional_Group])) %>%
  layout(xaxis = list(title = 'COG category', tickangle = -45),
         yaxis = list(title = 'Gene count'),
         margin = list(b = 100))
} else {
  cat('COG category data not available.')
}
```

Mobile Genetic Elements
=====================================

Column {.tabset data-width=400}
-------------------------------------

### Identified genomic islands by IslandPath

```{r islandpath_table}
islandpath_data <- safe_get('islandpath_results')
if (!is.null(islandpath_data)) {
    islandpath_data_adapted <- islandpath_data |>
    select(contig_tag_fasta, type_islandpath, start_position, end_position, length_bp, strand, source_file_islandpath) |>
    rename(
            `Contig` = contig_tag_fasta,
            `Type` = type_islandpath,
            `Start position` = start_position,
            `End position` = end_position,
            `Length (bp)` = length_bp,
            `Strand` = strand,
            `Source file` = source_file_islandpath
          )

  datatable(islandpath_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('IslandPath results not available.')
}
```

### Nucleotide code

```{r islandpath_table_nc}
islandpath_data <- safe_get('islandpath_results')
if (!is.null(islandpath_data)) {
    islandpath_data_adapted <- islandpath_data |>
    select(contig_tag_fasta, nucleotide_code) |>
    rename(
            `Contig` = contig_tag_fasta,
            `Nucleotide code` = nucleotide_code
          )

  datatable(islandpath_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('IslandPath results not available.')
}
```

  Column {.tabset data-width=400}
-------------------------------------
  
### Identified insertion sequence (IS) element by ISEScan
  
```{r isescan_table}
isescan_data <- safe_get('isescan_results')
if (!is.null(isescan_data)) {
  isescan_data_adapted <- isescan_data |>
    select(
      contig_tag_fasta, 
      tpase_cluster, 
      insertion_sequence_family, 
      insertion_sequence_start_position, 
      insertion_sequence_end_position, 
      insertion_sequence_length_bp, 
      strand, 
      insertion_sequence_copy_number, 
      orf_start_position, 
      orf_end_position, 
      orf_length_bp, 
      source_file_tsv
    ) |>
    dplyr::rename(
      `Contig` = contig_tag_fasta,
      `Transposase cluster` = tpase_cluster,
      `IS family` = insertion_sequence_family,
      `IS start position` = insertion_sequence_start_position,
      `IS end position` = insertion_sequence_end_position,
      `IS length (bp)` = insertion_sequence_length_bp,
      `Strand` = strand,
      `IS copy numbers` = insertion_sequence_copy_number,
      `ORF start position` = orf_start_position,
      `ORF end position` = orf_end_position,
      `ORF length (bp)` = orf_length_bp,
      `Source file` = source_file_tsv
    ) 

  datatable(isescan_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ISEScan results not available.')
}
```

### Inverted repeats (IR)

```{r isescan_table_ir}
isescan_data <- safe_get('isescan_results')
if (!is.null(isescan_data)) {
  isescan_data_adapted <- isescan_data |>
    select(
      contig_tag_fasta, 
      first_inverted_repeat_start_position,
      first_inverted_repeat_end_position,
      second_inverted_repeat_start_position,
      second_inverted_repeat_end_position,
      inverted_repeat_length_bp,
      inverted_repeat_identity_matches_bp,
      inverted_repeat_alignment_gaps_number,
      terminal_inverted_repeat_sequence
    ) |>
    rename(
      `Contig` = contig_tag_fasta,,
      `1. IR start position` = first_inverted_repeat_start_position,
      `1. IR end position` = first_inverted_repeat_end_position,
      `2. IR start position` = second_inverted_repeat_start_position,
      `2. IR end position` = second_inverted_repeat_end_position,
      `IR length (bp)` = inverted_repeat_length_bp,
      `IR nucleotide match (1. vs 2.)` = inverted_repeat_identity_matches_bp,
      `IR alignment gaps number` = inverted_repeat_alignment_gaps_number,
      `Terminal IR sequences` = terminal_inverted_repeat_sequence
    )
  
  datatable(isescan_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ISEScan results not available.')
}
```

### Nucleotide code (IS)

```{r isescan_table_is_nc}
isescan_data <- safe_get('isescan_results')
if (!is.null(isescan_data)) {
  isescan_data_adapted <- isescan_data |>
    select(
      contig_tag_fasta, 
      insertion_sequence_start_position,
      insertion_sequence_end_position,
      insertion_sequence_nucleotide_code
    ) |>
    rename(
      `Contig` = contig_tag_fasta,,
      `IS start position` = insertion_sequence_start_position,
      `IS end position` = insertion_sequence_end_position,
      `IS nucleotide code` = insertion_sequence_nucleotide_code
    )
  
  datatable(isescan_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ISEScan results not available.')
}
```

### Nucleotide code (ORF)

```{r isescan_table_orf_nc}
isescan_data <- safe_get('isescan_results')
if (!is.null(isescan_data)) {
  isescan_data_adapted <- isescan_data |>
    select(
      contig_tag_fasta, 
      orf_start_position,
      orf_end_position,
      orf_nucleotide_code
    ) |>
    rename(
      `Contig` = contig_tag_fasta,,
      `ORF start position` = orf_start_position,
      `ORF end position` = orf_end_position,
      `ORF nucleotide code` = orf_nucleotide_code
    )
  
  datatable(isescan_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ISEScan results not available.')
}
```

### Amino acid code (ORF)

```{r isescan_table_orf_aa}
isescan_data <- safe_get('isescan_results')
if (!is.null(isescan_data)) {
  isescan_data_adapted <- isescan_data |>
    select(
      contig_tag_fasta, 
      orf_start_position,
      orf_end_position,
      orf_amino_acid_code
    ) |>
    rename(
      `Contig` = contig_tag_fasta,,
      `ORF start position` = orf_start_position,
      `ORF end position` = orf_end_position,
      `ORF amino acid code` = orf_amino_acid_code
    )
  
  datatable(isescan_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ISEScan results not available.')
}
```

  Column {.tabset data-width=400}
-------------------------------------
  
### Identified viral sequences by VirSorter2
  
```{r virsorter2_table}
virsorter2_data <- safe_get('virsorter2_results')
if (!is.null(virsorter2_data)) {
  virsorter2_data_adapted <- virsorter2_data |>
    select(contig_tag_fasta, contig_tag_virsorter2, start_position, end_position, length_bp, score_virsorter2, max_score_group_virsorter2, hallmark_gene_count, viral, bacterial, archeal, eukaryotic, mixed, unaligned, source_file_virsorter2_tsv1, source_file_virsorter2_tsv2, source_file_virsorter2_fa) |>
    dplyr::rename(
      `Contig` = contig_tag_fasta,
      `Identification category` = contig_tag_virsorter2,
      `Start position` = start_position,
      `End position` = end_position,
      `Length (bp)` = length_bp,
      `VirSorter2 score` = score_virsorter2,
      `Type (best fit)` = max_score_group_virsorter2,
      `Hallmark gene count` = hallmark_gene_count,
      `viral gene %` = viral,
      `bacterial gene %` = bacterial,
      `archeal gene %` = archeal,
      `eukaryotic gene %` = eukaryotic,
      `mixed gene %` = mixed,
      `unaligned gene %` = unaligned,
      `Source file` = source_file_virsorter2_tsv1,
      `Source file2` = source_file_virsorter2_tsv2,
      `Source file3` = source_file_virsorter2_fa
    ) %>%
    mutate(`Identification category` = sub('.*_', '', `Identification category`))
  
  datatable(virsorter2_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('VirSorter2 results not available.')
}
```

### Nucleotide code

```{r virsorter2_table_nc}
virsorter2_data <- safe_get('virsorter2_results')
if (!is.null(virsorter2_data)) {
  virsorter2_data_adapted <- virsorter2_data |>
    select(contig_tag_fasta, start_position, end_position, nucleotide_code) |>
    rename(
      `Contig` = contig_tag_fasta,
      `Start position` = start_position,
      `End position` = end_position,
      `Nucleotide code` = nucleotide_code
    )
  
  datatable(virsorter2_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('VirSorter2 results not available.')
}
```

Functional Analyses
=====================================

Column {.tabset data-width=400}
-------------------------------------

### CAZyme analysis results dbCAN3

```{r dbcan3_table}
# Assuming safe_get is a function that retrieves the data
dbcan3_data <- safe_get('dbcan3_results')
annotation_data <- safe_get('annotation_results_integrated')

if (!is.null(dbcan3_data)) {
  merged <- FALSE
  if (!is.null(annotation_data)) {
      annotation_data <- annotation_data %>%
        mutate(
          start_position = as.numeric(start_position),
          end_position = as.numeric(end_position))

      dbcan3_data_adapted <- dbcan3_data %>%
        arrange(start_position) %>%
        select(CAZyme_gene_clusters_dbcan3, gene_type_dbcan3, protein_family_dbcan3, start_position, end_position, strand, contig_tag_fasta, nucleotide_code, source_file_dbcan3_out, source_file_dbcan3_txt) |>
        mutate(
          start_position = as.numeric(start_position),
          end_position = as.numeric(end_position)
        ) %>%
        left_join(annotation_data, by = c('start_position', 'end_position', 'strand', 'nucleotide_code', 'contig_tag_fasta'))
      
      merged <- TRUE
  }
  
  if (!merged) {
    dbcan3_data_adapted <- dbcan3_data %>%
      arrange(start_position) %>%
      select(CAZyme_gene_clusters_dbcan3, gene_type_dbcan3, protein_family_dbcan3, start_position, end_position, strand, source_file_dbcan3_out, source_file_dbcan3_txt)
  }

  # Rename columns based on whether merge occurred or not
  if (merged) {
    dbcan3_data_adapted <- dbcan3_data_adapted %>%
      select(CAZyme_gene_clusters_dbcan3, gene_SC, gene_product_SC, gene_type_dbcan3, protein_family_dbcan3, start_position, end_position, length_bp, strand, contig_tag_fasta, COG_number_SC, EC_number_SC_best, source_file_dbcan3_out, source_file_dbcan3_txt) %>%
      rename(
        `CAZyme gene cluster tag` = CAZyme_gene_clusters_dbcan3,
        `Gene` = gene_SC,
        `Gene product` = gene_product_SC,
        `Gene type` = gene_type_dbcan3,
        `Protein family` = protein_family_dbcan3,
        `Start position` = start_position,
        `End position` = end_position,
        `Length (bp)` = length_bp,
        `Strand` = strand,
        `Contig` = contig_tag_fasta,
        `COG number` = COG_number_SC,
        `EC number (best hit)` = EC_number_SC_best,
        `Source file 1` = source_file_dbcan3_out,
        `Source file 2` = source_file_dbcan3_txt
      )
  } else {
    dbcan3_data_adapted <- dbcan3_data_adapted %>%
      select(CAZyme_gene_clusters_dbcan3, gene_type_dbcan3, protein_family_dbcan3, start_position, end_position, length_bp, strand, contig_tag_fasta, source_file_dbcan3_out, source_file_dbcan3_txt) %>%
      mutate(length_bp = end_position - start_position + 1) %>%
      rename(
        `CAZyme gene cluster tag` = CAZyme_gene_clusters_dbcan3,
        `Gene type` = gene_type_dbcan3,
        `Protein family` = protein_family_dbcan3,
        `Start position` = start_position,
        `End position` = end_position,
        `Length (bp)` = length_bp,
        `Strand` = strand,
        `Contig` = contig_tag_fasta,
        `Source file 1` = source_file_dbcan3_out,
        `Source file 2` = source_file_dbcan3_txt
      )
  }

  datatable(dbcan3_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('dbCAN3 results not available.')
}
```

### Nucleotide and amino acid code 

```{r dbcan3_table_nc}
if (!is.null(dbcan3_data)) {
  dbcan3_data_adapted <- dbcan3_data |>
    arrange(start_position) |>
    select(contig_tag_fasta, CAZyme_gene_clusters_dbcan3, gene_type_dbcan3, protein_family_dbcan3, start_position, end_position, strand, nucleotide_code, source_file_dbcan3_out, source_file_dbcan3_txt) |>
    mutate(
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position)
    )

  merged <- FALSE
  if (!is.null(annotation_data)) {
      annotation_data_adapted <- annotation_data |>
        arrange(start_position) |>
        select(contig_tag_fasta, gene_SC, gene_product_SC, start_position, end_position, strand, length_bp, COG_number_SC, EC_number_SC_best, nucleotide_code, amino_acid_code) |>
        mutate(
          start_position = as.numeric(start_position),
          end_position = as.numeric(end_position)
        )

      dbcan3_data_adapted <- dbcan3_data_adapted |>
        left_join(annotation_data_adapted, by = c('start_position', 'end_position', 'strand', 'nucleotide_code', 'contig_tag_fasta'))
      
      merged <- TRUE
    
  }

  if (merged) {
    dbcan3_data_adapted <- dbcan3_data_adapted |>
      select(CAZyme_gene_clusters_dbcan3,
             gene_SC,
             gene_product_SC,
             gene_type_dbcan3,
             protein_family_dbcan3,
             nucleotide_code,
             amino_acid_code
      ) |>
      rename(
        `CAZyme gene cluster tag` = CAZyme_gene_clusters_dbcan3,
        `Gene` = gene_SC,
        `Gene product` = gene_product_SC,
        `Gene type` = gene_type_dbcan3,
        `Protein family` = protein_family_dbcan3,
        `Nucleotide code` = nucleotide_code,
        `Amino acid code` = amino_acid_code
      )
  } else {
    dbcan3_data_adapted <- dbcan3_data_adapted |>
      select(CAZyme_gene_clusters_dbcan3,
             gene_type_dbcan3,
             protein_family_dbcan3,
             nucleotide_code
      ) |>
      rename(
        `CAZyme gene cluster tag` = CAZyme_gene_clusters_dbcan3,
        `Gene type` = gene_type_dbcan3,
        `Protein family` = protein_family_dbcan3,
        `Nucleotide code` = nucleotide_code
      )
  }

  datatable(dbcan3_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('dbCAN3 nucleotide code not available.')
}
```

Column {.tabset data-width=400}
-------------------------------------

### AMR phenotype profile by ResFinder

```{r resfinder_amr_profile_table}
resfinder_profile <- safe_get('resfinder_AMR_profile')
if (!is.null(resfinder_profile)) {
  resfinder_profile_adapted <- resfinder_profile %>%
  mutate(genomic_AMR_phenotype_resfinder = factor(
    genomic_AMR_phenotype_resfinder,
    levels = c(\"Resistant\", setdiff(unique(genomic_AMR_phenotype_resfinder), \"Resistant\"))
    )) |>
    arrange(genomic_AMR_phenotype_resfinder) |>
    select(genomic_AMR_phenotype_resfinder, antimicrobial_agent_resfinder, antimicrobial_agent_class_resfinder, AMR_phenotype_genetic_background_resfinder, source_file_AMR_profile_resfinder) |>
    rename(
      `Resistance status` = genomic_AMR_phenotype_resfinder,
      `Antimicrobial agent` = antimicrobial_agent_resfinder,
      `Antimicrobial agent class` = antimicrobial_agent_class_resfinder,
      `Genetic background` = AMR_phenotype_genetic_background_resfinder,
      `Source file` = source_file_AMR_profile_resfinder
    )
  datatable(resfinder_profile_adapted,
            options = list(pageLength = 10))
} else {
  cat('AMR phenotype profile not available.')
}
```

### Identified AMR genes by ResFinder

```{r resfinder_table}
resfinder_data <- safe_get('resfinder_results')
if (!is.null(resfinder_data)) {
  resfinder_data_adapted <- resfinder_data |>
    select(AMR_gene_resfinder, AMR_phenotype_resfinder, start_position, end_position, nucleotide_identity_resfinder, percentage_length_of_reference_sequence_resfinder, contig_tag_fasta, source_file_results_resfinder) |>
    mutate(percentage_length_of_reference_sequence_resfinder = round(percentage_length_of_reference_sequence_resfinder, 2)) %>%
    rename(
            `Gene` = AMR_gene_resfinder,
            `AMR phenotype` = AMR_phenotype_resfinder,
            `Start position` = start_position,
            `End position` = end_position,
            `Nucleotide identity` = nucleotide_identity_resfinder,
            `Percentage length of reference` = percentage_length_of_reference_sequence_resfinder,
            `Contig tag` = contig_tag_fasta,
            `Source file` = source_file_results_resfinder
          )

  datatable(resfinder_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ResFinder results not available.')
}
```

### Nucleotide and amino acid code (ResFinder)

```{r resfinder_table_nc}
resfinder_data <- safe_get('resfinder_results')
annotation_data <- safe_get('annotation_results_integrated')

if (!is.null(resfinder_data) & !is.null(annotation_data)) {
  resfinder_data_adapted <- resfinder_data |>
    arrange(start_position) |>
    select(AMR_gene_resfinder, AMR_phenotype_resfinder, start_position, end_position)

  annotation_data_adapted <- annotation_data |>
    arrange(start_position) |>
    select(start_position, end_position, nucleotide_code, amino_acid_code) |>
    mutate(
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position)
    )
  
  resfinder_data_adapted <- resfinder_data_adapted |>
    mutate(
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position)
    ) |>
    left_join(annotation_data_adapted, by = c('start_position', 'end_position')) |>
    select(AMR_gene_resfinder, 
           AMR_phenotype_resfinder, 
           nucleotide_code,
           amino_acid_code
    ) |>
    rename(
      `AMR gene` = AMR_gene_resfinder,
      `AMR phenotype` = AMR_phenotype_resfinder,
      `Nucleotide code` = nucleotide_code,
      `Amino acid code` = amino_acid_code
    )

  datatable(resfinder_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('ResFinder nucleotide and amino acid codes not available.')
}
```

Column {.tabset data-width=400}
-------------------------------------

### Identified AMR genes by AMRFinderPlus

```{r amrfinderplus_table}
amrfinder_data <- safe_get('amrfinderplus_results')

if (!is.null(amrfinder_data)) {
  amrfinder_data_adapted <- amrfinder_data |>
    select(
      gene_symbol,
      sequence_name,
      resistance_class,
      resistance_subclass,
      element_type,
      element_subtype,
      start_position,
      end_position,
      strand,
      identity_percentage,
      coverage_percentage,
      contig_tag_fasta,
      source_file
    ) |>
    rename(
      `Gene` = gene_symbol,
      `Gene name` = sequence_name,
      `Resistance class` = resistance_class,
      `Resistance subclass` = resistance_subclass,
      `Element type` = element_type,
      `Element subtype` = element_subtype,
      `Start position` = start_position,
      `End position` = end_position,
      `Strand` = strand,
      `Identity (%)` = identity_percentage,
      `Coverage (%)` = coverage_percentage,
      `Contig tag` = contig_tag_fasta,
      `Source file` = source_file
    )
    
  datatable(amrfinder_data_adapted,
            options = list(pageLength = 10))
} else {
  cat('AMRFinderPlus results not available.')
}
```

### Nucleotide code and mutations (AMRFinderPlus)

```{r amrfinderplus_sequences}
amrfinder_data <- safe_get('amrfinderplus_results')

if (!is.null(amrfinder_data)) {
  amrfinder_sequences <- amrfinder_data |>
    arrange(start_position) |>
    select(
      gene_symbol,
      sequence_name,
      nucleotide_code,
      mutation_details,
      resistance_class
    ) |>
    rename(
      `Gene` = gene_symbol,
      `Gene name` = sequence_name,
      `Nucleotide code` = nucleotide_code,
      `Mutations` = mutation_details,
      `Resistance class` = resistance_class
    )
    
  datatable(amrfinder_sequences,
            options = list(pageLength = 10))
} else {
  cat('AMRFinderPlus sequence data not available.')
}
```

Genome Visualisation
=====================================
  
Column {data-width=900}
-------------------------------------
  
### Integrated genome visualisation
  
```{r circlize_plot, fig.width=8, fig.height=8, dpi=300}

# Check if fasta_file is available (it's obligatory)
fasta_file <- safe_get('selected_assembly_results')
if (!is.null(fasta_file)) {
  
  # Prepare circlize genome data
  circlize_genome <- data.frame(
    contig_tag_fasta = fasta_file$contig_tag_fasta,
    start_position = 1,
    stop_position = nchar(fasta_file$nucleotide_code),
    nucleid_acid_code = as.character(fasta_file$nucleotide_code)
  )
  
  prepare_circlize_data <- function(data_name, required_cols) {
    data <- safe_get(data_name)
    if (!is.null(data) && all(required_cols %in% colnames(data))) {
      result <- data %>%
        select(all_of(c('contig_tag_fasta', 'start_position', 'end_position', required_cols)))
      
      # Check if the strand column exists and has valid values
      if ('strand' %in% colnames(data)) {
        result <- result %>% mutate(strand = data$strand)
      } else {
        result <- result %>% mutate(strand = '?')  # Assign '?' if the strand column is missing
      }
      
      return(result)
    } else {
      message(paste(data_name, 'not found or missing required columns.'))
      return(NULL)
    }
  }
  
  prepare_circlize_data_isescan <- function(data_name, required_cols) {
    data <- safe_get(data_name)
    if (!is.null(data) && all(required_cols %in% colnames(data))) {
      result <- data %>%
        select(all_of(c('contig_tag_fasta', 'insertion_sequence_start_position', 'insertion_sequence_end_position', required_cols))) %>%
        mutate(
          start_position = insertion_sequence_start_position,
          end_position = insertion_sequence_end_position,
        ) %>%
        select(all_of(c('contig_tag_fasta', 'start_position', 'end_position', required_cols)))
      
      # Fixed strand handling
      if ('strand' %in% colnames(data)) {
        result <- result %>% 
          # Replace NA values with '?' in the strand column
          mutate(strand = ifelse(is.na(data$strand), '?', as.character(data$strand)))
      } else {
        result <- result %>% mutate(strand = '?')
      }
      
      return(result)
    } else {
      message(paste(data_name, 'not found or missing required columns.'))
      return(NULL)
    }
  }
  
  # Prepare (standardised) circlize data for each result type
  circlize_annotation <- prepare_circlize_data('annotation_results_integrated', c('gene_SC', 'gene_product_SC', 'strand'))
  circlize_islandpath <- prepare_circlize_data('islandpath_results', 'strand')
  circlize_isescan <- prepare_circlize_data_isescan('isescan_results', 'strand')
  circlize_virsorter2 <- prepare_circlize_data('virsorter2_results', character(0))
  circlize_dbcan3 <- prepare_circlize_data('dbcan3_results', 'strand')
  circlize_amrfinderplus <- prepare_circlize_data('amrfinderplus_results', 'strand')
  circlize_resfinder <- prepare_circlize_data('resfinder_results', character(0))
  #circlize_crisprcasfinder <- prepare_circlize_data('crisprcasfinder_results', 'strand')
  
  # Combine RGI and Resfinder results for AMR visualization
  circlize_AMR <- bind_rows(circlize_amrfinderplus, circlize_resfinder)
  
  # Circlize plot function
  circlize_plot <- function() {
    circos.clear()
    circos.par(start.degree = 90, cell.padding = c(0, 0, 0, 0))
    circos.genomicInitialize(circlize_genome)
    
    plot_sections <- function(genes, y_start, y_end, color) {
      if (!is.null(genes) && nrow(genes) > 0) {
        starts <- genes$start_position
        ends <- genes$end_position
        circos.rect(starts, rep(y_start, length(starts)), ends, rep(y_end, length(ends)), col = color, border = NA)
      }
    }
    
    # Origin of replication
    circos.track(ylim = c(0.25, 0.5), panel.fun = function(x, y) {
      if (!is.null(circlize_annotation)) {
        contig <- CELL_META$sector.index
        ori <- subset(circlize_annotation, contig_tag_fasta == contig & gene_product_SC == 'origin of replication')
        if(nrow(ori) > 0) {
          circos.points(ori$start_position, rep(0.5, nrow(ori)), pch = 18, col = '#A94322', cex = 1)
        }
      }
    }, bg.border = 0, track.height = 0.01)
    
    # Function to plot gene segments
    plot_gene_segments <- function(data, colors) {
      if (!is.null(data)) {
        circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
          contig <- CELL_META$sector.index
          genes_in_contig <- subset(data, contig_tag_fasta == contig)
          if(nrow(genes_in_contig) > 0) {
            plot_sections(genes_in_contig[genes_in_contig$strand == '+', ], 0.5, 1, colors[1])
            plot_sections(genes_in_contig[genes_in_contig$strand == '-', ], 0, 0.5, colors[2])
            if (any(genes_in_contig$strand == '?')) {
              plot_sections(genes_in_contig[genes_in_contig$strand == '?', ], 0.25, 0.75, colors[3])
            }
          }
        }, bg.border = 1, track.height = 0.1)
      }
    }
    
    # Plot gene segments for each data type
    plot_gene_segments(circlize_annotation, c('#699DC6', '#EA7E2E', '#4D4D4D'))
    plot_gene_segments(circlize_islandpath, c('#2B5C8A', '#9E3D22', '#4D4D4D'))
    plot_gene_segments(circlize_isescan, c('#699DC6', '#EA7E2E', '#4D4D4D'))
    plot_gene_segments(circlize_virsorter2, c('#2B5C8A', '#9E3D22', '#4D4D4D'))
    plot_gene_segments(circlize_dbcan3, c('#699DC6', '#EA7E2E', '#4D4D4D'))
    plot_gene_segments(circlize_AMR, c('#2B5C8A', '#9E3D22', '#4D4D4D'))
    #plot_gene_segments(circlize_crisprcasfinder, c('#699DC6', '#EA7E2E', '#4D4D4D'))
  }
  
  # Generate the plot
  circlize_plot()
  
} else {
  cat('Fasta file data not available.')
}

```

Column {data-width=300}
-------------------------------------
  
### Genome visualisation legend (oriented form outside to inside)
  
```{r circlize_plot_legend, fig.width=2, fig.height=6, dpi=300}
fasta_file <- safe_get('selected_assembly_results')
if (!is.null(fasta_file)) {
  
  plot_legend <- function() {
    # Initialize an empty list to store legends
    legend_list <- list()
    
    # Helper function to add legend if data is available
    add_legend_if_available <- function(data, legend) {
      if (!is.null(data) && nrow(data) > 0) {
        legend_list <<- c(legend_list, list(legend))
      }
    }
    
    # Origin of replication legend (always include if annotation data is available)
    if (!is.null(circlize_annotation)) {
      legend_list <- c(legend_list, list(
        Legend(at = c(''), type = 'points', 
               legend_gp = gpar(col = c('#A94322')), 
               title_position = 'topleft', title = 'Origin of replication', pch = 18)
      ))
    }
    
    # Genes legend
    add_legend_if_available(circlize_annotation, 
                            Legend(at = c('+ strand', '- strand'), type = 'points', 
                                   legend_gp = gpar(col = c('#699DC6', '#EA7E2E')), 
                                   title_position = 'topleft', title = 'Genes'))
    
    # Genomic Islands legend
    add_legend_if_available(circlize_islandpath, 
                            Legend(at = c('+ strand', '- strand'), type = 'points', 
                                   legend_gp = gpar(col = c('#2B5C8A', '#9E3D22')), 
                                   title_position = 'topleft', title = 'Genomic Islands'))
    
    # ISE legend
    add_legend_if_available(circlize_isescan, 
                            Legend(at = c('+ strand', '- strand'), type = 'points', 
                                   legend_gp = gpar(col = c('#699DC6', '#EA7E2E')), 
                                   title_position = 'topleft', title = 'Insertion sequence elements'))
    
    # Viral sequences legend
    if (!is.null(circlize_virsorter2) && nrow(circlize_virsorter2) > 0) {
      legend_list <- c(legend_list, list(
        Legend(at = c('+ strand', '- strand', 'unknown strand'), type = 'points', 
               legend_gp = gpar(col = c('#2B5C8A', '#9E3D22', '#4D4D4D')), 
               title_position = 'topleft', title = 'Viral sequences'))
      )
    }
    
    # CAZymes legend
    add_legend_if_available(circlize_dbcan3, 
                            Legend(at = c('+ strand', '- strand'), type = 'points', 
                                   legend_gp = gpar(col = c('#699DC6', '#EA7E2E')), 
                                   title_position = 'topleft', title = 'CAZymes'))
    
    # AMR legend
    if (!is.null(circlize_AMR) && nrow(circlize_AMR) > 0) {
      legend_list <- c(legend_list, list(
        Legend(at = c('+ strand', '- strand', 'unknown strand'), type = 'points', 
               legend_gp = gpar(col = c('#2B5C8A', '#9E3D22', '#4D4D4D')), 
               title_position = 'topleft', title = 'AMR'))
      )
    }
    
    # # CRISPRCas legend
    # add_legend_if_available(circlize_crisprcasfinder, 
    #                         Legend(at = c('+ strand', '- strand'), type = 'points', 
    #                                legend_gp = gpar(col = c('#699DC6', '#EA7E2E')), 
    #                                title_position = 'topleft', title = 'CRISPRCas sequences'))
    
    # Only create and plot the legend if there's at least one item
    if (length(legend_list) > 0) {
      # Pack legends into a single grob using do.call
      lgd_grob <- do.call(packLegend, legend_list)
      
      # Plot the legend
      grid.newpage()  # Clear the plot area
      grid.draw(lgd_grob)  # Use grid.draw to draw the grob
    } else {
      # If no data is available, print a message
      grid.newpage()
      grid.text('No data available for legend', x = 0.5, y = 0.5, gp = gpar(cex = 1.2))
    }
  }
  
  # Call the function to plot the legend
  plot_legend()
  
} else {
  cat('Fasta file data not available.')
}
```

"

return(rmd_content)
}

# Main script
main <- function() {
  # Generate RMarkdown content
  rmd_content <- generate_rmarkdown()
  
  # File names
  rmd_file <- file.path(output_dir, paste("StrainCascade_Analysis_Report_", sample_name, ".Rmd", sep = ""))
  html_file <- file.path(output_dir, paste("StrainCascade_Analysis_Report_", sample_name, ".html", sep = ""))
  
  # Write RMarkdown content to file
  writeLines(rmd_content, rmd_file)
  
  # Render the RMarkdown document
  tryCatch({
    rmarkdown::render(rmd_file, 
                      output_format = "flexdashboard::flex_dashboard", 
                      output_file = html_file,
                      params = list(rdata_file = rdata_file))
    cat("Interactive Markdown summary has been generated as '", basename(html_file), "'\n", sep = "")
  }, error = function(e) {
    cat("Error occurred while rendering the document:\n", conditionMessage(e), "\n")
    cat("Please check if all required libraries are installed and data files are present.\n")
  })
}

# Run the main function
main()