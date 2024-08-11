# BIOMES_WGseq_shiny_app.r - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# IMPORTANT: Before running this script, please make sure to open the corresponding Rproject 
# (*sample_name*_BIOMES_WGseq_Rproject.Rproj) in RStudio. This will set the working directory 
# to the correct location and ensure that the necessary data files are found.

# List of required packages
packages <- c("shiny", "DT", "tibble", "plotly", "dplyr")

# Function to check and install packages if necessary
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, lib = .libPaths()[1], dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Check and install required packages
sapply(packages, check_and_install)

# Load necessary libraries
library(shiny)
library(DT)
library(tibble)
library(plotly)
library(dplyr)

# List files matching the pattern in your working directory
RData_dir <- list.files(pattern = "*_BIOMES_WGseq.RData")

# Check if there is exactly one matching file
if (length(RData_dir) != 1) {
  stop("Expected exactly one *_BIOMES_WGseq.RData file in the working directory, but found ", length(RData_dir))
}

# Load the saved environment
load(RData_dir)

# Define UI for application
ui <- fluidPage(
  tags$head(
    tags$title("BIOMES WGseq summary"),
    tags$style(HTML(".title h1 { margin-bottom: 40px; } 
                    .subtitle h3 { margin-top: -25px; margin-bottom: 25px; font-style: italic; font-size: 16px; }
                     body {font-family: 'Helvetica', sans-serif}
                    "))
  ),
  tags$div(
    class = "title",
    tags$h1("BIOMES WGseq summary")
  ),
  tags$div(
    class = "subtitle",
    tags$h3(paste("Organism: ", gtdbtk_organism_name, " / ", "Sequencing file: ", sequencing_file, sep = ""))
  ),
  sidebarLayout(
    sidebarPanel(
      width = 2,  # Set the width by x to x/12 of the total width (screen width)
      selectInput("analysis_steps", "Analysis steps:", choices = c("Genome assembly", "Taxonomic classification", "Genome annotation", "Metabolic functions")),
      uiOutput("var_ui"),
      uiOutput("description_ui"),  # Add a UI output for the description dropdown menu
    ),
    mainPanel(
      width = 10,
      uiOutput("ui")
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$ui <- renderUI({
    if (input$analysis_steps == "Genome assembly") {
      tagList(
        tabsetPanel(
          tabPanel("Table", DT::dataTableOutput("table")),
        )
      )
    } else if (input$analysis_steps == "Taxonomic classification") {
      tagList(
        tabsetPanel(
          tabPanel("Table", DT::dataTableOutput("table")),
        )
      )
    } else if (input$analysis_steps == "Genome annotation") {
      if (is.numeric(annotation_summary[[input$var]])) {
        tagList(
          tabsetPanel(
            tabPanel("Table", DT::dataTableOutput("table")),
            tabPanel("Histogram - Gene lengths distribution", plotlyOutput("hist"))
          )
        )
      } else if (is.list(annotation_summary[[input$var]])) {
        tagList(
          selectInput("index", "Gene name:", choices = unique(annotation_summary$gene_name_consensus_version_appendix)),
          tabsetPanel(
            tabPanel("Table", tableOutput("list_df"))
          )
        )
      } else {
        tagList(
          tabsetPanel(
            tabPanel("Table", DT::dataTableOutput("table"))
          )
        )
      }
    } else if (input$analysis_steps == "Metabolic functions") {
      if (input$var == "pathway_group") {
        tagList(
          uiOutput("completeness_level_ui"),  # Add the completeness level radio button to the UI
          plotlyOutput("donut")
        )
      } else if (input$var == "module_completeness") {
        plotlyOutput("scatter")
      } else if (input$var == "description") {
        tagList(
          uiOutput("pathway_group_ui"),  # Add the pathway_group dropdown menu to the UI
          DT::dataTableOutput("table")  # Add a table output for the description variable
        )
      }
    }  
    })
  
  output$var_ui <- renderUI({
    if (input$analysis_steps == "Genome assembly") {
      selectInput("var", "Results:", choices = c("Assembly summary" = "assembly_summary_short", 
                                                  "Extended assembly results" = "assembly_summary"
      ))
    } else if (input$analysis_steps == "Taxonomic classification") {
      selectInput("var", "Results:", choices = c("Taxonomy summary" = "taxonomy_summary_short", 
                                                "Extended taxonomy results" = "taxonomy_summary", 
                                                "Other related organisms" = "other_related_organisms"
      ))
    } else if (input$analysis_steps == "Genome annotation") {
      selectInput("var", "Results:", choices = c("Annotation overview" = "length_bp", 
                                                  "Gene name (detailed)" = "gene_name_list", 
                                                  "Gene product (detailed)" = "gene_product_list", 
                                                  "K number (KEGG Orthology)" = "K_list", 
                                                  "EC number (Enzyme Commission)" = "EC_list",
                                                  "COG number (Clusters of Orthologous Groups)" = "COG_list",
                                                  "Detailed gene product info" = "gene_product_list"
      ))
    } else if (input$analysis_steps == "Metabolic functions") {
      selectInput("var", "Results:", choices = c("Pathway description" = "description", 
                                                  "Pathway group" = "pathway_group", 
                                                  "Pathway completeness" = "module_completeness"
      ))
    }  
  })
  
  
  # Add a server output for the description dropdown menu
  output$description_ui <- renderUI({
    if (input$analysis_steps == "Metabolic functions" && input$var == "module_completeness") {
      selectInput("selected_description", "Select Description:", choices = unique(metabolic_summary_microbeannotator$description))
    }
  })
  
  output$pathway_group_ui <- renderUI({
    if (input$analysis_steps == "Metabolic functions" && input$var == "description") {
      selectInput("selected_pathway_group", "Select pathway group:", choices = c("All", unique(metabolic_summary_microbeannotator$pathway_group)))
    }
  })
  
  # Add a server output for the completeness level radio button
  output$completeness_level_ui <- renderUI({
    if (input$analysis_steps == "Metabolic functions" && input$var == "pathway_group") {
      radioButtons("selected_completeness_level", "Select pathway completeness level:", choices = c("Pathway completeness > 0", "Pathway completeness > 50", "Pathway completeness > 80", "Pathway completeness = 100"))
    }
  })
  
  output$table <- DT::renderDataTable({
    if (input$analysis_steps == "Genome annotation") {
      df <- annotation_summary |>
        dplyr::rename("Consensus gene name" = gene_name_consensus,
                      "Confidence in consensus gene name" = gene_name_consensus_confidence,
                      "Consensus gene product" = gene_product_consensus,
                      "Confidence in consensus gene product" = gene_product_consensus_confidence,
                      "Gene start position" = start_position,
                      "Gene end position" = end_position,
                      "Gene length (bp)" = length_bp,
                      "Strand orientation" = strand
        )
      if (!is.list(annotation_summary[[input$var]])) {
        df <- df[ , !sapply(df, is.list)]  # Omit list columns when a non-list variable is selected
      }
      df <- df[ , !grepl("nucleotide|locus_tag|amino_acid|gene_name_consensus_version_appendix", names(df))]  # Omit columns that contain "nucleotide", "locus_tag", "amino_acid" or "gene_name_consensus_version_appendix" in their names
      df <- DT::datatable(df, options = list(pageLength = 25)) |>
        DT::formatStyle(columns = c("Consensus gene name",
                                    "Confidence in consensus gene name",
                                    "Consensus gene product",
                                    "Confidence in consensus gene product",
                                    "Gene start position",
                                    "Gene end position",
                                    "Gene length (bp)",
                                    "Strand orientation"),
                        valueColumns = 'Confidence in consensus gene name',
                        backgroundColor = styleEqual(c("complete_consensus", 
                                                       "high_probability_consensus", 
                                                       "only_one_result",
                                                       "medium_probability_consensus", 
                                                       "low_probability_consensus", 
                                                       "no_consensus", 
                                                       NA), 
                                                     c(rgb(0/255, 128/255, 0/255, alpha = 0.1), 
                                                       rgb(173/255, 255/255, 47/255, alpha = 0.1), 
                                                       rgb(173/255, 255/255, 47/255, alpha = 0.1), 
                                                       rgb(255/255, 255/255, 0/255, alpha = 0.1), 
                                                       rgb(255/255, 165/255, 0/255, alpha = 0.1), 
                                                       rgb(255/255, 0/255, 0/255, alpha = 0.1), 
                                                       rgb(255/255, 0/255, 0/255, alpha = 0.1)))) |>
        DT::formatStyle('Gene length (bp)',
                        background = styleColorBar(df$'Gene length (bp)','lightblue'),
                        backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',
                        backgroundPosition = 'center')
      
    } else if (input$analysis_steps == "Metabolic functions" && input$var == "description") {
      df <- metabolic_summary_microbeannotator
      df <- df[ , c("pathway_group", "description", "module_completeness", "M_number", "comment")]
      
      # Subset the data based on the selected pathway_group
      if (input$selected_pathway_group != "All") {
        df <- df[df$pathway_group == input$selected_pathway_group, ]
      }
      
      # Order by pathway_group and module_completeness
      df <- df |>
        dplyr::arrange(pathway_group, desc(module_completeness)) |> 
        dplyr::filter(module_completeness > 0) |>
        dplyr::mutate(module_completeness = as.numeric(module_completeness)) |>
        dplyr::rename("Pathway group" = pathway_group,
                                 "Pathway completeness" = module_completeness,
                                 "M number (KEGG Module)" = M_number,
                                 "Additional info" = comment,
                                 "Pathway description" = description)
      
      DT::datatable(df, options = list(pageLength = 25)) |>
        DT::formatStyle('Pathway completeness',
                        background = styleColorBar(range(1:100),'lightblue'),
                        backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',
                        backgroundPosition = 'center')
    } else if (input$analysis_steps == "Taxonomic classification") {
      if (input$var == "taxonomy_summary") {
        df <- taxonomy_summary |>
          dplyr::rename("GTDB-Tk result" = result,
                        "Data type" = data_type,
                        "Explanation" = explanation)
      } else if (input$var == "taxonomy_summary_short") {
        df <- taxonomy_summary_short |>
          dplyr::rename("GTDB-Tk result" = result,
                        "Data type" = data_type,
                        "Explanation" = explanation)
      } else if (input$var == "other_related_organisms") {
        df <- other_related_organisms |>
          dplyr::rename("Genome ID (RefSeq)" = genome_id,
                        "Species name" = species_name,
                        "Radius" = radius,)
      }
      df <- df[ , !sapply(df, is.list)]  # Omit list columns
      df <- DT::datatable(df, options = list(pageLength = 25))
    } else if (input$analysis_steps == "Genome assembly") {
      if (input$var == "assembly_summary") {
        df <- assembly_summary |>
          dplyr::rename("Data type" = data_type,
                        "LJA assembler" = lja_assembler,
                        "SPAdes assembler" = spades_assembler,
                        "Canu assembler" = canu_assembler,
                        "CISA contig integrator" = cisa_assembler)
        
      } else if (input$var == "assembly_summary_short") {
        df <- assembly_summary_short |>
          dplyr::rename("Data type" = data_type,
                        "LJA assembler" = lja_assembler,
                        "SPAdes assembler" = spades_assembler,
                        "Canu assembler" = canu_assembler,
                        "CISA contig integrator" = cisa_assembler)
      }
      df <- df[ , !sapply(df, is.list)]  # Omit list columns
      df <- DT::datatable(df, options = list(pageLength = 25))
    }
  })

  output$hist <- renderPlotly({
    if (input$analysis_steps == "Genome annotation" && is.numeric(annotation_summary[[input$var]]) && !input$var %in% c("start_position", "end_position")) {
      # Set the desired bin width
      bin_width <- 100  # Adjust this value according to your preference
      
      # Filter out infinite values and calculate the maximum
      max_value <- max(annotation_summary[[input$var]], na.rm = TRUE, finite = TRUE)
      
      # Calculate the breaks based on the bin width
      breaks <- seq(0, max_value + bin_width, by = bin_width)
      
      # Create the histogram with specified bin width
      data <- annotation_summary[[input$var]]
      
      plot_ly(x = ~data, type = "histogram", histnorm = "count", xbins = list(start = 0, end = max_value + bin_width, size = bin_width), 
              marker = list(color = '#005AB5', line = list(color = 'rgba(0, 0, 0, 1)', width = 2))) |>
        layout(title = input$var, xaxis = list(title = input$var))
    }
  })
  
  output$list_df <- renderTable({
    if (input$analysis_steps == "Genome annotation" && is.list(annotation_summary[[input$var]])) {
      df <- annotation_summary[[input$var]][[which(annotation_summary$gene_name_consensus_version_appendix == input$index)]]
      df <- df |>
        rownames_to_column("Row") |>  
        dplyr::mutate(Row = case_when(
          Row == "gene_name_consensus" ~ "Consensus gene name",
          Row == "gene_name_consensus_confidence" ~ "Confidence in consensus gene name",
          Row == "gene_bakta" ~ "Gene name by Bakta",
          Row == "gene_prokka" ~ "Gene name by Prokka",
          Row == "gene_name_consensus_version_appendix" ~ "Consensus gene name with version number",
          TRUE ~ Row
        )) |>
        dplyr::rename_with(~case_when(
          .x == "Row" ~ "Origin / type / method",
          .x == "gene_name" ~ "Gene name",
          .x == "gene_product" ~ "Gene product",
          .x == "locus_tag_prokka" ~ "Prokka locus tag",
          .x == "locus_tag_bakta" ~ "Bakta locus tag",
          .x == "nucleotide_code" ~ "Nucleotide sequence",
          .x == "K_number" ~ "K number (KEGG Orthology)",
          .x == "EC_number" ~ "EC number (Enzyme Commission)",
          .x == "COG_number" ~ "COG number (Clusters of Orthologous Groups)",
          TRUE ~ .x
        ), any_of(c("Row", "gene_name", "locus_tag_prokka", "locus_tag_bakta", "nucleotide_code")))
      df
    }
  })
  
  output$donut <- renderPlotly({
    if (input$analysis_steps == "Metabolic functions" && input$var == "pathway_group") {
      df <- metabolic_summary_microbeannotator
      df <- df[ , c(input$var, "description", "module_completeness")]
      
      # Subset the data based on the selected completeness level
      if (input$selected_completeness_level == "Pathway completeness > 0") {
        df <- df[df$module_completeness > 0, ]
      } else if (input$selected_completeness_level == "Pathway completeness > 50") {
        df <- df[df$module_completeness > 50, ]
      } else if (input$selected_completeness_level == "Pathway completeness > 80") {
        df <- df[df$module_completeness > 80, ]
      } else if (input$selected_completeness_level == "Pathway completeness = 100") {
        df <- df[df$module_completeness == 100, ]
      }
      
      p <- df |>
        group_by_at(input$var) |>
        summarise(description = length(description)) |>
        plot_ly(labels = ~get(input$var), values = ~description, type = 'pie', hole = 0.4) |>
        layout(title = input$selected_completeness_level, showlegend = T)  # Set the title to the selected completeness level
      p
    }
  })
  
  output$scatter <- renderPlotly({
    if (input$analysis_steps == "Metabolic functions" && input$var == "module_completeness") {
      df <- metabolic_summary_microbeannotator
      df <- df[ , c("description", "module_completeness")]
      
      # Add a new column that checks if the description matches the selected description
      df$selected <- df$description == input$selected_description
      
            # Order the data frame by description and module_completeness
      df <- df |> arrange(desc(description), desc(module_completeness))
      
      
      # Use the new column to conditionally set the color, size, and opacity of the markers
      p <- df |>
        plot_ly(x = ~module_completeness, y = ~description, type = 'scatter', mode = 'markers',
                marker = list(size = ~ifelse(selected, 15, 6), color = ~ifelse(selected, '#DC3220', '#005AB5'), opacity = ~ifelse(selected, 1, 0.6))) |>
        layout(title = "Metabolic module completeness vs pathway description", 
               xaxis = list(title = "Metabolic module completeness"), 
               yaxis = list(title = input$selected_description, titlefont = list(color = '#DC3220'), ticktext = list(), tickvals = list()))  # Set y-axis label to selected description and color it red
      p
    }
  })
  
}

# Run the application 
shinyApp(ui = ui, 
         server = server)
