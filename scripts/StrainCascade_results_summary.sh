#!/bin/bash

# BIOMES_results_summary.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

#!/bin/bash

# BIOMES_results_summary.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <output_directory> <sample_name> <miniconda_path> <genome_assembly_main_abs>"
    exit 1
fi
script_dir=$1
logs_dir=$2
output_dir=$3
sample_name=$4
miniconda_path=$5
genome_assembly_main_abs=$6
taxonomic_classification_main_abs=$7
genome_annotation_main_abs=$8

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Create output directory
create_output_directory "$output_dir"  # Ensure the output directory is created

# Define the paths to the TSV files
quast_assembly_evaluation_tsv=$(find "$genome_assembly_main_abs" -type f -name "*quast_assembly_evaluation.tsv")
taxonomy_gtdbtk_tsv=$(find "$taxonomic_classification_main_abs" -type f -name "*summary.tsv")
annotation_prokka_tsv=$(find "$genome_annotation_main_abs" -type f -name "*annotation_prokka.tsv")
annotation_bakta_tsv=$(find "$genome_annotation_main_abs" -type f -name "*annotation_bakta.tsv")
annotation_bakta_hypotheticals_tsv=$(find "$genome_annotation_main_abs" -type f -name "*annotation_bakta.hypotheticals.tsv")

# Activate the default Python environment
log "$logs_dir" "results_summary.log" "Activating Conda environment: BIOMES_WGseq_python"
source "$miniconda_path/etc/profile.d/conda.sh"
conda activate BIOMES_WGseq_python

# Run the BIOMES_TSV_summariser.py script
python "${script_dir}/BIOMES_TSV_summariser.py" "$quast_assembly_evaluation_tsv" "$taxonomy_gtdbtk_tsv" "$annotation_prokka_tsv" "$annotation_bakta_tsv" "$annotation_bakta_hypotheticals_tsv" "$(dirname $genome_assembly_main_abs)" "$sample_name"

# Deactivate the BIOMES_WGseq_python environment
conda deactivate