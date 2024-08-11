#!/bin/bash

# BIOMES_IslandPath_genomic_islands_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <output_directory> <sample_name> <miniconda_path> <genome_annotation_main_abs> <functional_analysis_main_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
output_dir=$3
sample_name=$4
miniconda_path=$5
genome_annotation_main_abs=$6
functional_analysis_main_abs=$7

# Load configuration from the script directory
config_file="$script_dir/config.sh"
if [ -f "$config_file" ]; then
  source "$config_file"
else
  echo "Error: config.sh not found in $script_dir"
  exit 1
fi

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Retrieve the chosen .gbff file from genome_annotation_main_abs
chosen_gbff_file=$(find "$genome_annotation_main_abs" -type f -name "*bakta.gbff")

# Activate the IslandPath-DIMOB miniconda environment
log "$logs_dir" "IslandPath_DIMOB_identification.log" "Activating Conda environment: IslandPath-DIMOB"
source "$miniconda_path/etc/profile.d/conda.sh"
conda activate IslandPath-DIMOB

# Run IslandPath-DIMOB annotation
islandpath_output_dir="$output_dir/IslandPath_DIMOB_results"
create_output_directory "$islandpath_output_dir"  # Ensure the main output directory exists
log "$logs_dir" "IslandPath_DIMOB_identification.log" "Running IslandPath-DIMOB annotation for $chosen_gbff_file in $islandpath_output_dir"

# Use the chosen .gbff file for IslandPath-DIMOB annotation
"$islandpath_path" "$chosen_gbff_file" "$islandpath_output_dir/${sample_name}_IslandPath_DIMOB_results.txt"

# Deactivate the conda environment
conda deactivate

# Copy specific files to functional_analysis_main_abs
output_files=$(find "$islandpath_output_dir" -mindepth 1 -name "*.txt" -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$functional_analysis_main_abs"
    done
else
    echo "Error: No (suitable) files found in $islandpath_output_dir"
fi

find "$genome_annotation_main_abs" -name "dimob_*" -exec rm -rf {} \;

# Remove Dimob.log from the parent directory of output_dir (base_dir)
dimob_log_file="$(dirname "$output_dir")/Dimob.log"
if [ -f "$dimob_log_file" ]; then
    rm "$dimob_log_file"
fi