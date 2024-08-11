#!/bin/bash

# BIOMES_result_processing.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <output_dir> <sample_name> <miniconda_path> <main_results_dir_abs> <results_integration_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
output_dir=$3
sample_name=$4
miniconda_path=$5
main_results_dir_abs=$6
results_integration_abs=$7

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Create a new directory for R analysis
R_analysis_dir="${results_integration_abs}/R_analysis_${sample_name}"
create_output_directory "$R_analysis_dir"

# Activate the BIOMES_WGseq_R miniconda environment
log "$logs_dir" "result_processing.log" "Activating Conda environment: BIOMES_WGseq_R"
source "$miniconda_path/etc/profile.d/conda.sh"
conda activate BIOMES_WGseq_R

# Run the R script
log "$logs_dir" "result_processing.log" "Running R script: BIOMES_result_processing_for_R.r"
Rscript_output=$(Rscript "${script_dir}/BIOMES_result_processing_for_R.r" "$main_results_dir_abs" "$R_analysis_dir" "$sample_name" 2>&1)
Rscript_exit_code=$?

# Check if the R script exited with an error
if [ $Rscript_exit_code -ne 0 ]; then
  # Log the error message
  log "$logs_dir" "result_processing.log" "Error in R script: $Rscript_output"
else
  # Log a success message
  log "$logs_dir" "result_processing.log" "R script completed successfully"
fi

# Deactivate the BIOMES_WGseq_R environment
log "$logs_dir" "result_processing.log" "Deactivating Conda environment: BIOMES_WGseq_R"
conda deactivate

# Move the BIOMES_WGseq_shiny_app.r script to the R_analysis_dir directory
log "$logs_dir" "result_processing.log" "Moving BIOMES_WGseq_shiny_app.r script to $R_analysis_dir"
mv "${script_dir}/BIOMES_WGseq_shiny_app.r" "$R_analysis_dir"

# Create an Rproject file
log "$logs_dir" "result_processing.log" "Creating Rproject file in $R_analysis_dir"
echo "Version: 1.0" > "${R_analysis_dir}/${sample_name}_BIOMES_WGseq_Rproject.Rproj"
echo "RestoreWorkspace: Default" >> "${R_analysis_dir}/${sample_name}_BIOMES_WGseq_Rproject.Rproj"
echo "SaveWorkspace: Default" >> "${R_analysis_dir}/${sample_name}_BIOMES_WGseq_Rproject.Rproj"
echo "AlwaysSaveHistory: Default" >> "${R_analysis_dir}/${sample_name}_BIOMES_WGseq_Rproject.Rproj"