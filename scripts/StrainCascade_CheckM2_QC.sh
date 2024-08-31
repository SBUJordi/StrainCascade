#!/bin/bash

# StrainCascade_CheckM2_QC.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_assembly_main_abs> <databases_dir>"
    exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
output_dir=$5
sample_name=$6
threads=$7
genome_assembly_main_abs=$8
databases_dir=$9

source "$utils_file"

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
matching_files=($(ls "$apptainer_images_dir"/straincascade_assembly_qc_refinement*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_assembly_qc_refinement=${matching_files[0]}

# Create output directory
checkm2_output_dir="$output_dir/CheckM2_coverage_results"
create_directory "$checkm2_output_dir"  # Ensure the output directory is created

# Retrieve the analysis assembly file from genome_assembly_main_abs using the new function
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (CheckM2 quality control) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for analysis."
fi

## Run CheckM2 predict for the analysis assembly file
log "$logs_dir" "CheckM2_predict.log" "Running CheckM2 predict for $analysis_assembly_file in $checkm2_output_dir"

apptainer exec \
    --bind "$analysis_assembly_file":/mnt/assembly_file \
    --bind "$checkm2_output_dir":/mnt/output \
    --bind "$databases_dir":/mnt/databases_dir \
    "$straincascade_assembly_qc_refinement" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate checkm2_env && \
                 checkm2 predict \
                    --threads $threads \
                    --input /mnt/assembly_file \
                    --output-directory /mnt/output \
                    --database_path /mnt/databases_dir/checkm2_db/uniref100.KO.1.dmnd" 2>&1

## Copy the the relevant output to the genome_assembly_main_abs directory
# Extract the basename without extension from analysis_assembly_file
analysis_assembly_basename=$(basename "$analysis_assembly_file" .fasta) # Assuming the extension is .fasta; adjust as necessary

# Find the quality report file
quality_report=$(find "$checkm2_output_dir" -type f -name "*quality_report.tsv")

if [ -n "$quality_report" ]; then
    # Construct the new filename with the prefix and an underscore
    new_filename="${analysis_assembly_basename}_checkm2_quality_report.tsv"
    # Copy the file to genome_assembly_main_abs with the new name
    cp "$quality_report" "${genome_assembly_main_abs}/${new_filename}"
else
    echo "Error: No CheckM2 quality report found in $checkm2_output_dir"
fi