#!/bin/bash

# StrainCascade_dbCAN3_CAZymes_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_assembly_main_abs> <functional_analysis_main_abs> <databases_dir>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
threads=$6
genome_assembly_main_abs=$7
functional_analysis_main_abs=$8
databases_dir=$9

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
matching_files=($(ls "$apptainer_images_dir"/straincascade_taxonomic_functional_analysis*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_taxonomic_functional_analysis=${matching_files[0]}

# Create output directory
dbcan3_output_dir="$output_dir/dbCAN3_CAZymes_identification_results"
create_directory "$dbcan3_output_dir"  

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (dbCAN3 identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for dbCAN3."
fi

# Run dbCAN3 CAZymes identification
log "$logs_dir" "dbCAN3_CAZymes_identification.log" "Running dbCAN3 CAZymes identification for $analysis_assembly_file in $dbcan3_output_dir"

# Prepare dbCAN3 command
# dbCAN3's standard cut-off values were used
dbcan3_cmd="run_dbcan \
             /mnt/input/temp_input_assembly_dbCAN3.fasta \
             prok \
             --out_dir /mnt/output \
             --db_dir /mnt/dbcan3_db \
             --dia_cpu $threads \
             --hmm_cpu $threads \
             --tf_cpu $threads \
             --stp_cpu $threads \
             --tools all \
             --dia_eval 1e-102 \
             --hmm_eval 1e-15 \
             --hmm_cov 0.35 \
             --tf_eval 1e-4 \
             --tf_cov 0.35 \
             --stp_eval 1e-4 \
             --stp_cov 0.3 \
             --cluster 1 \
             --cgc_dis 2 \
             --cgc_sig_genes all"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$dbcan3_output_dir":/mnt/output \
    --bind "$databases_dir/dbcan3_db/db":/mnt/dbcan3_db \
    "$straincascade_taxonomic_functional_analysis" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate dbcan_env && \
                  
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_dbCAN3.fasta && \
                  
                  $dbcan3_cmd && \
                  
                  rm /mnt/input/temp_input_assembly_dbCAN3.fasta" 2>&1

# Copy specific files to functional_analysis_main_abs
mv "${dbcan3_output_dir}/overview.txt" "${dbcan3_output_dir}/overview_${sample_name}_CAZymes_dbcan3.txt"
mv "${dbcan3_output_dir}/cgc_standard.out" "${dbcan3_output_dir}/cgc_standard_${sample_name}_CAZymes_dbcan3.out"

output_files=$(find "$dbcan3_output_dir" -mindepth 1 -name "*CAZymes_dbcan3*" -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$functional_analysis_main_abs"
    done
else
    echo "Error: No (suitable) files found in $dbcan3_output_dir"
fi