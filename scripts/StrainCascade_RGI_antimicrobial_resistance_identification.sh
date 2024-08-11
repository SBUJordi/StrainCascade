#!/bin/bash

# StrainCascade_RGI_antimicrobial_resistance_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_annotation_main_abs> <functional_analysis_main_abs> <databases_dir>"
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
rgi_output_dir="$output_dir/RGI_AMR_identification_results"
create_directory "$rgi_output_dir"  

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (RGI AMR identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for RGI."
fi

# Run RGI AMR identificatio
log "$logs_dir" "RGI_AMR_identification.log" "Running RGI AMR identificatio for $analysis_assembly_file in $rgi_output_dir"

# Prepare RGI command
output_file="${sample_name}_identification_rgi"
heatmap_file="${sample_name}_heatmap_identification_rgi"

rgi_cmd1="rgi main \
             -i /mnt/input/temp_input_assembly_RGI.fasta \
             -o /mnt/output/"$output_file" \
             -t contig \
             -a BLAST \
              -n $threads \
              --local\
              --clean"

rgi_cmd2="rgi heatmap \
             -i /mnt/output/ \
             -o /mnt/output/"$heatmap_file" \
             --cluster both"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$rgi_output_dir":/mnt/output \
    --bind "$databases_dir/rgi_db":/mnt/rgi_db \
    "$straincascade_taxonomic_functional_analysis" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate rgi_env && \
                  
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_RGI.fasta && \
                  
                  # Load the CARD database 
                  rgi load --card_json /mnt/rgi_db/card.json --local && \
                  
                  $rgi_cmd1 && \
                  $rgi_cmd2 && \
                  
                  rm /mnt/input/temp_input_assembly_RGI.fasta" 2>&1

# Copy specific files to functional_analysis_main_abs
output_files=$(find "$rgi_output_dir" -mindepth 1 \( -name "${sample_name}_*.txt" -o -name "${sample_name}_*.tsv" \) -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$functional_analysis_main_abs"
    done
else
    echo "Error: No (suitable) files found in $rgi_output_dir"
fi