#!/bin/bash

# StrainCascade_VirSorter2_phage_identification.sh - Version 1.0.0
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
matching_files=($(ls "$apptainer_images_dir"/straincascade_phage_detection*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_phage_detection=${matching_files[0]}

# Create output directory
virsorter2_output_dir="$output_dir/VirSorter2_phage_identification_results"
create_directory "$virsorter2_output_dir"  

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (VirSorter2 phage identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for VirSorter2."
fi

# Run VirSorter2 phage identification
log "$logs_dir" "VirSorter2_phage_identification.log" "Running VirSorter2 phage identificatio for $analysis_assembly_file in $virsorter2_output_dir"

# Prepare VirSorter2 command
virsorter2_cmd="virsorter run \
             -i /mnt/input/temp_input_assembly_VirSorter2.fasta \
             -w /mnt/output \
             -d /mnt/virsorter2_db \
             -j $threads \
             --min-length 200 \
             all"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$virsorter2_output_dir":/mnt/output \
    --bind "$databases_dir/virsorter2_db":/mnt/virsorter2_db \
    "$straincascade_phage_detection" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate virsorter2_env && \
                  
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_VirSorter2.fasta && \
                  
                  $virsorter2_cmd && \
                  
                  rm /mnt/input/temp_input_assembly_VirSorter2.fasta" 2>&1

# Find and rename specific files, then move them to functional_analysis_main_abs
output_files=$(find "$virsorter2_output_dir" -mindepth 1 \( -name "*combined.fa" -o -name "*score.tsv" \) -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        # Extract the directory, base name, and extension
        dir=$(dirname "$file")
        base=$(basename "$file")
        extension="${base##*.}"
        filename="${base%.*}"
        
        # Construct the new file name
        new_filename="${filename}_${sample_name}_virsorter2_identification.${extension}"
        
        # Move and rename the file
        cp "$file" "$functional_analysis_main_abs/$new_filename"
    done
else
    echo "Error: No suitable files found in $virsorter2_output_dir"
fi