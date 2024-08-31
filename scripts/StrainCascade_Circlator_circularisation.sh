#!/bin/bash

# StrainCascade_Circlator_circularisation.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <input_file> <output_dir> <sample_name> <genome_assembly_main_abs>"
    exit 1
fi
script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
input_file=$5
output_dir=$6
sample_name=$7
genome_assembly_main_abs=$8

source "$utils_file"

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'straincascade_assembly_qc_refinement*.sif' -print)

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
circlator_output_dir="$output_dir/Circlator_circularisation_results"
create_directory "$circlator_output_dir"  # Ensure the output directory is created

# Initialize assembly to empty
assembly=""

# Look for a file containing "best_ev2" pattern
best_ev2_file=$(find "$genome_assembly_main_abs" -type f -name "*best_ev2*.fasta")
if [ ! -z "$best_ev2_file" ]; then
    assembly=$best_ev2_file
    echo "Using assembly file with 'best_ev2' pattern: $assembly"
else
    # If no "best_ev2" file, look for a file containing "best_ev1" pattern
    best_ev1_file=$(find "$genome_assembly_main_abs" -type f -name "*best_ev1*.fasta")
    if [ ! -z "$best_ev1_file" ]; then
        assembly=$best_ev1_file
        echo "Using assembly file with 'best_ev1' pattern: $assembly"
    else
        # If no "best_ev1" file, just take the first .fasta file found (not acutally random)
        random_fasta_file=$(find "$genome_assembly_main_abs" -type f -name "*.fasta" | head -n 1)
        if [ ! -z "$random_fasta_file" ]; then
            assembly=$random_fasta_file
            echo "No 'best_ev1' or 'best_ev2' patterned assembly file found. Using a random .fasta file: $assembly"
        else
            echo "No assembly files found. Skipping this module (Circlator) and continuing with the next script in the pipeline."
            exit 0
        fi
    fi
fi

# Before running circlator remove previous circlator output directories and files
[ -d "${circlator_output_dir}/${sample_name}_circlator_output" ] && rm -rf "${circlator_output_dir}/${sample_name}_circlator_output"
find "${genome_assembly_main_abs}" -type f -name "*_circularised.fasta" -exec rm {} +

# Run Circlator
apptainer exec \
    --bind "$(dirname "$input_file")":/mnt/sequencing_file_dir \
    --bind "$(dirname "$assembly")":/mnt/assembly_file_dir \
    --bind "$circlator_output_dir":/mnt/output \
    "$straincascade_assembly_qc_refinement" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate circlator_env && \
                  circlator all \
                  /mnt/assembly_file_dir/$(basename "$assembly") \
                  /mnt/sequencing_file_dir/$(basename "$input_file") \
                  /mnt/output/"${sample_name}_circlator_output"" 2>&1
    
# Extract the basename without the extension
basename_without_extension=$(basename "$assembly")
basename_without_extension="${basename_without_extension%.*}"

# Move and rename the 06.fixstart.fasta file
mv "${circlator_output_dir}/${sample_name}_circlator_output/06.fixstart.fasta" "${circlator_output_dir}/${basename_without_extension}_circularised.fasta"

# Copy the renamed fasta file to genome_assembly_main_abs
cp "${circlator_output_dir}/${basename_without_extension}_circularised.fasta" "${genome_assembly_main_abs}/${basename_without_extension}_circularised.fasta"