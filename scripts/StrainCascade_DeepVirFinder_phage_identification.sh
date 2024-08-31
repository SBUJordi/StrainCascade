#!/bin/bash

# StrainCascade_DeepVirFinder_phage_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <input_file> <output_dir> <sample_name> <threads> <genome_annotation_main_abs> <functional_analysis_main_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
input_file=$5
output_dir=$6
sample_name=$7
threads=$8
genome_assembly_main_abs=$9
functional_analysis_main_abs=${10}

source "$utils_file"

# Define paths and variables
matching_files=($(ls "$apptainer_images_dir"/straincascade_phage_detection*.sif 2> /dev/null))

if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

straincascade_phage_detection=${matching_files[0]}

echo "Using $input_file for DeepVirFind detection."

# Create output directory
deepvirfinder_output_dir="$output_dir/DeepVirFinder_phage_identification_results"
create_directory "$deepvirfinder_output_dir"  

# Retrieve the analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (DeepVirFinder phage identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for alingment."
fi

# Log start time
echo "Started DeepVirFinder at $(date)" >> "$logs_dir/DeepVirFinder_phage_identification.log"

# Prepare DeepVirFinder command with Theano flags for better CPU utilization
deepvirfinder_cmd="THEANO_FLAGS='device=cpu,floatX=float32,openmp=True' python /opt/conda/envs/deepvirfinder_env/bin/dvf.py \
             -i /mnt/input/temp_input_reads_DeepVirFinder.fasta \
             -o /mnt/output \
             -c "$threads" \
             -l 200"

# Run DeepVirFinder
apptainer exec \
    --bind "$(dirname "$input_file")":/mnt/input \
    --bind "$deepvirfinder_output_dir":/mnt/output \
    "$straincascade_phage_detection" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate deepvirfinder_env && \

                  cp /mnt/input/$(basename "$input_file") /mnt/input/temp_input_reads_DeepVirFinder.fasta && \

                  $deepvirfinder_cmd && \

                  rm /mnt/input/temp_input_reads_DeepVirFinder.fasta" 2>&1 | tee -a "$logs_dir/DeepVirFinder_phage_identification.log"

# Log end time
echo "Finished DeepVirFinder at $(date)" >> "$logs_dir/DeepVirFinder_phage_identification.log"

# Check output and log file sizes
echo "Output directory size: $(du -sh "$deepvirfinder_output_dir" | cut -f1)" >> "$logs_dir/DeepVirFinder_phage_identification.log"
echo "Log file size: $(du -sh "$logs_dir/DeepVirFinder_phage_identification.log" | cut -f1)" >> "$logs_dir/DeepVirFinder_phage_identification.log"