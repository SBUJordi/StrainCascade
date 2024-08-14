#!/bin/bash

# StrainCascade_DeepVirFinder_phage_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_annotation_main_abs> <functional_analysis_main_abs>"
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

# Load utils
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Define paths and variables
matching_files=($(ls "$apptainer_images_dir"/straincascade_phage_detection*.sif 2> /dev/null))

if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

straincascade_phage_detection=${matching_files[0]}

# Create output directory
deepvirfinder_output_dir="$output_dir/DeepVirFinder_phage_identification_results"
create_directory "$deepvirfinder_output_dir"  

# Retrieve the analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (DeepVirFinder phage identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for DeepVirFinder."
fi

# Log start time
echo "Started DeepVirFinder at $(date)" >> "$logs_dir/DeepVirFinder_phage_identification.log"

# Prepare DeepVirFinder command with Theano flags for better CPU utilization
deepvirfinder_cmd="THEANO_FLAGS='device=cpu,floatX=float32,openmp=True' python /opt/conda/envs/deepvirfinder_env/bin/dvf.py \
             -i /mnt/input/temp_input_assembly_DeepVirFinder.fasta \
             -o /mnt/output \
             -c "$threads" \
             -l 200"

# Function to monitor CPU usage
monitor_cpu() {
    while true; do
        cpu_usage=$(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}')
        echo "$(date): CPU Usage: $cpu_usage" >> "$logs_dir/DeepVirFinder_cpu_usage.log"
        sleep 60
    done
}

# Start CPU monitoring in background
monitor_cpu &
monitor_pid=$!

# Run DeepVirFinder
apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$deepvirfinder_output_dir":/mnt/output \
    "$straincascade_phage_detection" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate deepvirfinder_env && \
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_DeepVirFinder.fasta && \
                  $deepvirfinder_cmd && \
                  rm /mnt/input/temp_input_assembly_DeepVirFinder.fasta" 2>&1 | tee -a "$logs_dir/DeepVirFinder_phage_identification.log"

# Stop CPU monitoring
kill $monitor_pid

# Log end time
echo "Finished DeepVirFinder at $(date)" >> "$logs_dir/DeepVirFinder_phage_identification.log"

# Check output and log file sizes
echo "Output directory size: $(du -sh "$deepvirfinder_output_dir" | cut -f1)" >> "$logs_dir/DeepVirFinder_phage_identification.log"
echo "Log file size: $(du -sh "$logs_dir/DeepVirFinder_phage_identification.log" | cut -f1)" >> "$logs_dir/DeepVirFinder_phage_identification.log"