#!/bin/bash

# StrainCascade_arrow_medaka_polishing.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <input_file> <bam_file> <output_dir> <sample_name> <sequencing_type> <threads> <genome_assembly_main_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
input_file=$4
bam_file=$5
output_dir=$6
sample_name=$7
sequencing_type=$8
threads=$9
genome_assembly_main_abs=${10}

# Determine polishing tool based on sequencing type
if [[ "$sequencing_type" == *"pacbio"* ]]; then
    if [[ -n "$bam_file" && "$bam_file" != "not_available" ]]; then
        polishing_tool="arrow"
    else
        echo "BAM file is required for polishing of PacBio data by arrow but is not available. Continuing with the next module in the pipeline."
        exit 0  # Exit gracefully, allowing the pipeline to continue
    fi
elif [[ "$sequencing_type" == *"nano"* ]]; then
    polishing_tool="medaka"
else
    echo "No valid sequencing type identified. Continuing with the next module in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
fi

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
polishing_output_dir="$output_dir/${polishing_tool}_polishing_results"
create_directory "$polishing_output_dir"  # Ensure the output directory is created

# Retrieve the analysis assembly file from genome_assembly_main_abs using the new function
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (${polishing_tool} polishing) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file assembly for analysis."
fi

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file does not exist."
    log "$logs_dir" "${polishing_tool}_polishing.log" "Error: Input file $input_file does not exist."
    exit 1
else
    echo "Using $input_file as sequencing input file for analysis"
fi

# Set up polishing rounds
polishing_rounds=2
current_assembly="$analysis_assembly_file"

# Perform polishing rounds
for ((round=1; round<=polishing_rounds; round++)); do
    echo "Starting polishing round $round"
    log "$logs_dir" "${polishing_tool}_polishing.log" "Starting polishing round $round"

    round_output_dir="${polishing_output_dir}/round_${round}"
    create_directory "$round_output_dir"

    if [ "$polishing_tool" = "medaka" ]; then
        polishing_cmd="medaka_consensus -i /mnt/input_file/$(basename "$input_file") \
                       -d /mnt/input_assembly/$(basename "$current_assembly") \
                       -o /mnt/output \
                       -t $threads"
        conda_env="medaka_env"
    elif [ "$polishing_tool" = "arrow" ]; then
        polishing_cmd="variantCaller --algorithm=arrow \
                       -j $threads \
                       -r /mnt/input_assembly/$(basename "$current_assembly") \
                       /mnt/input_file/$(basename "$input_file") \
                       -o /mnt/output/consensus.fasta"
        conda_env="genomicconsensus_env"
    else
        echo "Error: Unknown polishing tool: $polishing_tool"
        log "$logs_dir" "polishing.log" "Error: Unknown polishing tool: $polishing_tool"
        exit 1
    fi

    # Run the polishing command
    apptainer exec \
        --bind "$(dirname "$input_file")":/mnt/input_file \
        --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
        --bind "$round_output_dir":/mnt/output \
        "$straincascade_assembly_qc_refinement" \
        /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                      conda activate $conda_env && \
                      $polishing_cmd" 2>&1

    # Check if polishing ran successfully
    if [ $? -eq 0 ]; then
        echo "Round $round of ${polishing_tool} polishing completed successfully."
        log "$logs_dir" "${polishing_tool}_polishing.log" "Round $round of ${polishing_tool} polishing completed successfully."
        
        if [ "$polishing_tool" = "medaka" ] || [ "$polishing_tool" = "arrow" ]; then
            current_assembly="$round_output_dir/consensus.fasta"
        fi
    else
        echo "Error: Round $round of ${polishing_tool} polishing failed."
        log "$logs_dir" "${polishing_tool}_polishing.log" "Error: Round $round of ${polishing_tool} polishing failed."
        exit 1
    fi
done

echo "Final polished assembly: $current_assembly"
log "$logs_dir" "${polishing_tool}_polishing.log" "Final polished assembly: $current_assembly"

# Clean up
rm -f "$polishing_output_dir"/*/temp_input_*