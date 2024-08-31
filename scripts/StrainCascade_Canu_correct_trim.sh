#!/bin/bash

# StrainCascade_Canu_correct_trim.sh - Version 1.0.1
# Author: Sebastian Bruno Ulrich Jordi

if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <input_file> <output_dir> <sample_name> <sequencing_type> <genome_assembly_main_abs> <threads>"
    exit 1
fi

# Assign the command line arguments to named variables
script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
input_file=$5
output_dir=$6
sample_name=$7
sequencing_type=$8
threads=$9
genome_assembly_main_abs=${10}

source "$utils_file"

# Define canu_sequencing_type and operation based on the sequencing_type argument
case $sequencing_type in
  pacbio-raw)
    canu_sequencing_type="-pacbio"
    operation="-correct -trim"
    ;;
  pacbio-corr)
    canu_sequencing_type="-pacbio -corrected"
    operation="-trim"
    ;;
  pacbio-hifi)
    canu_sequencing_type="-pacbio-hifi"
    operation=""
    ;;
  nano-raw)
    canu_sequencing_type="-nanopore"
    operation="-correct -trim"
    ;;
  nano-corr)
    canu_sequencing_type="-nanopore -corrected"
    operation="-trim"
    ;;
  nano-hq)
    canu_sequencing_type="-nanopore -corrected"
    operation="-trim"
    ;;
  *)
    echo "Error: Invalid sequencing type. Please use one of the StrainCascade sequencing types: pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw | nano-corr | nano-hq"
    exit 1
    ;;
esac

# If no operation is needed (e.g., for PacBio HiFi), copy the input file to the output and exit
if [ -z "$operation" ]; then
    echo "Note: PacBio HiFi data is typically already corrected and trimmed. Skipping correction and trimming."
    exit 0
fi

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'straincascade_genome_assembly*.sif' -print)

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Exiting."
    exit 1
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_genome_assembly=${matching_files[0]}

# Create output directory
canu_output_dir="$output_dir/Canu_correct_trim_results"
create_directory "$canu_output_dir"

# Set genome size (use a default value or read from a file if available)
if [ -f "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" ]; then
    genome_size=$(awk '{printf "%.1fm\n", $1/1000000}' "${genome_assembly_main_abs}/informed_genome_size_estimation.txt")
else
    genome_size="4.5m"  # Default value
fi

## Start with the actual processing ##
log "$logs_dir" "Canu_correct_trim.log" "Running Canu for $operation on $input_file in $canu_output_dir with sequencing type $sequencing_type"

apptainer exec \
    --bind "$(dirname "$input_file")":/mnt/input \
    --bind "$canu_output_dir":/mnt/output \
    "$straincascade_genome_assembly" canu \
    $operation \
    -p "$sample_name" \
    -d /mnt/output \
    genomeSize="$genome_size" \
    $canu_sequencing_type /mnt/input/$(basename "$input_file") \
    useGrid=false \
    corThreads="$threads" \
    corConcurrency=1

# Find the processed reads file
if [ "$operation" = "-correct -trim" ]; then
    processed_file=$(find "$canu_output_dir" -name "${sample_name}*trimmedReads.fasta.gz")
elif [ "$operation" = "-trim" ]; then
    processed_file=$(find "$canu_output_dir" -name "${sample_name}*trimmedReads.fasta.gz")
elif [ "$operation" = "-correct" ]; then
    processed_file=$(find "$canu_output_dir" -name "${sample_name}*correctedReads.fasta.gz")
else
    echo "Unexpected operation: $operation"
    exit 1
fi

if [ -f "$processed_file" ]; then
    # Decompress the file
    gunzip -c "$processed_file" > "${canu_output_dir}/${sample_name}_processed.fasta"
    
    # Overwrite the input file with the processed reads
    mv "${canu_output_dir}/${sample_name}_processed.fasta" "$input_file"
    log "$logs_dir" "Canu_correct_trim.log" "Replaced $input_file with processed reads"
else
    log "$logs_dir" "Canu_correct_trim.log" "Error: Processed reads file not found"
    exit 1
fi