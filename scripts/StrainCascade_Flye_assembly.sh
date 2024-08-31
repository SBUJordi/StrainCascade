#!/bin/bash

# StrainCascade_Flye_assembly.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

## External inputs ##
# Check for the correct number of command line arguments
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <input_file> <output_dir> <sample_name> <genome_assembly_main_abs>" "<threads>" "<sequencing_type>"
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

## Define paths and variables for this script ##
# Find the necessary .sif file
straincascade_genome_assembly=$(find_apptainer_sif_file "$apptainer_images_dir" 'straincascade_genome_assembly*.sif')

# Create output directory
flye_output_dir="$output_dir/Flye_assembly_results"
create_directory "$flye_output_dir"  # Ensure the output directory is created

# Set counter variable to 0 (to have a starting estimate incase the genome size is not known at this point)
counter="0"
genome_size="4.5m"

## Start with the actual processing ##
# Initialize estimate_genome_size
if [ -f "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" ]; then
    estimate_genome_size="no"
    # Read the file, divide by 1,000,000, and append 'm'
    genome_size=$(awk '{printf "%.1fm\n", $1/1000000}' "${genome_assembly_main_abs}/informed_genome_size_estimation.txt")
else
    estimate_genome_size="yes"
fi

while [ "$counter" -eq 0 ] || ([ "$counter" -eq 1 ] && [ "$estimate_genome_size" = "yes" ]); do
  # Run the Flye Assembler
  log "$logs_dir" "Flye_assembly.log" "Running Flye Assembler for $input_file in $flye_output_dir"
  apptainer exec \
    --bind "$(dirname "$input_file")":/mnt/input \
    --bind "$flye_output_dir":/mnt/output \
    "$straincascade_genome_assembly" flye \
    --"$sequencing_type" /mnt/input/$(basename "$input_file") \
    --out-dir /mnt/output \
    --genome-size "$genome_size" \
    --threads "$threads" \
    --iterations 2

  # Rename and copy the assembly.fasta file in the output directory
  cd "$flye_output_dir" || { log "$logs_dir" "Flye_assembly.log" "Error changing directory to $flye_output_dir."; exit 1; }

  # Define the paths
  assembly_file="${flye_output_dir}/assembly.fasta"
  new_name="${flye_output_dir}/${sample_name}_assembly_flye.fasta"

  # Rename the assembly.fasta file
  mv "$assembly_file" "$new_name"
  log "$logs_dir" "Flye_assembly.log" "Renamed $assembly_file to $new_name"

  # Copy the renamed assembly file to genome_assembly_main_abs
  cp "$new_name" "${genome_assembly_main_abs}/$(basename "$new_name")"
  log "$logs_dir" "Flye_assembly.log" "Copied $new_name to ${genome_assembly_main_abs}/$(basename "$new_name")"

  # Estimate genome size if necessary
  if ([ "$counter" -eq 0 ] && [ "$estimate_genome_size" = "yes" ]); then
    # Use the newly generated assembly file
    genome_size_estimation=$(grep -v ">" "$new_name" | tr -d '\n' | tr -cd 'ACGTNacgtn' | wc -c)
    
    if (( genome_size_estimation > 1300000 )); then
      if (( genome_size_estimation > 10000000 )); then
        echo "Warning: Genome size estimation is $genome_size_estimation bp, which is implausibly large for gut bacteria (https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a)."
      else
        echo "Genome size estimation is $genome_size_estimation bp. This is a valid estimation."
        echo "$genome_size_estimation" > "${genome_assembly_main_abs}/informed_genome_size_estimation.txt"
        genome_size="$(awk '{printf "%.1fm\n", $1/1000000}' <<< "$genome_size_estimation")"
      fi
    elif (( genome_size_estimation <= 1300000 )); then
      echo "Warning: Genome size estimation is $genome_size_estimation bp, which may be too small for a cell to replicate independently in nature and could indicate an issue (https://www.science.org/doi/10.1126/science.1114057, Fig. 1)."
    fi
  fi

  # Increment counter
  counter=$((counter + 1))
done

# Log completion of the script
log "$logs_dir" "Flye_assembly.log" "Flye assembly completed for $sample_name"