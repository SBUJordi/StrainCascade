#!/bin/bash

# StrainCascade_SPAdes_assembly.sh - Version 1.0.0
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

# Determine the appropriate SPAdes flag based on the sequencing type
case "$sequencing_type" in
    "pacbio-hifi"|"pacbio-corr"|"nano-corr"|"nano-hq")
        spades_sequencing_type="-s"
        ;;
    *)
        log "$logs_dir" "SPAdes_assembly.log" "Skipping SPAdes due to incompatible sequencing type: $sequencing_type. SPAdes 'pacbio-hifi', 'pacbio-corr', 'nano-corr', or 'nano-hq'."
        echo "Skipping due to incompatible sequencing type: $sequencing_type. SPAdes requires 'pacbio-hifi', 'pacbio-corr', 'nano-corr', or 'nano-hq'."
        exit 0  # Exit gracefully, allowing the pipeline to continue
        ;;
esac

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'straincascade_genome_assembly*.sif' -print)

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_genome_assembly=${matching_files[0]}

# Create output directory
spades_output_dir="$output_dir/SPAdes_assembly_results"
create_directory "$spades_output_dir"  # Ensure the output directory is created

## Start with the actual processing ##
# Check if the file informed_genome_size_estimation.txt is in the genome_assembly_main_abs directory
if [ -f "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" ]; then
    estimate_genome_size="no"
else
    estimate_genome_size="yes"
fi

# Run SPAdes with the appropriate option based on the sequencing type
log "$logs_dir" "SPAdes_assembly.log" "Running SPAdes for $input_file as sequencing type $sequencing_type in $spades_output_dir"

apptainer exec \
    --bind "$(dirname "$input_file")":/mnt/input \
    --bind "$spades_output_dir":/mnt/output \
    "$straincascade_genome_assembly" spades.py \
    -t "$threads" \
    --isolate \
    -o /mnt/output \
    $spades_sequencing_type /mnt/input/$(basename "$input_file")
    
# Rename contigs.fasta files in the output directory
cd "$spades_output_dir" || { log "$logs_dir" "SPAdes_assembly.log" "Error changing directory to $spades_output_dir."; exit 1; }
for assembly_file in contigs.fasta; do
  prefix="${sample_name}_assembly_spades"  # Use the provided sample name
  new_name="${prefix}.fasta"
  mv "$assembly_file" "$new_name"
  log "$logs_dir" "SPAdes_assembly.log" "Renamed $assembly_file to $new_name"
  
  # Copy the renamed assembly file to genome_assembly_main_abs
  cp "$new_name" "${genome_assembly_main_abs}/$new_name"
  log "$logs_dir" "SPAdes_assembly.log" "Copied $new_name to $genome_assembly_main_abs"
done

# Estimate genome size if necessary
if [ "$estimate_genome_size" = "yes" ]; then

  # Get a list of all .fasta files
  fasta_files=(${genome_assembly_main_abs}/*.fasta)

  # Select a random .fasta file
  random_fasta_file="${fasta_files[RANDOM % ${#fasta_files[@]}]}"

  # Estimate genome size
  genome_size_estimation=$(grep -v ">" "$random_fasta_file" | tr -d '\n' | tr -cd 'ACGTNacgtn' | wc -c)
  
  if (( genome_size_estimation > 1300000 )); then
      if (( genome_size_estimation > 10000000 )); then
          echo "Warning: Genome size estimation is $genome_size_estimation bp, which is implausibly large for gut bacteria (https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a)."
      else
          echo "Genome size estimation is $genome_size_estimation bp. This is a valid estimation."
          echo "$genome_size_estimation" > "${genome_assembly_main_abs}/informed_genome_size_estimation.txt"
      fi
  elif (( genome_size_estimation <= 1300000 )); then
      echo "Warning: Genome size estimation is $genome_size_estimation bp, which may be too small for a cell to replicate independently in nature and could indicate an issue (https://www.science.org/doi/10.1126/science.1114057, Fig. 1)."
  fi
fi