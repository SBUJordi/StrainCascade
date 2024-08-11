#!/bin/bash

# StrainCascade_Canu_assembly.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <input_file> <output_dir> <sample_name> <sequencing_type> <genome_assembly_main_abs> <threads>"
    exit 1
fi

# Assign the command line arguments to named variables
script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
input_file=$4
output_dir=$5
sample_name=$6
sequencing_type=$7
threads=$8
genome_assembly_main_abs=$9

# Define canu_sequencing_type based on the sequencing_type argument
case $sequencing_type in
  pacbio-raw)
    canu_sequencing_type="-raw -pacbio"
    ;;
  pacbio-corr)
    canu_sequencing_type="-corrected -pacbio"
    ;;
  pacbio-hifi)
    canu_sequencing_type="-pacbio-hifi"
    ;;
  nano-raw)
    canu_sequencing_type="-raw -nanopore"
    ;;
  nano-corr)
    canu_sequencing_type="-corrected -nanopore"
    ;;
  nano-hq)
    canu_sequencing_type="-corrected -nanopore"
    ;;
  *)
    echo "Error: Invalid sequencing type. Please use one of the StrainCascade sequencing types."
    exit 1
    ;;
esac

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
canu_output_dir="$output_dir/Canu_assembly_results"
create_directory "$canu_output_dir"  # Ensure the output directory is created

# Set counter variable to 0 (to have a starting estimate incase the genome size is not known at this point)
counter="0"
genome_size="4.5m"


## Start with the actual processing ##
# Initialize estimate_genome_size
if [ -f "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" ]; then
    estimate_genome_size="no"
else
    estimate_genome_size="yes"
fi

while [ "$counter" -eq 0 ] || ([ "$counter" -eq 1 ] && [ "$estimate_genome_size" = "yes" ]); do
  # Check if the file informed_genome_size_estimation.txt is in the genome_assembly_main_abs directory and retrieve if so (effectively it will initiate genome size estimation if this is the first genome assembly module)
  if [ -f "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" ]; then
      estimate_genome_size="no"
      # Read the file, divide by 1,000,000, and append 'm'
      genome_size=$(awk '{printf "%.1fm\n", $1/1000000}' "${genome_assembly_main_abs}/informed_genome_size_estimation.txt")
  else
      estimate_genome_size="yes"
  fi

  # Run the Canu Assembler
  log "$logs_dir" "Canu_assembly.log" "Running Canu Assembler for $input_file in $canu_output_dir with sequencing type $sequencing_type"
apptainer exec \
    --bind "$(dirname "$input_file")":/mnt/input \
    --bind "$canu_output_dir":/mnt/output \
    "$straincascade_genome_assembly" canu \
    -p "$sample_name" \
    -d /mnt/output \
    genomeSize="$genome_size" \
    $canu_sequencing_type /mnt/input/$(basename "$input_file") \
    useGrid=false \
    corThreads="$threads" \
    corConcurrency=1
      
  # Rename *.contigs.fasta files in the output directory
  cd "$canu_output_dir" || { log "$logs_dir" "Canu_assembly.log" "Error changing directory to $canu_output_dir."; exit 1; }
  for assembly_file in *.contigs.fasta; do
    
    # Get the prefix from the directory name (remove the "output_" part)
    prefix="${sample_name}_assembly_canu"  # Use the provided sample name
    new_name="${prefix}.fasta"
    mv "$assembly_file" "$new_name"
    log "$logs_dir" "Canu_assembly.log" "Renamed $assembly_file to $new_name"

    # Copy the renamed assembly file to genome_assembly_main_abs
    cp "$new_name" "${genome_assembly_main_abs}/$new_name"
    log "$logs_dir" "Canu_assembly.log" "Copied $new_name to $genome_assembly_main_abs"
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

  # Increment counter
  counter=$((counter + 1))
done