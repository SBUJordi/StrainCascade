#!/bin/bash

# BIOMES_CISA_merging.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <base_dir> <script_dir> <logs_dir> <output_directory> <sample_name> <miniconda_path> <genome_assembly_main_abs>"
    exit 1
fi

base_dir=$1
script_dir=$2
logs_dir=$3
output_dir=$4
sample_name=$5
miniconda_path=$6
genome_assembly_main_abs=$7

# Load configuration from the script directory
config_file="$script_dir/config.sh"
if [ -f "$config_file" ]; then
  source "$config_file"
else
  echo "Error: config.sh not found in $script_dir"
  exit 1
fi

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Check if there are any files in genome_assembly_main_abs with name pattern *assembly*
if [ -z "$(find "$genome_assembly_main_abs" -type f -name "*assembly*")" ]; then
  log "$logs_dir" "CISA_reconciliation.log" "No assembly to process with CISA. Exiting."
  exit 0
fi

# Read genome size from "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" if it exists
# If it doesn't exist, set genome_size to 4500000
if [ -f "${genome_assembly_main_abs}/informed_genome_size_estimation.txt" ]; then
  genome_size=$(cat "${genome_assembly_main_abs}/informed_genome_size_estimation.txt")
else
  genome_size=4500000
fi

# Create output directory
create_output_directory "$output_dir"  # Ensure the output directory is created

# Create an array to store the paths to the copied assembly files
copied_assemblies=()

# Find all assembly files in genome_assembly_main_abs and copy them to output_dir
for original_assembly in $(find "$genome_assembly_main_abs" -type f -name "*assembly*")
do
  # Copy the assembly file to the output directory
  file=$(basename "$original_assembly")
  assembly="${output_dir}/${file}"
  cp "$original_assembly" "$assembly"
  assembly_abs=$(realpath "$assembly")  # Convert to absolute path
  copied_assemblies+=("$assembly_abs")  # Store the absolute path
done

# Merge assemblies using CISA if there are more than one
if [ ${#copied_assemblies[@]} -gt 1 ]; then
  log "$logs_dir" "CISA_reconciliation.log" "Merging assemblies using CISA"
  merged_assembly="${output_dir}/${sample_name}_merged_assembly_cisa.fasta"
  merged_assembly_abs=$(realpath "$merged_assembly")  # Convert to absolute path  
  
  # Activate the CISA miniconda environment
  log "$logs_dir" "CISA_reconciliation.log" "Activating Conda environment: cisa"
  source "$miniconda_path/etc/profile.d/conda.sh"
  conda activate cisa

  # Convert output_dir to an absolute path
  output_dir_abs=$(realpath "$output_dir")

  # Change the current working directory to output_dir_abs
  pushd "$output_dir_abs"  

  # Prepare the input file for Merge.py
  merge_config="Merge.config"
  echo "count=${#copied_assemblies[@]}" > "$merge_config"
  for i in "${!copied_assemblies[@]}"; do
    echo "data=${copied_assemblies[$i]},title=Contig_m$(($i+1))" >> "$merge_config"
  done
  echo "Master_file=${merged_assembly_abs}" >> "$merge_config"  # Use merged_assembly_abs
  echo "min_length=100" >> "$merge_config"
  echo "Gap=11" >> "$merge_config"

  # Run Merge.py
  python "$cisa_merge_py_path" "$merge_config"

  # Prepare the input file for CISA.py
  cisa_config="CISA.config"
  echo "genome=${genome_size}" > "$cisa_config"
  echo "infile=${merged_assembly_abs}" >> "$cisa_config"  # Use merged_assembly_abs
  echo "outfile=${merged_assembly_abs}" >> "$cisa_config"  # Use merged_assembly_abs
  echo "nucmer=${miniconda_path}/envs/cisa/MUMmer3.23/nucmer" >> "$cisa_config"
  echo "makeblastdb=${miniconda_path}/envs/cisa/bin/makeblastdb" >> "$cisa_config"
  echo "blastn=${miniconda_path}/envs/cisa/bin/blastn" >> "$cisa_config"
  echo "CISA=${miniconda_path}/envs/cisa/CISA1.3" >> "$cisa_config"
  echo "R2_Gap=0.95" >> "$cisa_config"

  # Run CISA.py
  python "$cisa_cisa_py_path" "$cisa_config"

  # Restore the original working directory
  popd  

  # Deactivate the CISA environment
  conda deactivate
fi

# Copy the merged assembly to genome_assembly_main_abs
cp "$merged_assembly_abs" "${genome_assembly_main_abs}/${sample_name}_merged_assembly_cisa.fasta"

# Delete all .p.fa files from the output directory
find "$output_dir" -type f -name "*.p.fa" -delete

# Delete all files in the output directory that match any of the assembly patterns
for pattern in "${assembly_patterns[@]}"
do
  find "$output_dir" -type f -name "*${pattern}" -delete
done