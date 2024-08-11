#!/bin/bash

# StrainCascade_PlasmidFinder_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <genome_assembly_main_abs> <databases_dir>" 
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
genome_assembly_main_abs=$6
databases_dir=$7

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
plasmidfinder_output_dir="$output_dir/PlasmidFinder_plasmid_identification_results"
create_directory "$plasmidfinder_output_dir"  

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (PlasmidFinder identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for PlasmidFinder."
fi

# Run PlasmidFinder plasmid identification
log "$logs_dir" "PlasmidFinder_plasmid_identification.log" "Running PlasmidFinder plasmid identification for $analysis_assembly_file in $plasmidfinder_output_dir"



# Prepare Bakta command
plasmidfinder_cmd='python /opt/conda/envs/plasmidfinder_env/bin/plasmidfinder.py \
                          -i /mnt/input/temp_input_assembly_PlasmidFinder.fasta \
                          -o /mnt/output \
                          -mp /opt/conda/envs/plasmidfinder_env/bin/blastn \
                          -p /mnt/plasmidfinder_db' 

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$plasmidfinder_output_dir":/mnt/output \
    --bind "$databases_dir/plasmidfinder_db":/mnt/plasmidfinder_db \
    "$straincascade_assembly_qc_refinement" \
    /bin/bash -c "
                  source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate plasmidfinder_env && \
                  
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_PlasmidFinder.fasta && \
                  
                  cd /mnt/plasmidfinder_db && \
                  python3 INSTALL.py kma_index && \
                  
                  $plasmidfinder_cmd && \
                  
                  rm /mnt/input/temp_input_assembly_PlasmidFinder.fasta" 2>&1

# Check exit status of Bakta
if [ $? -ne 0 ]; then
    log "$logs_dir" "PlasmidFinder_plasmid_identification.log" "Error: PlasmidFinder identification failed for $analysis_assembly_file"
    exit 1
fi

# Rename and copy specific files to genome_assembly_main_abs
for file in "$plasmidfinder_output_dir"/*; do
    base=$(basename "$file")
    mv "$file" "$plasmidfinder_output_dir/${sample_name}_plasmidfinder_$base" # Add a prefix to the file name
done

output_files=$(find "$plasmidfinder_output_dir" -mindepth 1 \( -name "*.txt" -o -name "*.json" \) -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$genome_assembly_main_abs"
    done
else
    echo "Error: No (suitable) files found in $plasmidfinder_output_dir"
fi