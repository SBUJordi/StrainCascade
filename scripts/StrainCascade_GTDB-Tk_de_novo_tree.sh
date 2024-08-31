#!/bin/bash

# StrainCascade_GTDB-Tk_de_novo_tree.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <output_directory> <sample_name> <threads> <genome_assembly_main_abs> <taxonomic_classification_main_abs> <databases_dir>" 
    exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
output_dir=$5
sample_name=$6
threads=$7
genome_assembly_main_abs=$8
taxonomic_classification_main_abs=$9
databases_dir=${10}

source "$utils_file"

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
matching_files=($(ls "$apptainer_images_dir"/straincascade_taxonomic_functional_analysis*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_taxonomic_functional_analysis=${matching_files[0]}

# Create output directory
de_novo_tree_output_dir="$output_dir/GTDB-Tk_de_novo_tree_results"
create_directory "$de_novo_tree_output_dir"

# Retrieve the analysis assembly file from genome_assembly_main_abs using the new function
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (GTDB-Tk de novo tree) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for analysis."
fi

# Define the input assembly directory and copy the analysis_assembly_file into it 
de_novo_tree_input_dir="$de_novo_tree_output_dir/input"
create_directory "$de_novo_tree_input_dir"  # Ensure the input directory is created
cp "$analysis_assembly_file" "$de_novo_tree_input_dir" # Copy the analysis assembly file to the input directory

# Read the phylum description from the *_gtdbtk_phylum_name.txt file
phylum_file=$(find "$output_dir/GTDB-Tk_taxonomy_results" -type f -name "*gtdbtk_phylum_name.txt")
if [ -f "$phylum_file" ]; then
    outgroup=$(cat "$phylum_file")
    if [ -z "$outgroup" ]; then
        echo "Warning: Phylum file is empty. Skipping this module (GTDB-Tk de novo tree) and continuing with the next script in the pipeline."
        exit 0
    else
        echo "Phylum file found. Continuing with this module (GTDB-Tk de novo tree)."
    fi
else
    echo "Warning: Phylum file not found. Skipping this module (GTDB-Tk de novo tree) and continuing with the next script in the pipeline."
    exit 0
fi

# Run GTDB-Tk de_novo_wf with the new input assembly directory
log "$logs_dir" "GTDB-Tk_de_novo_tree.log" "Running GTDB-Tk de_novo_wf for $analysis_assembly_file in $de_novo_tree_output_dir"

# Add the --extension flag based on the extension of your assembly files (e.g., .fasta)
assembly_extension=".fasta"

# Find directories starting with "release" and sort them to get the latest
latest_db_dir=$(find "$databases_dir/gtdbtk_db" -type d -name 'release*' | sort -V | tail -n 1)

# Check if the directory exists and is not empty
if [[ -z "$latest_db_dir" ]]; then
    echo "No GTDB-Tk database found."
    exit 1
else
    echo "Using GTDB-Tk database: $latest_db_dir"
fi

# Execute the GTDB-Tk de_novo_wf command inside the Apptainer container
apptainer exec \
    --bind "$de_novo_tree_input_dir:/mnt/input" \
    --bind "$de_novo_tree_output_dir:/mnt/output" \
    --bind "$latest_db_dir:/mnt/gtdbtk_db" \
    "$straincascade_taxonomic_functional_analysis" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate gtdbtk_env && \
                  export GTDBTK_DATA_PATH='/mnt/gtdbtk_db' && \
                  gtdbtk de_novo_wf \
                      --genome_dir /mnt/input \
                      --outgroup_taxon \"$outgroup\" \
                      --out_dir /mnt/output \
                      --prefix \"${sample_name}_de_novo_tree_gtdbtk\" \
                      --extension $assembly_extension \
                      --cpus $threads \
                      --bacteria" 2>&1

# Check exit status of GTDB-Tk de_novo_wf
if [ $? -ne 0 ]; then
    log "$logs_dir" "GTDB-Tk_de_novo_tree.log" "Error: GTDB-Tk de_novo_wf failed for $de_novo_tree_input_dir"
    exit 1
fi

## Copy the GTDB-Tk output files to taxonomic_classification_main_abs
cp "$de_novo_tree_output_dir"/*.decorated.tree* "$taxonomic_classification_main_abs"

## Clean up the GTDB-Tk input directory
rm -rf "$de_novo_tree_input_dir"