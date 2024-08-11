#!/bin/bash

# StrainCascade_GTDB-Tk_taxonomy.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_directory> <sample_name> <threads> <genome_assembly_main_abs> <taxonomic_classification_main_abs> <databases_dir>" 
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
threads=$6
genome_assembly_main_abs=$7
taxonomic_classification_main_abs=$8
databases_dir=$9

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

# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'python_3.12.4*.sif' -print)

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_python=${matching_files[0]}

# Create output directory
gtdbtk_output_dir="$output_dir/GTDB-Tk_taxonomy_results"
create_directory "$gtdbtk_output_dir"  # Ensure the output directory is created

# Retrieve the analysis assembly file from genome_assembly_main_abs using the new function
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (GTDB-Tk taxonomy) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for analysis."
fi

# Define the input assembly directory and copy the analysis_assembly_file into it 
gtdbtk_input_dir="$gtdbtk_output_dir/input"
create_directory "$gtdbtk_input_dir"  # Ensure the input directory is created
cp "$analysis_assembly_file" "$gtdbtk_input_dir" # Copy the analysis assembly file to the input directory

## Run GTDB-Tk classify_wf for the analysis assembly file 
log "$logs_dir" "GTDB-Tk_taxonomy.log" "Running GTDB-Tk classification for $analysis_assembly_file in $gtdbtk_output_dir"

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

# Execute the GTDB-Tk classify_wf command inside the Apptainer container
apptainer exec \
    --bind "$gtdbtk_input_dir:/mnt/input" \
    --bind "$gtdbtk_output_dir:/mnt/output" \
    --bind "$latest_db_dir:/mnt/gtdbtk_db" \
    "$straincascade_taxonomic_functional_analysis" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate gtdbtk_env && \
                  export GTDBTK_DATA_PATH='/mnt/gtdbtk_db' && \
                  gtdbtk classify_wf \
                      --genome_dir /mnt/input \
                      --out_dir /mnt/output \
                      --prefix \"${sample_name}_taxonomy_gtdbtk\" \
                      --skip_ani_screen \
                      --extension .fasta \
                      --cpus $threads" 2>&1

# Check exit status of GTDB-Tk classify_wf
if [ $? -ne 0 ]; then
    # Find the first .fasta file in the directory
    analysis_file_failed=$(find "$gtdbtk_input_dir" -type f -name "*.fasta" | head -n 1)
    if [ ! -z "$analysis_file_failed" ]; then
        # Extract just the filename from the path
        analysis_file_failed=$(basename "$analysis_file_failed")
        log "$logs_dir" "GTDB-Tk_taxonomy.log" "Error: GTDB-Tk classify_wf failed for file $analysis_file_failed in directory $gtdbtk_input_dir"
    else
        log "$logs_dir" "GTDB-Tk_taxonomy.log" "Error: GTDB-Tk classify_wf failed for $gtdbtk_input_dir, but no .fasta files were found."
    fi
    exit 0
fi

# Extract the organism name & taxonomic level from the GTDB-Tk classification results
# Run the Python script and pass the output_dir and sample_name as arguments
apptainer exec \
    --bind "${script_dir}:/mnt/script_dir" \
    --bind "${gtdbtk_output_dir}:/mnt/gtdbtk_output_dir" \
    "$straincascade_python" \
    python "/mnt/script_dir/StrainCascade_extract_organism_name_gtdbtk.py" --output_dir "/mnt/gtdbtk_output_dir" --sample_name "${sample_name}"

## Copy the GTDB-Tk classification results to taxonomic_classification_main_abs
# Copy the *summary.tsv files to taxonomic_classification_main_abs
summary_files=$(find "$gtdbtk_output_dir/classify" -type f -name "*taxonomy_gtdbtk*summary.tsv")
if [ -n "$summary_files" ]; then
    for file in $summary_files; do
        # Extract the filename from the full path
        filename=$(basename "$file")
        # Extract the prefix before "taxonomy_gtdbtk"
        prefix="${filename%%taxonomy_gtdbtk*}"
        # Construct the new filename
        new_filename="${prefix}taxonomy_gtdbtk_summary.tsv"
        # Copy the file to the new location with the new name
        cp "$file" "$taxonomic_classification_main_abs/$new_filename"
    done
else
    echo "Error: No summary.tsv files found in $gtdbtk_output_dir/classify"
fi

# Copy the *.tree files individually in for loop to taxonomic_classification_main_abs
tree_files=$(find "$gtdbtk_output_dir/classify" -type f -name "*.tree")
if [ -n "$tree_files" ]; then
    for file in $tree_files; do
        cp "$file" "$taxonomic_classification_main_abs"
    done
else
    echo "Error: No .tree files found in $gtdbtk_output_dir/classify"
fi

# Copy the *_gtdbtk_organism_name.txt files to taxonomic_classification_main_abs
organism_name_files=$(find "$gtdbtk_output_dir" -type f -name "*gtdbtk_organism_name.txt")
if [ -n "$organism_name_files" ]; then
    for file in $organism_name_files; do
        cp "$file" "$taxonomic_classification_main_abs"
    done
else
    echo "Error: No gtdbtk_organism_name.txt files found in $gtdbtk_output_dir"
fi

# Copy the *level_organism_name.txt files to taxonomic_classification_main_abs
level_organism_name_files=$(find "$gtdbtk_output_dir" -type f -name "*level_organism_name.txt")
if [ -n "$level_organism_name_files" ]; then
    for file in $level_organism_name_files; do
        cp "$file" "$taxonomic_classification_main_abs"
    done
else
    echo "Error: No level_organism_name.txt files found in $gtdbtk_output_dir"
fi

## Clean up the GTDB-Tk input directory
rm -rf "$gtdbtk_input_dir"