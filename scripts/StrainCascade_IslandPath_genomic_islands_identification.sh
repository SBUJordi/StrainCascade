#!/bin/bash

# StrainCascade_IslandPath_genomic_islands_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <genome_annotation_main_abs> <functional_analysis_main_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
genome_annotation_main_abs=$6
functional_analysis_main_abs=$7

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

# Create output directory
islandpath_output_dir="$output_dir/IslandPath_genomic_island_identification_results"
create_directory "$islandpath_output_dir"  

# Retrieve the chosen .gbff file from genome_annotation_main_abs
bakta_genbank_file=$(find "$genome_annotation_main_abs" -type f -name "*annotation_bakta.gbff" | head -n 1)

if [ -z "$bakta_genbank_file" ]; then
    # If no Bakta .gbff file is found, check for Prokka .gbk file
    prokka_genbank_file=$(find "$genome_annotation_main_abs" -type f -name "*annotation_prokka.gbk" | head -n 1)
    
    if [ -z "$prokka_genbank_file" ]; then
        # If no Prokka .gbk file is found, check for any .gbk, .gbff, or .embl file
        genbank_file=$(find "$genome_annotation_main_abs" -type f \( -name "*.gbk" -o -name "*.gbff" -o -name "*.embl" \) | head -n 1)

        if [ -z "$genbank_file" ]; then
            echo "Error: No suitable input file (with .gbk, .gbff, or .embl extension) found. Skipping this module (IslandPath annotation) and continuing with the next script in the pipeline. (Consider running the Bakta (or Prokka) annotation module before running IslandPath)"
            exit 0
        else
            echo "No presumptive Bakta or Prokka .gbk or .gbff file found. Trying to use $genbank_file for analysis. (Consider running the Bakta annotation module before running IslandPath)"
            genbank_file_to_use="$genbank_file"
            file_origin="non_StrainCascade_tool"
        fi
    else
        echo "Using $prokka_genbank_file (presumptive Prokka generated .gbk file) for analysis."
        genbank_file_to_use="$prokka_genbank_file"
        file_origin="prokka"
    fi
else
    echo "Using $bakta_genbank_file (presumptive Bakta generated .gbff file) for analysis."
    genbank_file_to_use="$bakta_genbank_file"
    file_origin="bakta"
fi

# Run IslandPath annotation
log "$logs_dir" "IslandPath_genomic_island_identification.log" "Running IslandPath annotation for $genbank_file_to_use (presumptively generated with $file_origin) in $islandpath_output_dir"



# Prepare IslandPath command 
islandpath_cmd="/opt/conda/envs/islandpath_env/opt/islandpath/Dimob.pl \
                                      /mnt/input/${sample_name}_${file_origin}_annotation_islandpath.gbk \
                                      /mnt/output/${sample_name}_results_genomic_island_identification_islandpath.gff3"

# Run Prokka using Apptainer
apptainer exec \
    --bind "$(dirname "$genbank_file_to_use")":/mnt/input \
    --bind "$islandpath_output_dir":/mnt/output \
    "$straincascade_taxonomic_functional_analysis" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate islandpath_env && \
                  
                  # Copy and rename the input file (too long names may cause problems; this way the name of the input files is standardised)
                  cp /mnt/input/$(basename "$genbank_file_to_use") /mnt/input/${sample_name}_${file_origin}_annotation_islandpath.gbk && \

                  # Run IslandPath command
                  $islandpath_cmd && \
                  
                  # Remove the temporary file
                  rm /mnt/input/${sample_name}_${file_origin}_annotation_islandpath.gbk" 2>&1

# Copy specific files to functional_analysis_main_abs
output_files=$(find "$islandpath_output_dir" -mindepth 1 -name "*.gff3" -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$functional_analysis_main_abs"
    done
else
    echo "Error: No (suitable) files found in $islandpath_output_dir"
fi

find "$genome_annotation_main_abs" -name "dimob_*" -exec rm -rf {} \;

# Remove Dimob.log from the parent directory of output_dir (base_dir)
dimob_log_file="$(dirname "$output_dir")/Dimob.log"
if [ -f "$dimob_log_file" ]; then
    rm "$dimob_log_file"
fi