#!/bin/bash

# StrainCascade_MicrobeAnnotator_annotation.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_annotation_main_abs> <functional_analysis_main_abs> <databases_dir>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
threads=$6
genome_annotation_main_abs=$7
functional_analysis_main_abs=$8
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
matching_files=($(ls "$apptainer_images_dir"/straincascade_genome_annotation*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_genome_annotation=${matching_files[0]}

# Create output directory
microbeannotator_output_dir="$output_dir/MicrobeAnnotator_annotation_results"
create_directory "$microbeannotator_output_dir"  

# Retrieve the chosen .faa file from genome_annotation_main_abs
bakta_faa_file=$(find "$genome_annotation_main_abs" -type f -name "*annotation_bakta.faa" | head -n 1)

if [ -z "$bakta_faa_file" ]; then
    # If no Bakta .faa file is found, check for Prokka .faa file
    prokka_faa_file=$(find "$genome_annotation_main_abs" -type f -name "*annotation_prokka.faa" | head -n 1)
    
    if [ -z "$prokka_faa_file" ]; then
        # If no Prokka .faa file is found, check for any .faa file
        faa_file=$(find "$genome_annotation_main_abs" -type f -name "*.faa" | head -n 1)
        
        if [ -z "$faa_file" ]; then
            echo "Error: No suitable input file (with .faa extension) found. Skipping this module (MicrobeAnnotator annotation) and continuing with the next script in the pipeline. (Consider running the Bakta (or Prokka) annotation module before running MicrobeAnnotator)"
            exit 0
        else
            echo "No presumptive Bakta or Prokka .faa file found. Trying to use $faa_file for analysis. (Consider running the Bakta annotation module before running MicrobeAnnotator)"
            faa_file_to_use="$faa_file"
            file_origin="non_StainCascade_tool"
        fi
    else
        echo "Using $prokka_faa_file (presumptive Prokka generated .faa file) for analysis."
        faa_file_to_use="$prokka_faa_file"
        file_origin="prokka"
    fi
else
    echo "Using $bakta_faa_file (presumptive Bakta generated .faa file) for analysis."
    faa_file_to_use="$bakta_faa_file"
    file_origin="bakta"
fi

# Run MicrobeAnnotator annotation
log "$logs_dir" "MicrobeAnnotator_annotation.log" "Running MicrobeAnnotator annotation for $faa_file_to_use (presumtively generated with $file_origin) in $microbeannotator_output_dir"





# Prepare MicrobeAnnotator command 
microbeannotator_cmd="microbeannotator -i /mnt/input/${sample_name}_${file_origin}_annotation_microbeannotator.faa \
                                       -d /mnt/microbeannotator_db \
                                       -o /mnt/output \
                                       -m blast \
                                       -p 1 \
                                       -t $threads \
                                       --refine"

# Run Prokka using Apptainer
apptainer exec \
    --bind "$(dirname "$faa_file_to_use")":/mnt/input \
    --bind "$microbeannotator_output_dir":/mnt/output \
    --bind "$databases_dir/microbeannotator_db":/mnt/microbeannotator_db \
    "$straincascade_genome_annotation" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate microbeannotator_env && \
                  
                  # Copy and rename the input file (too long names may cause problems; this way the name of the input files is standardised)
                  cp /mnt/input/$(basename "$faa_file_to_use") /mnt/input/${sample_name}_${file_origin}_annotation_microbeannotator.faa && \
                  
                  # Set TMPDIR environment variable
                  export TMPDIR=/tmp && \

                  # Run MicrobeAnnotator command
                  $microbeannotator_cmd && \
                  
                  # Remove the temporary file
                  rm /mnt/input/${sample_name}_${file_origin}_annotation_microbeannotator.faa" 2>&1

# Check exit status of MicrobeAnnotator
if [ $? -ne 0 ]; then
    log "$logs_dir" "MicrobeAnnotator_annotation.log" "Error: MicrobeAnnotator annotation failed for $faa_file_to_use"
    exit 1
fi

# Rename directories
find "$microbeannotator_output_dir" -mindepth 1 -type d | while read -r dir; do
    mv "$dir" "${dir}_microbeannotator"
done

# Replace all double underscores with single underscores in filenames
find "$microbeannotator_output_dir" -name '*__*' | while read -r file; do
    new_file=$(echo "$file" | sed 's/__/_/g')
    mv "$file" "$new_file"
done

# Copy specific files to genome_annotation_main_abs and functional_analysis_main_abs
output_files=$(find "$microbeannotator_output_dir" -mindepth 1 \( -name "*kofam*" -o -name "*.annot" -o -name "*.ko" \) -type f)
metabolic_files=$(find "$microbeannotator_output_dir" -mindepth 1 -name "metabolic_summary*" -type f)

if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$genome_annotation_main_abs"
    done
else
    echo "Error: No (suitable) output files found in $microbeannotator_output_dir"
fi

if [ -n "$metabolic_files" ]; then
    for file in $metabolic_files; do
        cp "$file" "$functional_analysis_main_abs"
    done
else
    echo "Error: No (suitable) metabolic_summary files found in $microbeannotator_output_dir"
fi