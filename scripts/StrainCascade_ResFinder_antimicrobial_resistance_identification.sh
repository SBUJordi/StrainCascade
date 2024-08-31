#!/bin/bash

# StrainCascade_ResFinder_antimicrobial_resistance_identification.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 11 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <output_dir> <sample_name> <threads> <sequencing_type> <genome_assembly_main_abs> <functional_analysis_main_abs> <databases_dir>"
    exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
output_dir=$5
sample_name=$6
threads=$7
sequencing_type=$8
genome_assembly_main_abs=$9
functional_analysis_main_abs=${10}
databases_dir=${11}

# Define resfinder_sequencing_type based on the sequencing_type argument
case $sequencing_type in
  pacbio-raw)
    resfinder_sequencing_type=""
    ;;
  pacbio-corr)
    resfinder_sequencing_type=""
    ;;
  pacbio-hifi)
    resfinder_sequencing_type=""
    ;;
  nano-raw)
    resfinder_sequencing_type="--nanopore"
    ;;
  nano-corr)
    resfinder_sequencing_type="--nanopore"
    ;;
  nano-hq)
    resfinder_sequencing_type="--nanopore"
    ;;
  *)
    echo "Error: Invalid sequencing type. Please use one of the StrainCascade sequencing types."
    exit 1
    ;;
esac

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
resfinder_output_dir="$output_dir/ResFinder_AMR_identification_results"
create_directory "$resfinder_output_dir"  

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (ResFinder identification) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for ResFinder."
fi

# Run ResFinder annotation
log "$logs_dir" "ResFinder_annotation.log" "Running ResFinder AMR identification for $analysis_assembly_file in $resfinder_output_dir"

# Prepare ResFinder command
# ResFinder's standard coverage and identity cut-off values were used
resfinder_cmd="python /opt/conda/envs/resfinder_env/lib/python3.8/site-packages/resfinder/run_resfinder.py \
    -ifa /mnt/input/temp_input_assembly_ResFinder.fasta \
    -o /mnt/output \
    -s 'Other' \
    --acquired \
    --point \
    --ignore_missing_species \
    -l 0.6 \
    -t 0.8 \
    -l_p 0.6 \
    -t_p 0.8 \
    -db_res /mnt/resfinder_db \
    -db_point /mnt/pointfinder_db \
    -db_disinf /mnt/disinfinder_db \
    --output_aln \
    -acq \
    -d \
    -u \
    $resfinder_sequencing_type"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$resfinder_output_dir":/mnt/output \
    --bind "$databases_dir/":/mnt \
    "$straincascade_taxonomic_functional_analysis" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate resfinder_env && \
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_ResFinder.fasta && \
                  $resfinder_cmd && \
                  rm /mnt/input/temp_input_assembly_ResFinder.fasta" 2>&1

# Rename files and copy specific files to functional_analysis_main_abs
# Add suffix sample_name_AMR_ResFinder before the last dot in the filenames
for file in "${resfinder_output_dir}"/*; do
    dir=$(dirname "$file")
    base=$(basename "$file")
    extension="${base##*.}"
    filename="${base%.*}"
    new_base="${filename}_${sample_name}_AMR_identification_ResFinder.${extension}"
    mv "$file" "$dir/$new_base"
done

# Find and copy .txt and .fsa files to functional_analysis_main_abs
output_files=$(find "$resfinder_output_dir" -mindepth 1 \( -name "*.txt" -o -name "*.fsa" \) -type f)
if [ -n "$output_files" ]; then
    for file in $output_files; do
        cp "$file" "$functional_analysis_main_abs"
    done
else
    echo "Error: No (suitable) files found in $resfinder_output_dir"
fi