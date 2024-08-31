#!/bin/bash

# StrainCascade_NGMLR_BBMap_coverage.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 10 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <input_file> <output_dir> <sample_name> <sequencing_type> <threads> <genome_assembly_main_abs>"
    exit 1
fi

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

# Check if sequencing type contains or is the pattern "pacbio" or "nano"; this covers all sequencing types supported by StrainCascade
if [[ "$sequencing_type" == *"pacbio"* ]]; then
    ngmlr_sequencing_type="pacbio"
elif [[ "$sequencing_type" == *"nano"* ]]; then
    ngmlr_sequencing_type="ont"
else
    log "$logs_dir" "NGMLR_BBMap_coverage.log" "Skipping BBMap due to sequencing type: $sequencing_type. StrainCascade was developped for 'pacbio' or 'nanopore' data."
    echo "Skipping due to incompatible sequencing type: $sequencing_type. StrainCascade was developped for 'pacbio' or 'nanopore' data."
    exit 0  # Exit gracefully, allowing the pipeline to continue
fi

source "$utils_file"

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
ngmlr_bbmap_output_dir="$output_dir/NGMLR_BBMap_coverage_results"
create_directory "$ngmlr_bbmap_output_dir"  # Ensure the output directory is created

# Retrieve the analysis assembly file from genome_assembly_main_abs using the new function
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (BBMap coverage) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file" for analysis.
fi

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file does not exist."
    log "$logs_dir" "NGMLR_BBMap_coverage.log" "Error: Input file $input_file does not exist."
    exit 1
else
    echo "Calculating maximum read length from $input_file"

    # Determine if the file is gzipped
    if [[ "$input_file" == *.gz ]]; then
        # For gzipped files
        max_read_length=$(gunzip -c "$input_file" | awk '{if(NR%4==2) print length}' | sort -nr | head -1)
    else
        # For non-gzipped files
        max_read_length=$(awk '{if(NR%4==2) print length}' "$input_file" | sort -nr | head -1)
    fi

    # Check if max_read_length was calculated correctly
    if [ -z "$max_read_length" ]; then
        echo "Error: Failed to calculate maximum read length."
        log "$logs_dir" "NGMLR_BBMap_coverage.log" "Error: Failed to calculate maximum read length."
        exit 1
    else
        echo "Maximum read length: $max_read_length"

        # Add a buffer to the maximum read length
        max_read_length=$((max_read_length + 100))

        echo "Maximum read length including buffer (+100bp): $max_read_length"
    fi
fi

# Run BBMap or NGMLR based on the maximum read length
if [ "$max_read_length" -le 600 ]; then
    log "$logs_dir" "NGMLR_BBMap_coverage.log" "Running BBMap for $input_file and $analysis_assembly_file"

    apptainer exec \
        --bind "$(dirname "$input_file")":/mnt/sequencing_file_dir \
        --bind "$(dirname "$analysis_assembly_file")":/mnt/assembly_file_dir \
        --bind "$ngmlr_bbmap_output_dir":/mnt/output \
        "$straincascade_assembly_qc_refinement" \
        /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                    conda activate bbmap_ngmlr_env && \
                    bbmap.sh \
                    in=/mnt/sequencing_file_dir/$(basename "$input_file") \
                    ref=/mnt/assembly_file_dir/$(basename "$analysis_assembly_file") \
                    out=/mnt/output/"${sample_name}_mapped.sam" \
                    fastareadlen="$max_read_length" \
                    overwrite=true \
                    path=/mnt/output" 2>&1

else
    log "$logs_dir" "NGMLR_BBMap_coverage.log" "Running NGMLR for $input_file and $analysis_assembly_file"
    apptainer exec \
        --bind "$(dirname "$input_file")":/mnt/sequencing_file_dir \
        --bind "$(dirname "$analysis_assembly_file")":/mnt/assembly_file_dir \
        --bind "$ngmlr_bbmap_output_dir":/mnt/output \
        "$straincascade_assembly_qc_refinement" \
        /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                    conda activate bbmap_ngmlr_env && \
                    ngmlr \
                    -q /mnt/sequencing_file_dir/$(basename "$input_file") \
                    -r /mnt/assembly_file_dir/$(basename "$analysis_assembly_file") \
                    -o /mnt/output/"${sample_name}_mapped.sam" \
                    -x "$ngmlr_sequencing_type" \
                    -t $threads" 2>&1
fi

# Generate coverage statistics
apptainer exec \
    --bind "$ngmlr_bbmap_output_dir":/mnt/input_output \
    "$straincascade_assembly_qc_refinement" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                conda activate bbmap_ngmlr_env && \
                pileup.sh \
                in=/mnt/input_output/$(basename "${sample_name}_mapped.sam") \
                out=/mnt/input_output/$(basename "${analysis_assembly_file%.*}_coverage.txt") \
                overwrite=true" 2>&1
log "$logs_dir" "NGMLR_BBMap_coverage.log" "Generated coverage statistics for $analysis_assembly_file"

# Copy the *_coverage.txt file to genome_assembly_main_abs
coverage_file=$(find "$ngmlr_bbmap_output_dir" -type f -name "*_coverage.txt")
if [ -n "$coverage_file" ]; then
    cp "$coverage_file" "$genome_assembly_main_abs"
else
    echo "Error: No _coverage.txt file found in $ngmlr_bbmap_output_dir"
fi

# Remove all files ending with .ngm from genome_assembly_main_abs
rm -f "$genome_assembly_main_abs"/*.ngm
echo "All .ngm files have been removed from $genome_assembly_main_abs"