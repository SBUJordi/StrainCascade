#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_DeepVirFinder_phage_identification.sh
# Description: Identifies phage sequences using DeepVirFinder as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <input_file> 
          <output_dir> <sample_name> <threads> <genome_annotation_main_abs> 
          <functional_analysis_main_abs>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 11 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly INPUT_FILE="$6"
readonly OUTPUT_DIR="$7"
readonly SAMPLE_NAME="$8"
readonly THREADS="$9"
readonly GENOME_ASSEMBLY_DIR="${10}"
readonly FUNCTIONAL_ANALYSIS_DIR="${11}"

# Derived constants
readonly DVF_OUTPUT_DIR="$OUTPUT_DIR/DeepVirFinder_phage_identification_results"
readonly TEMP_INPUT="temp_input_reads_DeepVirFinder.fasta"
readonly MIN_CONTIG_LENGTH=200

# Create compilation directory structure with unique path per job
readonly THEANO_BASE_DIR="$DVF_OUTPUT_DIR/theano_compile_${SAMPLE_NAME}_$$"
mkdir -p "$THEANO_BASE_DIR"

# Source utility functions
source "$UTILS_FILE"

# Create output directory
create_directory "$DVF_OUTPUT_DIR"

# Find appropriate Apptainer image
straincascade_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_crisprcas_phage_is_elements*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
[[ -z "$analysis_assembly_file" ]] && {
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping DeepVirFinder module."
    rm -rf "$THEANO_BASE_DIR"  # Cleanup
    exit 0
}

# Log start of processing
log "$LOGS_DIR" "$LOG_NAME" "Starting DeepVirFinder phage identification"
log "$LOGS_DIR" "$LOG_NAME" "Assembly file: $analysis_assembly_file"
log "$LOGS_DIR" "$LOG_NAME" "Minimum contig length: $MIN_CONTIG_LENGTH"
log "$LOGS_DIR" "$LOG_NAME" "Using Theano compile directory: $THEANO_BASE_DIR"

# Run DeepVirFinder with modified environment
log "$LOGS_DIR" "$LOG_NAME" "Running DeepVirFinder annotation"
apptainer exec \
    --bind "$(dirname "$INPUT_FILE")":/mnt/input \
    --bind "$DVF_OUTPUT_DIR":/mnt/output \
    --bind "$THEANO_BASE_DIR":/theano_compile \
    "$straincascade_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate deepvirfinder_env && \
             cp /mnt/input/$(basename "$INPUT_FILE") /mnt/input/$TEMP_INPUT && \
             THEANO_FLAGS='device=cpu,floatX=float32,openmp=True,base_compiledir=/theano_compile' \
             python /opt/conda/envs/deepvirfinder_env/bin/dvf.py \
             -i /mnt/input/$TEMP_INPUT \
             -o /mnt/output \
             -c $THREADS \
             -l $MIN_CONTIG_LENGTH && \
             rm /mnt/input/$TEMP_INPUT" 2> >(while read -r line; do log "$LOGS_DIR" "$LOG_NAME" "$line"; done)

# Process output files
log "$LOGS_DIR" "$LOG_NAME" "Processing DeepVirFinder output files"
if output_files=$(find "$DVF_OUTPUT_DIR" -mindepth 1 -name "*dvfpred.txt" -type f); then
    while IFS= read -r file; do
        [[ -f "$file" ]] || continue
        
        # Extract filename components and create new name
        base=$(basename "$file")
        extension="${base##*.}"
        filename="${base%.*}"
        new_filename="${filename}_${SAMPLE_NAME}_deepvirFinder_identification.${extension}"
        
        # Copy to functional analysis directory
        cp "$file" "$FUNCTIONAL_ANALYSIS_DIR/$new_filename"
        log "$LOGS_DIR" "$LOG_NAME" "Processed output file: $new_filename"
    done <<< "$output_files"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: No prediction files found in $DVF_OUTPUT_DIR"
fi

# Cleanup Theano compile directory
rm -rf "$THEANO_BASE_DIR"

# Log completion statistics
log "$LOGS_DIR" "$LOG_NAME" "DeepVirFinder phage analysis completed successfully"