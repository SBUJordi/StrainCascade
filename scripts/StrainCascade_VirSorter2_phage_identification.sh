#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_VirSorter2_phage_identification.sh
# Description: Performs VirSorter2 phage identification as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <functional_analysis_main_abs>
          <results_integration_abs> <databases_dir> <version>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 13 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly OUTPUT_DIR="$6"
readonly SAMPLE_NAME="$7"
readonly THREADS="$8"
readonly GENOME_ASSEMBLY_DIR="$9"
readonly FUNCTIONAL_ANALYSIS_DIR="${10}"
readonly RESULTS_INTEGRATION_DIR="${11}"
readonly DATABASES_DIR="${12}"
readonly VERSION="${13}"

# Derived constants
readonly VIRSORTER2_OUTPUT_DIR="$OUTPUT_DIR/VirSorter2_phage_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly TEMP_INPUT_FILE="temp_input_assembly_VirSorter2.fasta"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
mkdir -p "$VIRSORTER2_OUTPUT_DIR" "$QS_FILES_DIR"

# Find Apptainer images
straincascade_crisprcas_phage_is_elements_sif=$(find "$APPTAINER_DIR" -name 'straincascade_crisprcas_phage_is_elements*.sif' -print -quit)
r_sif=$(find "$APPTAINER_DIR" -name 'r_4.4.1*.sif' -print -quit)

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Run VirSorter2 analysis
echo "[$(date '+%Y-%m-%d %H:%M')] Running VirSorter2 phage identification for $analysis_assembly_file" >> "$LOGS_DIR/$LOG_NAME"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$VIRSORTER2_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/virsorter2_db":/mnt/virsorter2_db \
    --env VIRSORTER_CONFIG_PATH=/mnt/virsorter2_db/config \
    "$straincascade_crisprcas_phage_is_elements_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate virsorter2_env && \
             cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/$TEMP_INPUT_FILE && \
             virsorter run -i /mnt/input/$TEMP_INPUT_FILE -w /mnt/output -d /mnt/virsorter2_db -j $THREADS --min-length 200 all && \
             rm /mnt/input/$TEMP_INPUT_FILE" || {
    echo "[$(date '+%Y-%m-%d %H:%M')] Error: VirSorter2 analysis failed for $analysis_assembly_file" >> "$LOGS_DIR/$LOG_NAME"
    exit 1
}

# Copy output files to functional analysis directory
echo "[$(date '+%Y-%m-%d %H:%M')] Processing VirSorter2 output files" >> "$LOGS_DIR/$LOG_NAME"

# Simple file copying without loops
for base_name in "final-viral-combined.fa" "final-viral-score.tsv" "final-viral-boundary.tsv"; do
    if [[ -f "$VIRSORTER2_OUTPUT_DIR/$base_name" ]]; then
        extension="${base_name##*.}"
        filename="${base_name%.*}"
        cp "$VIRSORTER2_OUTPUT_DIR/$base_name" "$FUNCTIONAL_ANALYSIS_DIR/${filename}_${SAMPLE_NAME}_virsorter2_identification.${extension}"
    else
        echo "[$(date '+%Y-%m-%d %H:%M')] Warning: Expected output file not found: $base_name" >> "$LOGS_DIR/$LOG_NAME"
    fi
done

# Process results with R
R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Find input files for R processing
combined_fa=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "final-viral-combined_${SAMPLE_NAME}_virsorter2_identification.fa" -print -quit)
score_tsv=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "final-viral-score_${SAMPLE_NAME}_virsorter2_identification.tsv" -print -quit)
boundary_tsv=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "final-viral-boundary_${SAMPLE_NAME}_virsorter2_identification.tsv" -print -quit)

# Check if required files exist
if [[ ! -f "$combined_fa" ]] || [[ ! -f "$score_tsv" ]] || [[ ! -f "$boundary_tsv" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M')] Error: One or more required VirSorter2 output files missing" >> "$LOGS_DIR/$LOG_NAME"
    exit 1
fi

# Log R processing configuration
log "$LOGS_DIR" "$LOG_NAME" "Processing VirSorter2 results with R"


# Run R script for processing
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$FUNCTIONAL_ANALYSIS_DIR":/mnt/input \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_virsorter2.R" \
    --output_dir "/mnt/output" \
    --combined_fa "/mnt/input/$(basename "$combined_fa")" \
    --score_tsv "/mnt/input/$(basename "$score_tsv")" \
    --boundary_tsv "/mnt/input/$(basename "$boundary_tsv")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "VirSorter2 analysis completed successfully"