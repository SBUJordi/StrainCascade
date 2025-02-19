#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_CheckM2_QC.sh
# Description: Performs CheckM2 quality control analysis as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <results_integration_abs>
          <databases_dir> <version>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 12 ]] || show_usage

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
readonly RESULTS_INTEGRATION_DIR="${10}"
readonly DATABASES_DIR="${11}"
readonly VERSION="${12}"

# Derived constants
readonly CHECKM2_OUTPUT_DIR="$OUTPUT_DIR/CheckM2_QC_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly CHECKM2_DB_PATH="$DATABASES_DIR/checkm2_db/uniref100.KO.1.dmnd"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$CHECKM2_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Working variables
straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$analysis_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found. Skipping CheckM2 quality control."
    exit 0
fi

log "$LOGS_DIR" "$LOG_NAME" "Using assembly file: $analysis_assembly_file"

# Run CheckM2 predict
log "$LOGS_DIR" "$LOG_NAME" "Running CheckM2 predict analysis"

apptainer exec \
    --bind "$analysis_assembly_file":/mnt/assembly_file \
    --bind "$CHECKM2_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR":/mnt/databases_dir \
    "$straincascade_assembly_qc_refinement_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate checkm2_env && \
             checkm2 predict \
                --threads $THREADS \
                --input /mnt/assembly_file \
                --output-directory /mnt/output \
                --database_path /mnt/databases_dir/checkm2_db/uniref100.KO.1.dmnd" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: CheckM2 predict analysis failed"
    exit 1
}

# Process CheckM2 output
analysis_assembly_basename=$(basename "${analysis_assembly_file%.*}")
quality_report=$(find "$CHECKM2_OUTPUT_DIR" -type f -name "*quality_report.tsv")

if [[ -n "$quality_report" ]]; then
    new_filename="${analysis_assembly_basename}_checkm2_quality_report.tsv"
    cp "$quality_report" "$GENOME_ASSEMBLY_DIR/$new_filename"
    log "$LOGS_DIR" "$LOG_NAME" "Copied quality report to: $GENOME_ASSEMBLY_DIR/$new_filename"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: No CheckM2 quality report found in $CHECKM2_OUTPUT_DIR"
    exit 1
fi

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"
tsv_file=$(find "$GENOME_ASSEMBLY_DIR" -name "*checkm2_quality_report.tsv" -print -quit)

# Validate required file exists
if [[ ! -f "$tsv_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Required TSV file not found"
    exit 1
fi

# Log R processing configuration
log "$LOGS_DIR" "$LOG_NAME" "Processing CheckM2 results with R:
Assembly file: $analysis_assembly_file
TSV file: $tsv_file"

# Run R script for processing
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$GENOME_ASSEMBLY_DIR":/mnt/input \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_checkm2.R" \
    --output_dir "/mnt/output" \
    --tsv "/mnt/input/$(basename "$tsv_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "CheckM2 quality control completed successfully"