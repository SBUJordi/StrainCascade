#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_data_integration.sh
# Description: Integrates and processes genome annotation results from multiple tools
#             (Bakta, Prokka, MicrobeAnnotator) and generates comprehensive reports

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> 
          <results_integration_abs> <genome_assembly_main_abs> <version> <sample_name>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 9 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly GENOME_ASSEMBLY_DIR="$6"
readonly RESULTS_DIR="$7"
readonly VERSION="$8"
readonly SAMPLE_NAME="$9"

# Derived constants
readonly QS_FILES_DIR="$RESULTS_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$QS_FILES_DIR"

# Find required Apptainer image
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Log script initialization
log "$LOGS_DIR" "$LOG_NAME" "Starting StrainCascade data integration for sample: $SAMPLE_NAME"

# Retrieve and validate assembly file
readonly ANALYSIS_ASSEMBLY_FILE=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$ANALYSIS_ASSEMBLY_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping data integration."
    exit 0
fi
log "$LOGS_DIR" "$LOG_NAME" "Using assembly file: $ANALYSIS_ASSEMBLY_FILE"

# Find annotation result files
readonly BAKTA_QS=$(find "$QS_FILES_DIR" -name 'bakta_results.qs' -print -quit)
readonly PROKKA_QS=$(find "$QS_FILES_DIR" -name 'prokka_results.qs' -print -quit)
readonly MICROBEANNOTATOR_QS=$(find "$QS_FILES_DIR" -name 'microbeannotator_results.qs' -print -quit)

# Build R script arguments
r_args=""
[[ -n "$BAKTA_QS" ]] && r_args+=" --bakta $(basename "$BAKTA_QS")"
[[ -n "$PROKKA_QS" ]] && r_args+=" --prokka $(basename "$PROKKA_QS")"
[[ -n "$MICROBEANNOTATOR_QS" ]] && r_args+=" --microbeannotator $(basename "$MICROBEANNOTATOR_QS")"

# Log annotation files found
log "$LOGS_DIR" "$LOG_NAME" "Found annotation files:
Bakta: ${BAKTA_QS:-not available}
Prokka: ${PROKKA_QS:-not available}
MicrobeAnnotator: ${MICROBEANNOTATOR_QS:-not available}"

# Process annotation results if available
if [[ -n "$r_args" ]]; then
    # Run genome annotation aggregation
    log "$LOGS_DIR" "$LOG_NAME" "Running genome annotation aggregation"
    
    apptainer exec \
        --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
        --bind "${QS_FILES_DIR}:/mnt/input_output" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_genome_annotation_aggregation.R" \
        --output_dir "/mnt/input_output" \
        $r_args \
        --version "$VERSION"

    # Run genome annotation integration if aggregation successful
    if [[ -f "${QS_FILES_DIR}/annotation_results_aggregated.qs" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Running genome annotation integration"
        
        apptainer exec \
            --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
            --bind "${QS_FILES_DIR}:/mnt/input_output" \
            "$r_sif" \
            Rscript "/mnt/r_script_dir/R_genome_annotation_integration.R" \
            --input_file "/mnt/input_output/annotation_results_aggregated.qs" \
            --output_dir "/mnt/input_output" \
            --version "$VERSION"
    else
        log "$LOGS_DIR" "$LOG_NAME" "Error: Aggregation failed. Skipping integration step."
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "No input files available for aggregation and integration"
fi

# Clean up previous runs
rm -f "${RESULTS_DIR}/*.RData"

# Generate final RData file
log "$LOGS_DIR" "$LOG_NAME" "Generating final RData file"

apptainer exec \
    --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
    --bind "${QS_FILES_DIR}:/mnt/input" \
    --bind "${RESULTS_DIR}:/mnt/output" \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_qs2RData.R" \
    --sample_name "${SAMPLE_NAME}" \
    --input_dir "/mnt/input" \
    --output_dir "/mnt/output" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Final data integration failed"
    exit 0
}

# Locate and validate RData file
readonly RDATA_FILE=$(find "$RESULTS_DIR" -name "StrainCascade_Results_${SAMPLE_NAME}*.RData" -print -quit)
if [[ -z "$RDATA_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: RData file not found for sample ${SAMPLE_NAME}"
    exit 0
fi

# Generate final report
log "$LOGS_DIR" "$LOG_NAME" "Generating StrainCascade report"

apptainer exec \
    --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
    --bind "${RESULTS_DIR}:/mnt/input_output" \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/StrainCascade_report_generation.R" \
    --sample_name "${SAMPLE_NAME}" \
    --RData "/mnt/input_output/$(basename "$RDATA_FILE")" \
    --output_dir "/mnt/input_output" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: StrainCascade report generation failed"
    rm -f "$RESULTS_DIR"/StrainCascade_Analysis_Report_*.Rmd 2>/dev/null || true
    exit 0
}

# Remove any existing StrainCascade Analysis Report Rmd files
rm -f "$RESULTS_DIR"/StrainCascade_Analysis_Report_*.Rmd 2>/dev/null || true

log "$LOGS_DIR" "$LOG_NAME" "StrainCascade data integration completed successfully"