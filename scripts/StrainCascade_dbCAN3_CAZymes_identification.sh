#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_dbCAN3_CAZymes_identification.sh
# Description: Identifies CAZymes using dbCAN3 and processes results for StrainCascade pipeline

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
readonly DBCAN3_OUTPUT_DIR="$OUTPUT_DIR/dbCAN3_CAZymes_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$DBCAN3_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Find required Apptainer images
readonly straincascade_taxonomic_functional_analysis_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_taxonomic_functional_analysis*.sif')
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Initialize logging
log "$LOGS_DIR" "$LOG_NAME" "Starting dbCAN3 CAZymes identification for sample: $SAMPLE_NAME"

# Retrieve and validate assembly file
readonly ANALYSIS_ASSEMBLY_FILE=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$ANALYSIS_ASSEMBLY_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping dbCAN3 identification."
    exit 0
fi
log "$LOGS_DIR" "$LOG_NAME" "Using assembly file: $ANALYSIS_ASSEMBLY_FILE"

# Define dbCAN3 parameters
readonly DBCAN3_PARAMS=(
    "--dia_eval 1e-102"
    "--hmm_eval 1e-15"
    "--hmm_cov 0.35"
    "--tf_eval 1e-4"
    "--tf_cov 0.35"
    "--stp_eval 1e-4"
    "--stp_cov 0.3"
    "--cluster 1"
    "--cgc_dis 2"
    "--cgc_sig_genes all"
)

# Construct dbCAN3 command
readonly DBCAN3_CMD="run_dbcan \
    /tmp/temp_input_assembly_dbCAN3.fasta \
    prok \
    --out_dir /mnt/output \
    --db_dir /mnt/dbcan3_db \
    --dia_cpu $THREADS \
    --hmm_cpu $THREADS \
    --tf_cpu $THREADS \
    --stp_cpu $THREADS \
    --tools all \
    ${DBCAN3_PARAMS[*]}"

# Run dbCAN3 analysis
log "$LOGS_DIR" "$LOG_NAME" "Running dbCAN3 analysis"

apptainer exec \
    --bind "$(dirname "$ANALYSIS_ASSEMBLY_FILE")":/mnt/input \
    --bind "$DBCAN3_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/dbcan3_db/db":/mnt/dbcan3_db \
    "$straincascade_taxonomic_functional_analysis_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate dbcan_env && \
             export TMPDIR=/tmp && \
             mkdir -p \$TMPDIR && \
             cp /mnt/input/$(basename "$ANALYSIS_ASSEMBLY_FILE") /tmp/temp_input_assembly_dbCAN3.fasta && \
             $DBCAN3_CMD" 2>&1 | tee "$LOG_NAME"

# Check for analysis errors
if grep -q "FileNotFoundError: \[Errno 2\] No such file or directory:" "$LOG_NAME"; then
    log "$LOGS_DIR" "$LOG_NAME" "WARNING: No sequences passed filtering criteria. Consider adjusting parameters: ${DBCAN3_PARAMS[*]}"
    exit 0
fi

# Process output files
log "$LOGS_DIR" "$LOG_NAME" "Processing output files"

# Rename output files
mv "$DBCAN3_OUTPUT_DIR/overview.txt" "$DBCAN3_OUTPUT_DIR/overview_${SAMPLE_NAME}_CAZymes_dbcan3.txt"
mv "$DBCAN3_OUTPUT_DIR/cgc_standard.out" "$DBCAN3_OUTPUT_DIR/cgc_standard_${SAMPLE_NAME}_CAZymes_dbcan3.out"

# Copy results to functional analysis directory
readonly OUTPUT_FILES=$(find "$DBCAN3_OUTPUT_DIR" -mindepth 1 -name "*CAZymes_dbcan3*" -type f)
if [[ -n "$OUTPUT_FILES" ]]; then
    for file in $OUTPUT_FILES; do
        cp "$file" "$FUNCTIONAL_ANALYSIS_DIR"
    done
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: No output files found"
    exit 1
fi

# Locate required files for R processing
readonly OUT_FILE=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name '*dbcan3.out' -print -quit)
readonly TXT_FILE=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name '*_dbcan3.txt' -print -quit)

# Validate R input files
if [[ -z "$OUT_FILE" || -z "$TXT_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Missing required files for R processing"
    exit 1
fi

# Process results with R
log "$LOGS_DIR" "$LOG_NAME" "Processing dbCAN3 results with R"

apptainer exec \
    --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
    --bind "${FUNCTIONAL_ANALYSIS_DIR}:/mnt/input" \
    --bind "$(dirname "$ANALYSIS_ASSEMBLY_FILE")":/mnt/input_fasta \
    --bind "${QS_FILES_DIR}:/mnt/output" \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_dbcan3.R" \
    --output_dir "/mnt/output" \
    --out_file "/mnt/input/$(basename "$OUT_FILE")" \
    --txt_file "/mnt/input/$(basename "$TXT_FILE")" \
    --fasta "/mnt/input_fasta/$(basename "$ANALYSIS_ASSEMBLY_FILE")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "dbCAN3 analysis completed successfully"