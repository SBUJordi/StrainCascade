#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_Circlator_circularisation.sh
# Description: Performs bacterial genome assembly circularization using Circlator as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    echo "Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <input_file> \
<output_dir> <sample_name> <genome_assembly_main_abs> <results_integration_abs> <version>"
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
readonly GENOME_ASSEMBLY_DIR="$9"
readonly RESULTS_INTEGRATION_DIR="${10}"
readonly VERSION="${11}"

# Derived constants
readonly CIRCLATOR_OUTPUT_DIR="$OUTPUT_DIR/Circlator_circularisation_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly MAX_RUNTIME_SECONDS=720000  # 200 hours -> adjust if you want real time-constraints

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$CIRCLATOR_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Clean up empty FASTA files
log "$LOGS_DIR" "$LOG_NAME" "Checking for empty FASTA files in $GENOME_ASSEMBLY_DIR"
while IFS= read -r -d '' file; do
    if [[ -f "$file" && ! -s "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Removing empty file: $(basename "$file")"
        rm "$file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

# Find required Apptainer images
readonly straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file - prioritizing best_ev2, then best_ev1, then any .fasta
log "$LOGS_DIR" "$LOG_NAME" "Looking for assembly file to process"

assembly=""
if best_ev2_file=$(find "$GENOME_ASSEMBLY_DIR" -type f -name "*best_ev2*.fasta" -print -quit); then
    assembly=$best_ev2_file
    log "$LOGS_DIR" "$LOG_NAME" "Using best_ev2 assembly: $assembly"
elif best_ev1_file=$(find "$GENOME_ASSEMBLY_DIR" -type f -name "*best_ev1*.fasta" -print -quit); then
    assembly=$best_ev1_file
    log "$LOGS_DIR" "$LOG_NAME" "Using best_ev1 assembly: $assembly"
elif any_fasta=$(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print -quit); then
    assembly=$any_fasta
    log "$LOGS_DIR" "$LOG_NAME" "Using available assembly: $assembly"
else
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found. Skipping Circlator."
    exit 0
fi

# Clean previous Circlator output
[[ -d "${CIRCLATOR_OUTPUT_DIR}/${SAMPLE_NAME}_circlator_output" ]] && rm -rf "${CIRCLATOR_OUTPUT_DIR}/${SAMPLE_NAME}_circlator_output"
find "${GENOME_ASSEMBLY_DIR}" -type f -name "*_circularised.fasta" -delete

# Create deterministic entropy source
readonly ENTROPY_FILE="$CIRCLATOR_OUTPUT_DIR/deterministic_entropy_file"
if [[ ! -f "$ENTROPY_FILE" ]]; then
    dd if=/dev/zero bs=1024 count=100 > "$ENTROPY_FILE"
    log "$LOGS_DIR" "$LOG_NAME" "Deterministic entropy file created at $ENTROPY_FILE"
fi

# Run Circlator with time limit
log "$LOGS_DIR" "$LOG_NAME" "Starting Circlator processing"

start_time=$(date +%s)
apptainer exec \
    --bind "$(dirname "$INPUT_FILE")":/mnt/sequencing_file_dir \
    --bind "$(dirname "$assembly")":/mnt/assembly_file_dir \
    --bind "$CIRCLATOR_OUTPUT_DIR":/mnt/output \
    --bind "$ENTROPY_FILE":/dev/random \
    --bind "$ENTROPY_FILE":/dev/urandom \
    "$straincascade_assembly_qc_refinement_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate circlator_env && \
             circlator all \
             /mnt/assembly_file_dir/$(basename "$assembly") \
             /mnt/sequencing_file_dir/$(basename "$INPUT_FILE") \
             /mnt/output/${SAMPLE_NAME}_circlator_output" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Circlator processing failed. Skipping Circlator circularisation."
    exit 0
}

current_time=$(date +%s)
elapsed_time=$((current_time - start_time))
if [ "$elapsed_time" -gt "$MAX_RUNTIME_SECONDS" ]; then
    log "$LOGS_DIR" "$LOG_NAME" "Circlator process running >$MAX_RUNTIME_SECONDS seconds without finishing. Exiting."
    exit 0
fi

# Process output files
basename_no_ext="${assembly%.*}"
basename_no_ext="$(basename "$basename_no_ext")"

output_file="${CIRCLATOR_OUTPUT_DIR}/${basename_no_ext}_circularised.fasta"
mv "${CIRCLATOR_OUTPUT_DIR}/${SAMPLE_NAME}_circlator_output/06.fixstart.fasta" "$output_file"
cp "$output_file" "${GENOME_ASSEMBLY_DIR}/${basename_no_ext}_circularised.fasta"

# Clean up empty FASTA files
log "$LOGS_DIR" "$LOG_NAME" "Checking for empty FASTA files in $GENOME_ASSEMBLY_DIR"
while IFS= read -r -d '' file; do
    if [[ -f "$file" && ! -s "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Removing empty file: $(basename "$file")"
        rm "$file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

log "$LOGS_DIR" "$LOG_NAME" "Circlator processing complete. Processing results with R"

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"
circlator_log=$(find "$CIRCLATOR_OUTPUT_DIR" -name '04.merge.circularise.log' -print -quit)

if [[ -n "$circlator_log" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Processing Circlator results with R."

    # Create a temporary directory within the output directory
    TMP_DIR="${OUTPUT_DIR}/R_temp"
    mkdir -p "$TMP_DIR"

    # Process results with R
    apptainer exec \
        --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
        --bind "$(dirname "$circlator_log")":/mnt/input \
        --bind "${QS_FILES_DIR}:/mnt/output" \
        --bind "${TMP_DIR}:/tmp" \
        "$r_sif" \
        bash -c "export TMPDIR=/tmp && \
                 export R_SESSION_TMPDIR=/tmp && \
                 Rscript /mnt/r_script_dir/R_process_circlator.R \
                 --output_dir '/mnt/output' \
                 --log_file '/mnt/input/$(basename "$circlator_log")' \
                 --version '$VERSION'" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed. Skipping this step."
        exit 0
    }

    # Clean up temporary directory
    rm -rf "$TMP_DIR"
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: Circlator log file not found"
fi