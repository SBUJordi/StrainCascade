#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_CRISPRCasFinder_identification.sh
# Description: Performs CRISPR-Cas identification analysis as part of StrainCascade

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <functional_analysis_main_abs>
          <results_integration_abs> <version>
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
readonly FUNCTIONAL_ANALYSIS_DIR="${10}"
readonly RESULTS_INTEGRATION_DIR="${11}"
readonly VERSION="${12}"

# Source utility functions
source "$UTILS_FILE"

# Create output directory
readonly CRISPRCAS_OUTPUT_DIR="$OUTPUT_DIR/CRISPRCasFinder_identification_results"
create_directory "$CRISPRCAS_OUTPUT_DIR"

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$analysis_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found. Skipping CRISPRCasFinder identification."
    exit 0
fi

log "$LOGS_DIR" "$LOG_NAME" "Using assembly file: $analysis_assembly_file"

# Find the necessary .sif file
straincascade_crisprcas_phage_is_elements=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_crisprcas_phage_is_elements*.sif')

# Run CRISPRCasFinder identification
log "$LOGS_DIR" "$LOG_NAME" "Running CRISPRCasFinder identification for $analysis_assembly_file in $CRISPRCAS_OUTPUT_DIR"

# Setup temporary files and directories
readonly INPUT_DIR=$(dirname "$analysis_assembly_file")
readonly INPUT_FILENAME=$(basename "$analysis_assembly_file")
readonly TEMP_INPUT_FILE="temp_input_assembly_CRISPRCasFinder.fasta"
readonly TEMP_OUTPUT_DIR="temp_crisprcasfinder_output"

# Prepare CRISPRCasFinder command - removed -outdir parameter
readonly CRISPRCAS_CMD="perl /opt/CRISPRCasFinder/CRISPRCasFinder.pl \
    -in /mnt/input/$TEMP_INPUT_FILE \
    -so /opt/CRISPRCasFinder/sel392v2.so \
    -cas \
    -keep \
    -cpuM $THREADS"

# Create temporary directory to catch intermediate files
readonly ORIGINAL_DIR=$(pwd)
readonly TEMP_DIR="$(get_absolute_path "$(pwd)/crisprcasfinder_temp${SAMPLE_NAME}")"
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

# Execute CRISPRCasFinder in the container
apptainer exec \
    --bind "$INPUT_DIR":/mnt/input \
    --bind "$CRISPRCAS_OUTPUT_DIR":/mnt/output \
    --bind "$TEMP_DIR":/workdir \
    "$straincascade_crisprcas_phage_is_elements" \
    bash -c "cd /workdir && \
             cp /mnt/input/$INPUT_FILENAME /mnt/input/$TEMP_INPUT_FILE && \
             source /opt/conda/etc/profile.d/conda.sh && \
             conda activate crisprcasfinder && \
             $CRISPRCAS_CMD && \
             # Move all output files from working directory to output directory
             mv Result_* /mnt/output/ && \
             rm /mnt/input/$TEMP_INPUT_FILE" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: CRISPRCasFinder identification failed"
    cd "$ORIGINAL_DIR"
    rm -rf "$TEMP_DIR"
    exit 1
}

# Copy only the most relevant output files
log "$LOGS_DIR" "$LOG_NAME" "Copying output files"

# Copy main result files
if [[ -f "$CRISPRCAS_OUTPUT_DIR/result.json" ]]; then
    cp "$CRISPRCAS_OUTPUT_DIR/result.json" "$FUNCTIONAL_ANALYSIS_DIR/result_${SAMPLE_NAME}_crisprcasfinder_identification.json"
    log "$LOGS_DIR" "$LOG_NAME" "Copied result.json"
fi

if [[ -f "$CRISPRCAS_OUTPUT_DIR/TSV/CRISPR-Cas_summary.tsv" ]]; then
    cp "$CRISPRCAS_OUTPUT_DIR/TSV/CRISPR-Cas_summary.tsv" "$FUNCTIONAL_ANALYSIS_DIR/CRISPR-Cas_summary_${SAMPLE_NAME}_crisprcasfinder_identification.tsv"
    log "$LOGS_DIR" "$LOG_NAME" "Copied CRISPR-Cas_summary.tsv"
fi

if [[ -f "$CRISPRCAS_OUTPUT_DIR/TSV/Cas_REPORT.tsv" ]]; then
    cp "$CRISPRCAS_OUTPUT_DIR/TSV/Cas_REPORT.tsv" "$FUNCTIONAL_ANALYSIS_DIR/Cas_REPORT_${SAMPLE_NAME}_crisprcasfinder_identification.tsv"
    log "$LOGS_DIR" "$LOG_NAME" "Copied Cas_REPORT.tsv"
fi

# Copy GFF files (if they exist)
if [[ -d "$CRISPRCAS_OUTPUT_DIR/GFF" ]]; then
    find "$CRISPRCAS_OUTPUT_DIR/GFF" -type f -name "*.gff" -exec bash -c '
        for file; do
            base=$(basename "${file%.*}")
            cp "$file" "'"$FUNCTIONAL_ANALYSIS_DIR"'/${base}_'"$SAMPLE_NAME"'_crisprcasfinder_identification.gff"
        done' bash {} +
    log "$LOGS_DIR" "$LOG_NAME" "Copied GFF files"
fi

# Cleanup
cd "$ORIGINAL_DIR"
rm -rf "$TEMP_DIR"

log "$LOGS_DIR" "$LOG_NAME" "CRISPRCasFinder analysis completed successfully"