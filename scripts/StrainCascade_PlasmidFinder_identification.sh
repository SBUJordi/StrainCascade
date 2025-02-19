#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_PlasmidFinder_identification.sh
# Description: Performs plasmid identification using PlasmidFinder and processes results as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <genome_assembly_main_abs> <results_integration_abs> 
          <databases_dir> <version>
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
readonly OUTPUT_DIR="$6"
readonly SAMPLE_NAME="$7"
readonly GENOME_ASSEMBLY_DIR="$8"
readonly RESULTS_INTEGRATION_DIR="$9"
readonly DATABASES_DIR="${10}"
readonly VERSION="${11}"

# Derived constants
readonly PLASMIDFINDER_OUTPUT_DIR="$OUTPUT_DIR/PlasmidFinder_plasmid_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"
readonly TEMP_INPUT_FILE="temp_input_assembly_PlasmidFinder.fasta"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$PLASMIDFINDER_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Find container images
straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
[[ -z "$analysis_assembly_file" ]] && {
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found. Skipping PlasmidFinder identification."
    exit 0
}

# Build PlasmidFinder command
readonly PLASMIDFINDER_CMD='python /opt/conda/envs/plasmidfinder_env/bin/plasmidfinder.py \
    -i /mnt/input/'"$TEMP_INPUT_FILE"' \
    -o /mnt/output \
    -mp /opt/conda/envs/plasmidfinder_env/bin/blastn \
    -p /mnt/plasmidfinder_db'

# Run PlasmidFinder
log "$LOGS_DIR" "$LOG_NAME" "Running PlasmidFinder for $analysis_assembly_file"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$PLASMIDFINDER_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/plasmidfinder_db":/mnt/plasmidfinder_db \
    "$straincascade_assembly_qc_refinement_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate plasmidfinder_env && \
             cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/$TEMP_INPUT_FILE && \
             cd /mnt/plasmidfinder_db && \
             python3 INSTALL.py kma_index && \
             $PLASMIDFINDER_CMD && \
             rm /mnt/input/$TEMP_INPUT_FILE" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: PlasmidFinder failed for $analysis_assembly_file"
    exit 1
}

# Process output files
log "$LOGS_DIR" "$LOG_NAME" "Processing PlasmidFinder output files"

# Rename output files with sample prefix
find "$PLASMIDFINDER_OUTPUT_DIR" -type f -exec bash -c '
    for file do
        mv "$file" "$(dirname "$file")/'${SAMPLE_NAME}'_plasmidfinder_$(basename "$file")"
    done
' bash {} +

# Copy relevant files to assembly directory
if output_files=$(find "$PLASMIDFINDER_OUTPUT_DIR" -type f \( -name "*.txt" -o -name "*.json" \)); then
    for file in $output_files; do
        cp "$file" "$GENOME_ASSEMBLY_DIR"
    done
    log "$LOGS_DIR" "$LOG_NAME" "Output files copied to assembly directory"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: No output files found"
    exit 1
fi

# Process results with R
log "$LOGS_DIR" "$LOG_NAME" "Processing results with R"

# Find json output file
json_file=$(find "$GENOME_ASSEMBLY_DIR" -name '*plasmidfinder_data.json' -print -quit)
[[ -z "$json_file" ]] && {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Required JSON file not found"
    exit 1
}

# Log R processing configuration
cat << EOF
Running data processing of PlasmidFinder results with R with the following configuration:
.fasta file: $analysis_assembly_file
.json file: $json_file
Output directory: $QS_FILES_DIR
StrainCascade version v$VERSION
EOF

# Run R script
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$GENOME_ASSEMBLY_DIR":/mnt/input_plasmidfinder \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input_fasta \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_plasmidfinder.R" \
    --output_dir "/mnt/output" \
    --json "/mnt/input_plasmidfinder/$(basename "$json_file")" \
    --fasta "/mnt/input_fasta/$(basename "$analysis_assembly_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "PlasmidFinder analysis completed successfully"