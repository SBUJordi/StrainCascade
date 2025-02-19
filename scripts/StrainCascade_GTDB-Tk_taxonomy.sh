#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_GTDB-Tk_taxonomy.sh
# Description: Performs GTDB-Tk taxonomic classification as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <taxonomic_classification_main_abs>
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
readonly TAXONOMIC_CLASS_DIR="${10}"
readonly RESULTS_INTEGRATION_DIR="${11}"
readonly DATABASES_DIR="${12}"
readonly VERSION="${13}"

# Derived constants
readonly GTDBTK_OUTPUT_DIR="$OUTPUT_DIR/GTDB-Tk_taxonomy_results"
readonly GTDBTK_INPUT_DIR="$GTDBTK_OUTPUT_DIR/input"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"
readonly CLASSIFY_DIR="$GTDBTK_OUTPUT_DIR/classify"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
for dir in "$GTDBTK_OUTPUT_DIR" "$GTDBTK_INPUT_DIR" "$QS_FILES_DIR"; do
    create_directory "$dir"
done

# Find required containers
straincascade_taxonomic_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_taxonomic_functional_analysis*.sif')
python_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'python_3.12.4*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find and validate analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$analysis_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found in $GENOME_ASSEMBLY_DIR"
    exit 0
fi

# Copy assembly file to input directory
cp "$analysis_assembly_file" "$GTDBTK_INPUT_DIR/"

# Find and validate GTDB-Tk database
latest_db_dir=$(find "$DATABASES_DIR/gtdbtk_db" -type d -name 'release*' | sort -V | tail -n 1)
if [[ -z "$latest_db_dir" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No GTDB-Tk database found in $DATABASES_DIR/gtdbtk_db"
    exit 1
fi

# Build GTDB-Tk command
gtdbtk_base_cmd="gtdbtk classify_wf \
    --genome_dir /mnt/input \
    --out_dir /mnt/output \
    --prefix ${SAMPLE_NAME}_taxonomy_gtdbtk \
    --skip_ani_screen \
    --extension .fasta \
    --cpus $THREADS"

# Log configuration
log "$LOGS_DIR" "$LOG_NAME" "Starting GTDB-Tk classification with configuration:
Input assembly: $analysis_assembly_file
Database: $latest_db_dir"

# Run GTDB-Tk classify_wf
apptainer exec \
    --bind "$GTDBTK_INPUT_DIR:/mnt/input" \
    --bind "$GTDBTK_OUTPUT_DIR:/mnt/output" \
    --bind "$latest_db_dir:/mnt/gtdbtk_db" \
    "$straincascade_taxonomic_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate gtdbtk_env && \
             export GTDBTK_DATA_PATH='/mnt/gtdbtk_db' && \
             $gtdbtk_base_cmd" || {
    analysis_file_failed=$(basename "$(find "$GTDBTK_INPUT_DIR" -type f -name "*.fasta" | head -n 1)")
    log "$LOGS_DIR" "$LOG_NAME" "Error: GTDB-Tk classify_wf failed for ${analysis_file_failed:-'no .fasta files found'}"
    exit 0
}

# Extract organism name and taxonomic level
apptainer exec \
    --bind "$SCRIPT_DIR:/mnt/script_dir" \
    --bind "$GTDBTK_OUTPUT_DIR:/mnt/gtdbtk_output_dir" \
    "$python_sif" \
    python "/mnt/script_dir/StrainCascade_extract_organism_name_gtdbtk.py" \
    --output_dir "/mnt/gtdbtk_output_dir" \
    --sample_name "$SAMPLE_NAME"

# Copy output files to taxonomic classification directory
log "$LOGS_DIR" "$LOG_NAME" "Copying classification files to $TAXONOMIC_CLASS_DIR"

# Copy summary files
if cp "$CLASSIFY_DIR"/*summary.tsv "$TAXONOMIC_CLASS_DIR/" 2>/dev/null; then
    log "$LOGS_DIR" "$LOG_NAME" "Copied summary files (.tsv) from $CLASSIFY_DIR"
else
    log "$LOGS_DIR" "$LOG_NAME" "No summary files found in $CLASSIFY_DIR"
fi

# Copy tree files
if cp "$CLASSIFY_DIR"/*.tree* "$TAXONOMIC_CLASS_DIR/" 2>/dev/null; then
    log "$LOGS_DIR" "$LOG_NAME" "Copied tree files (.tree) from $CLASSIFY_DIR"
else
    log "$LOGS_DIR" "$LOG_NAME" "No tree files found in $CLASSIFY_DIR"
fi

# Copy organism name files
if cp "$GTDBTK_OUTPUT_DIR"/*organism_name.txt "$TAXONOMIC_CLASS_DIR/" 2>/dev/null; then
    log "$LOGS_DIR" "$LOG_NAME" "Copied organism name files from $GTDBTK_OUTPUT_DIR"
else
    log "$LOGS_DIR" "$LOG_NAME" "No organism name files found in $GTDBTK_OUTPUT_DIR"
fi

# Copy level organism name file
if cp "$GTDBTK_OUTPUT_DIR"/level_organism_name.txt "$TAXONOMIC_CLASS_DIR/" 2>/dev/null; then
    log "$LOGS_DIR" "$LOG_NAME" "Copied level organism name file from $GTDBTK_OUTPUT_DIR"
else
    log "$LOGS_DIR" "$LOG_NAME" "No level organism name file found in $GTDBTK_OUTPUT_DIR"
fi

# Clean up
rm -rf "$GTDBTK_INPUT_DIR"

# Process results with R
tree_file=$(find "$TAXONOMIC_CLASS_DIR" -name '*gtdbtk.bac120.classify*' -print -quit)
tsv_file=$(find "$TAXONOMIC_CLASS_DIR" -name '*gtdbtk.bac120.summary.tsv' -print -quit)

if [[ -n "$tree_file" && -n "$tsv_file" ]]; then
    
    # Log R processing configuration
    log "$LOGS_DIR" "$LOG_NAME" "Processing GTDB-Tk taxonomy results with R"

    apptainer exec \
        --bind "$R_SCRIPT_DIR:/mnt/r_script_dir" \
        --bind "$TAXONOMIC_CLASS_DIR:/mnt/input_gtdbtk" \
        --bind "$QS_FILES_DIR:/mnt/output" \
        --bind "$RESULTS_INTEGRATION_DIR:/mnt/output2" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_gtdbtk_taxonomy.R" \
        --output_dir "/mnt/output" \
        --output_dir2 "/mnt/output2" \
        --tree "/mnt/input_gtdbtk/$(basename "$tree_file")" \
        --tsv "/mnt/input_gtdbtk/$(basename "$tsv_file")" \
        --sample_name "$SAMPLE_NAME" \
        --version "$VERSION" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
        exit 1
    }
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: Missing required files for R processing"
fi

log "$LOGS_DIR" "$LOG_NAME" "GTDB-Tk taxonomy analysis completed successfully"