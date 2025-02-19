#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_GTDB-Tk_de_novo_tree.sh
# Description: Performs GTDB-Tk de novo phylogenetic tree construction as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <external_assembly_dir> <threads> <genome_assembly_main_abs>
          <taxonomic_classification_main_abs> <databases_dir>
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
readonly EXTERNAL_ASSEMBLY_DIR="$8"
readonly THREADS="$9"
readonly GENOME_ASSEMBLY_DIR="${10}"
readonly TAXONOMIC_CLASS_DIR="${11}"
readonly DATABASES_DIR="${12}"

# Derived constants
readonly DE_NOVO_OUTPUT_DIR="$OUTPUT_DIR/GTDB-Tk_de_novo_tree_results"
readonly PROGRESS_FILE="$DE_NOVO_OUTPUT_DIR/progress.log"
readonly DE_NOVO_INPUT_DIR="$DE_NOVO_OUTPUT_DIR/input"
readonly ASSEMBLY_EXTENSION=".fasta"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$DE_NOVO_OUTPUT_DIR"
create_directory "$DE_NOVO_INPUT_DIR"

# Find required container
straincascade_taxonomic_functional_analysis_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_taxonomic_functional_analysis*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$analysis_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found in $GENOME_ASSEMBLY_DIR"
    exit 1
fi

# Copy analysis assembly file
cp "$analysis_assembly_file" "$DE_NOVO_INPUT_DIR/" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to copy analysis assembly file"
    exit 1
}

# Copy additional user-provided fasta files if directory exists
if [[ -d "$EXTERNAL_ASSEMBLY_DIR" ]]; then
    find "$EXTERNAL_ASSEMBLY_DIR" -name "*.fasta" -type f -exec cp {} "$DE_NOVO_INPUT_DIR/" \; || {
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Failed to copy some user-provided fasta files for de novo tree generation"
    }
    log "$LOGS_DIR" "$LOG_NAME" "Copied user-provided fasta files from $EXTERNAL_ASSEMBLY_DIR"
fi

# Get phylum for outgroup from GTDB-Tk results
readonly DEFAULT_PHYLUM="p__Bacillota"
readonly PHYLUM_FILE=$(find "$OUTPUT_DIR/GTDB-Tk_taxonomy_results" -type f -name "*gtdbtk_phylum_name.txt")

if [[ ! -f "$PHYLUM_FILE" ]] || [[ ! -s "$PHYLUM_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Warning: Valid phylum file not found or empty at $PHYLUM_FILE, using default: $DEFAULT_PHYLUM"
    outgroup="$DEFAULT_PHYLUM"
else
    outgroup=$(cat "$PHYLUM_FILE")
fi

# Find latest GTDB-Tk database
latest_db_dir=$(find "$DATABASES_DIR/gtdbtk_db" -type d -name 'release*' | sort -V | tail -n 1)
if [[ -z "$latest_db_dir" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No GTDB-Tk database found in $DATABASES_DIR/gtdbtk_db. Skipping GTDB-Tk de novo tree generation."
    exit 0
fi

# Log configuration
log "$LOGS_DIR" "$LOG_NAME" "Starting GTDB-Tk de novo tree analysis with configuration:
Input assembly: ${analysis_assembly_file}
Database: ${latest_db_dir}
Outgroup taxon: ${outgroup}
Track progress in: ${PROGRESS_FILE}"


# Execute GTDB-Tk de_novo_wf
apptainer exec \
    --bind "$DE_NOVO_INPUT_DIR":/mnt/input \
    --bind "$DE_NOVO_OUTPUT_DIR":/mnt/output \
    --bind "$latest_db_dir":/mnt/gtdbtk_db \
    "$straincascade_taxonomic_functional_analysis_sif" \
    bash -c "set -euo pipefail && \
             source /opt/conda/etc/profile.d/conda.sh && \
             conda activate gtdbtk_env && \
             export GTDBTK_DATA_PATH='/mnt/gtdbtk_db' && \
             gtdbtk de_novo_wf \
                 --genome_dir /mnt/input \
                 --outgroup_taxon '${outgroup}' \
                 --out_dir /mnt/output \
                 --prefix ${SAMPLE_NAME}_de_novo_tree_gtdbtk \
                 --extension ${ASSEMBLY_EXTENSION} \
                 --cpus ${THREADS} \
                 --bacteria > '/mnt/output/progress.log' 2>&1" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: GTDB-Tk de_novo_wf failed. Skipping GTDB-Tk de novo tree generation."
    exit 0
}

# Copy results to taxonomic classification directory
if ! cp "$DE_NOVO_OUTPUT_DIR"/*.decorated.tree* "$TAXONOMIC_CLASS_DIR/" 2>/dev/null; then
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No decorated tree files found to copy. Skipping GTDB-Tk de novo tree generation."
    exit 0
fi

if ! cp "$DE_NOVO_OUTPUT_DIR"/*.decorated.tree-taxonomy* "$TAXONOMIC_CLASS_DIR/" 2>/dev/null; then
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No taxonomy files for decorated tree found to copy. Skipping GTDB-Tk de novo tree generation."
    exit 0
fi

# Clean up only if decorated tree files exist in destination
if ls "$TAXONOMIC_CLASS_DIR"/*.decorated.tree* >/dev/null 2>&1; then
    rm -rf "$DE_NOVO_INPUT_DIR"
    log "$LOGS_DIR" "$LOG_NAME" "Cleaned up input directory"
fi

log "$LOGS_DIR" "$LOG_NAME" "GTDB-Tk de novo tree analysis completed successfully"