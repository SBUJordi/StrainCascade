#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_Bakta_annotation.sh
# Description: Performs Bakta genome annotation as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <taxonomic_classification_main_abs> 
          <genome_annotation_main_abs> <results_integration_abs> <databases_dir> <locus_tag> 
          <force_overwrite> <version> <reproducibility_mode>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 17 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly OUTPUT_DIR="$6"
readonly SAMPLE_NAME="$7"
readonly INITIAL_THREADS="$8"
readonly GENOME_ASSEMBLY_DIR="$9"
readonly TAXONOMIC_CLASS_DIR="${10}"
readonly GENOME_ANNOTATION_DIR="${11}"
readonly RESULTS_INTEGRATION_DIR="${12}"
readonly DATABASES_DIR="${13}"
readonly LOCUS_TAG="${14}"
readonly FORCE_OVERWRITE="${15}"
readonly VERSION="${16}"
readonly REPRODUCIBILITY_MODE="${17}"

# Source utility functions
source "$UTILS_FILE"

# Set threads based on algorithm type
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    readonly THREADS=1
else
    readonly THREADS="${INITIAL_THREADS}"
fi

# Create deterministic entropy source
readonly ENTROPY_FILE="$OUTPUT_DIR/deterministic_entropy_file"
if [[ ! -f "$ENTROPY_FILE" ]]; then
    dd if=/dev/zero bs=1024 count=100 > "$ENTROPY_FILE"
    log "$LOGS_DIR" "StrainCascade.log" "Deterministic entropy file created at $ENTROPY_FILE"
fi

# Derived constants
readonly BAKTA_OUTPUT_DIR="$OUTPUT_DIR/Bakta_annotation_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly OUTPUT_EXTENSIONS=(_annotation_bakta.{tsv,ffn,faa,gbff,png,svg})
readonly TEMP_INPUT_FILE="temp_input_assembly_Bakta.fasta"

# Create required directories
create_directory "$BAKTA_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Working variables for processing
genome_annotation_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_annotation*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Get taxonomic information
# Get taxonomic level file
level_file=$(find "$TAXONOMIC_CLASS_DIR" -name '*level_organism_name.txt' -print -quit)
if [[ -f "$level_file" ]]; then
    taxonomic_level=$(cat "$level_file")
    log "$LOGS_DIR" "$LOG_NAME" "Using taxonomic level from: $(basename "$level_file")"
    
    # Check for additional matching files and log warning if found
    level_files=$(find "$TAXONOMIC_CLASS_DIR" -name '*level_organism_name.txt' | wc -l)
    if (( level_files > 1 )); then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Found $level_files taxonomic level files. Using the first one found."
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Could not find level_organism_name.txt file in $TAXONOMIC_CLASS_DIR. Skipping Bakta annotation."
    exit 0
fi

# Get organism name file
organism_name=""
organism_name_file=$(find "$TAXONOMIC_CLASS_DIR" -name '*gtdbtk_organism_name.txt' -print -quit)
if [[ -f "$organism_name_file" ]]; then
    organism_name=$(cat "$organism_name_file")
    log "$LOGS_DIR" "$LOG_NAME" "Using organism name from: $(basename "$organism_name_file")"
    
    # Check for additional matching files and log warning if found
    organism_files=$(find "$TAXONOMIC_CLASS_DIR" -name '*gtdbtk_organism_name.txt' | wc -l)
    if (( organism_files > 1 )); then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Found $organism_files organism name files. Using the first one found."
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No organism name file found in $TAXONOMIC_CLASS_DIR"
fi

organism_name=""
[[ -s "$organism_name_file" ]] && organism_name=$(cat "$organism_name_file")

# Generate locus tag if set to automatic
final_locus_tag="$LOCUS_TAG"
if [[ "$LOCUS_TAG" == "automatic" ]]; then
    current_date=$(date +%y%m%d)
    random_digits=$(shuf -i 0-9 -n 3 | tr -d '\n')
    final_locus_tag="SC${random_digits}B${current_date}"
fi

# Build Bakta command
bakta_base_cmd="bakta --db /mnt/bakta_db/db \
                      --output /mnt/output \
                      --prefix ${SAMPLE_NAME}_annotation_bakta \
                      --compliant \
                      --threads $THREADS \
                      --locus-tag $final_locus_tag \
                      --tmp-dir /tmp"

# Add taxonomic info if available
if [[ -n "$organism_name" ]]; then
    if [[ "$taxonomic_level" == "genus" ]]; then
        bakta_base_cmd+=" --genus \"$organism_name\""
    elif [[ "$taxonomic_level" == "species" ]]; then
        bakta_base_cmd+=" --species \"$organism_name\""
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No organism name found. Skipping taxonomic info."
fi

[[ "$FORCE_OVERWRITE" == "yes" ]] && bakta_base_cmd+=" --force"
bakta_cmd="$bakta_base_cmd /mnt/input/$TEMP_INPUT_FILE"

# Run Bakta annotation
log "$LOGS_DIR" "$LOG_NAME" "Running Bakta annotation for $analysis_assembly_file"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$BAKTA_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/bakta_db/db":/mnt/bakta_db/db \
    --bind "$ENTROPY_FILE":/dev/random \
    --bind "$ENTROPY_FILE":/dev/urandom \
    "$genome_annotation_sif" \
    bash -c "set -e; \
             source /opt/conda/etc/profile.d/conda.sh && \
             conda activate bakta_env && \
             cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/$TEMP_INPUT_FILE && \
             export TMPDIR=/tmp && \
             $bakta_cmd && \
             rm /mnt/input/$TEMP_INPUT_FILE" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Bakta annotation failed for $analysis_assembly_file. Skipping Bakta annotation."
    exit 0
}

# Copy output files to genome annotation directory
for ext in "${OUTPUT_EXTENSIONS[@]}"; do
    if files=$(find "$BAKTA_OUTPUT_DIR" -type f -name "*$ext"); then
        for file in $files; do
            cp "$file" "$GENOME_ANNOTATION_DIR"
        done
    else
        log "$LOGS_DIR" "$LOG_NAME" "Warning: No *$ext files found in $BAKTA_OUTPUT_DIR"
    fi
done

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Find required output files
tsv_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_bakta.tsv' -print -quit)
ffn_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_bakta.ffn' -print -quit)
faa_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_bakta.faa' -print -quit)

# Validate required files exist
for file in "$tsv_file" "$ffn_file" "$faa_file"; do
    if [[ ! -f "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Required file $file not found. Skipping processing of Bakta annotation."
        exit 0
    fi
done

# Log R processing configuration
log "$LOGS_DIR" "$LOG_NAME" "Processing Bakta results with R"

# Run R script for processing
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$GENOME_ANNOTATION_DIR":/mnt/input_bakta \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input_fasta \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_bakta.R" \
    --output_dir "/mnt/output" \
    --tsv "/mnt/input_bakta/$(basename "$tsv_file")" \
    --ffn "/mnt/input_bakta/$(basename "$ffn_file")" \
    --faa "/mnt/input_bakta/$(basename "$faa_file")" \
    --fasta "/mnt/input_fasta/$(basename "$analysis_assembly_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed. Skipping processing of Bakta annotation."
    exit 0
}

log "$LOGS_DIR" "$LOG_NAME" "Bakta analysis completed successfully"