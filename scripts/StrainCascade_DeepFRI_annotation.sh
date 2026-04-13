#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_DeepFRI_annotation.sh
# Description: Performs DeepFRI functional annotation using structure-based predictions

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <genome_annotation_main_abs> <results_integration_abs> 
          <version> <reproducibility_mode>
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
readonly GENOME_ANNOTATION_DIR="$8"
readonly RESULTS_INTEGRATION_DIR="$9"
readonly VERSION="${10}"
readonly REPRODUCIBILITY_MODE="${11}"

# Source utility functions
source "$UTILS_FILE"

# Derived constants
readonly DEEPFRI_OUTPUT_DIR="$OUTPUT_DIR/DeepFRI_annotation_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"

# Create required directories
create_directory "$DEEPFRI_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Find required Apptainer image (combined genome annotation container)
genome_annotation_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_annotation*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find input .faa file (prefer Bakta, fallback to Prokka)
faa_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_bakta.faa' -print -quit)
if [[ -z "$faa_file" ]]; then
    faa_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_prokka.faa' -print -quit)
fi

if [[ -z "$faa_file" || ! -f "$faa_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No protein FASTA file (.faa) found. Skipping DeepFRI annotation."
    exit 0
fi

log "$LOGS_DIR" "$LOG_NAME" "Using protein file: $faa_file"

# Validate .faa file is non-empty
if [[ ! -s "$faa_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Protein FASTA file is empty. Skipping DeepFRI annotation."
    exit 0
fi

# Set random seed for deterministic mode
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    export PYTHONHASHSEED=42
    log "$LOGS_DIR" "$LOG_NAME" "Running in deterministic mode with PYTHONHASHSEED=42"
fi

# Define ontologies to process
readonly ONTOLOGIES=("mf" "bp" "cc" "ec")
readonly ONTOLOGY_NAMES=("MolecularFunction" "BiologicalProcess" "CellularComponent" "EnzymeCommission")

# Run DeepFRI for each ontology
log "$LOGS_DIR" "$LOG_NAME" "Running DeepFRI annotation for all ontologies"

for i in "${!ONTOLOGIES[@]}"; do
    ont="${ONTOLOGIES[$i]}"
    ont_name="${ONTOLOGY_NAMES[$i]}"
    
    log "$LOGS_DIR" "$LOG_NAME" "Processing ontology: $ont_name ($ont)"
    
    output_prefix="${SAMPLE_NAME}_deepfri_${ont}"
    
    apptainer exec \
        --bind "$(dirname "$faa_file")":/mnt/input \
        --bind "$DEEPFRI_OUTPUT_DIR":/mnt/output \
        "$genome_annotation_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate deepfri_env && \
                 cd /opt/DeepFRI && \
                 python predict.py \
                     --model_config /opt/DeepFRI/trained_models/model_config.json \
                     --fasta_fn '/mnt/input/$(basename "$faa_file")' \
                     -ont '$ont' \
                     --output_fn_prefix '/mnt/output/${output_prefix}'" || {
        log "$LOGS_DIR" "$LOG_NAME" "Warning: DeepFRI failed for ontology $ont. Continuing with other ontologies."
        continue
    }
    
    # Check if output was generated (DeepFRI uses uppercase ontology in output filename)
    ont_upper=$(echo "$ont" | tr '[:lower:]' '[:upper:]')
    output_csv="${DEEPFRI_OUTPUT_DIR}/${output_prefix}_${ont_upper}_predictions.csv"
    if [[ ! -f "$output_csv" || ! -s "$output_csv" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: No predictions generated for ontology $ont"
    else
        log "$LOGS_DIR" "$LOG_NAME" "Successfully generated predictions for ontology $ont"
    fi
done

# Copy output files to genome annotation directory
find "$DEEPFRI_OUTPUT_DIR" -type f -name "${SAMPLE_NAME}_deepfri_*.csv" -exec cp {} "$GENOME_ANNOTATION_DIR" \;

# Count generated predictions
prediction_files=$(find "$GENOME_ANNOTATION_DIR" -name "${SAMPLE_NAME}_deepfri_*_predictions.csv" 2>/dev/null | wc -l)

if [[ $prediction_files -eq 0 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No DeepFRI predictions were generated. Skipping R processing."
    exit 0
fi

log "$LOGS_DIR" "$LOG_NAME" "Generated predictions for $prediction_files ontologies"

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

log "$LOGS_DIR" "$LOG_NAME" "Processing DeepFRI results with R"

apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$GENOME_ANNOTATION_DIR":/mnt/input \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_deepfri.R" \
    --output_dir "/mnt/output" \
    --input_dir "/mnt/input" \
    --sample_name "$SAMPLE_NAME" \
    --faa "/mnt/input/$(basename "$faa_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed. Skipping processing of DeepFRI annotation."
    exit 0
}

log "$LOGS_DIR" "$LOG_NAME" "DeepFRI analysis completed successfully"
