#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_geNomad_phage_identification.sh
# Description: Performs geNomad phage and plasmid identification as part of StrainCascade

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
readonly GENOMAD_OUTPUT_DIR="$OUTPUT_DIR/geNomad_phage_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly TEMP_INPUT_FILE="temp_input_assembly_geNomad.fasta"
readonly MIN_LENGTH=200

# Source utility functions
source "$UTILS_FILE"

# Create required directories
mkdir -p "$GENOMAD_OUTPUT_DIR" "$QS_FILES_DIR"

# Find Apptainer images
straincascade_crisprcas_phage_is_elements_sif=$(find "$APPTAINER_DIR" -name 'straincascade_crisprcas_phage_is_elements*.sif' -print -quit)
r_sif=$(find "$APPTAINER_DIR" -name 'r_4.4.1*.sif' -print -quit)

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

[[ -z "$analysis_assembly_file" ]] && {
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping geNomad module."
    exit 0
}

# Run geNomad analysis
log "$LOGS_DIR" "$LOG_NAME" "Running geNomad phage and plasmid identification for $analysis_assembly_file"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$GENOMAD_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/genomad_db":/mnt/genomad_db \
    "$straincascade_crisprcas_phage_is_elements_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate genomad_env && \
             cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/$TEMP_INPUT_FILE && \
             genomad end-to-end \
                --cleanup \
                --splits 8 \
                --threads $THREADS \
                --min-score 0.7 \
                /mnt/input/$TEMP_INPUT_FILE \
                /mnt/output \
                /mnt/genomad_db && \
             rm /mnt/input/$TEMP_INPUT_FILE" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: geNomad analysis failed for $analysis_assembly_file"
    exit 0
}

# Process output files
log "$LOGS_DIR" "$LOG_NAME" "Processing geNomad output files"

# Get the output prefix (based on TEMP_INPUT_FILE without extension, since that's what geNomad uses)
input_base="${TEMP_INPUT_FILE%.fasta}"

# Define summary directory
summary_dir="$GENOMAD_OUTPUT_DIR/${input_base}_summary"

# Copy virus-related files to functional analysis directory
if [[ -d "$summary_dir" ]]; then
    for base_name in "${input_base}_virus_summary.tsv" "${input_base}_virus_genes.tsv" "${input_base}_virus.fna" "${input_base}_virus_proteins.faa"; do
        if [[ -f "$summary_dir/$base_name" ]]; then
            extension="${base_name##*.}"
            filename="${base_name%.*}"
            cp "$summary_dir/$base_name" "$FUNCTIONAL_ANALYSIS_DIR/${filename}_${SAMPLE_NAME}_genomad_identification.${extension}"
            log "$LOGS_DIR" "$LOG_NAME" "Processed virus output file: ${filename}_${SAMPLE_NAME}_genomad_identification.${extension}"
        else
            log "$LOGS_DIR" "$LOG_NAME" "Warning: Expected virus output file not found: $base_name"
        fi
    done

    # Copy plasmid-related files to functional analysis directory
    for base_name in "${input_base}_plasmid_summary.tsv" "${input_base}_plasmid_genes.tsv" "${input_base}_plasmid.fna" "${input_base}_plasmid_proteins.faa"; do
        if [[ -f "$summary_dir/$base_name" ]]; then
            extension="${base_name##*.}"
            filename="${base_name%.*}"
            cp "$summary_dir/$base_name" "$FUNCTIONAL_ANALYSIS_DIR/${filename}_${SAMPLE_NAME}_genomad_identification.${extension}"
            log "$LOGS_DIR" "$LOG_NAME" "Processed plasmid output file: ${filename}_${SAMPLE_NAME}_genomad_identification.${extension}"
        else
            log "$LOGS_DIR" "$LOG_NAME" "Warning: Expected plasmid output file not found: $base_name"
        fi
    done
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: geNomad summary directory not found: $summary_dir"
    exit 0
fi

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Find input files for R processing
virus_summary_tsv=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "*_virus_summary_${SAMPLE_NAME}_genomad_identification.tsv" -print -quit)
virus_genes_tsv=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "*_virus_genes_${SAMPLE_NAME}_genomad_identification.tsv" -print -quit)
plasmid_summary_tsv=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "*_plasmid_summary_${SAMPLE_NAME}_genomad_identification.tsv" -print -quit)
plasmid_genes_tsv=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name "*_plasmid_genes_${SAMPLE_NAME}_genomad_identification.tsv" -print -quit)

# Check if required files exist
if [[ ! -f "$virus_summary_tsv" ]] || [[ ! -f "$virus_genes_tsv" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Warning: Required virus output files missing. Continuing without R processing."
else
    log "$LOGS_DIR" "$LOG_NAME" "Processing geNomad results with R"

    # Run R script for processing
    apptainer exec \
        --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
        --bind "$FUNCTIONAL_ANALYSIS_DIR":/mnt/input \
        --bind "$QS_FILES_DIR":/mnt/output \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_genomad.R" \
        --output_dir "/mnt/output" \
        --virus_summary "/mnt/input/$(basename "$virus_summary_tsv")" \
        --virus_genes "/mnt/input/$(basename "$virus_genes_tsv")" \
        --plasmid_summary "${plasmid_summary_tsv:+/mnt/input/$(basename "$plasmid_summary_tsv")}" \
        --plasmid_genes "${plasmid_genes_tsv:+/mnt/input/$(basename "$plasmid_genes_tsv")}" \
        --version "$VERSION" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
        exit 0
    }
fi

log "$LOGS_DIR" "$LOG_NAME" "geNomad analysis completed successfully"
