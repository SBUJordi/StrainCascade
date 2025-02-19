#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_Prokka_annotation.sh
# Description: Performs Prokka genome annotation as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <genome_annotation_main_abs>
          <results_integration_abs> <locus_tag> <force_overwrite> <version> <reproducibility_mode>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 15 ]] || show_usage

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
readonly GENOME_ANNOTATION_DIR="${10}"
readonly RESULTS_INTEGRATION_DIR="${11}"
readonly LOCUS_TAG="${12}"
readonly FORCE_OVERWRITE="${13}"
readonly VERSION="${14}"
readonly REPRODUCIBILITY_MODE="${15}"

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
readonly PROKKA_OUTPUT_DIR="$OUTPUT_DIR/Prokka_annotation_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly OUTPUT_EXTENSIONS=(_annotation_prokka.{gff,faa,ffn,tsv})

# Create required directories
create_directory "$PROKKA_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Working variables for processing
temp_input_file="temp_input_assembly_Prokka.fasta"
straincascade_genome_annotation_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_annotation*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Generate locus tag if set to automatic
final_locus_tag="$LOCUS_TAG"
if [[ "$LOCUS_TAG" == "automatic" ]]; then
    current_date=$(date +%y%m%d)
    random_digits=$(shuf -i 0-9 -n 3 | tr -d '\n')
    final_locus_tag="SC${random_digits}P${current_date}"
fi

# Build Prokka command
prokka_base_cmd="prokka --outdir /mnt/output \
                       --prefix ${SAMPLE_NAME}_annotation_prokka \
                       --cpus $THREADS \
                       --compliant \
                       --addgenes \
                       --rfam \
                       --mincontiglen 200 \
                       --locustag $final_locus_tag"

[[ "$FORCE_OVERWRITE" == "yes" ]] && prokka_base_cmd+=" --force"
prokka_cmd="$prokka_base_cmd /mnt/input/$temp_input_file"

# Run Prokka annotation
log "$LOGS_DIR" "$LOG_NAME" "Running Prokka annotation for $analysis_assembly_file"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$PROKKA_OUTPUT_DIR":/mnt/output \
    --bind "$ENTROPY_FILE":/dev/random \
    --bind "$ENTROPY_FILE":/dev/urandom \
    "$straincascade_genome_annotation_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate prokka_env && \
             cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/$temp_input_file && \
             export TMPDIR=/tmp && \
             $prokka_cmd && \
             rm /mnt/input/$temp_input_file" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Prokka annotation failed for $analysis_assembly_file. Skipping Prokka annotation."
    exit 0
}

# Copy output files to genome annotation directory
for ext in "${OUTPUT_EXTENSIONS[@]}"; do
    if files=$(find "$PROKKA_OUTPUT_DIR" -type f -name "*$ext"); then
        for file in $files; do
            cp "$file" "$GENOME_ANNOTATION_DIR"
        done
    else
        echo "Warning: No *$ext files found in $PROKKA_OUTPUT_DIR"
    fi
done

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Find required output files
gff_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_prokka.gff' -print -quit)
ffn_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_prokka.ffn' -print -quit)
faa_file=$(find "$GENOME_ANNOTATION_DIR" -name '*_annotation_prokka.faa' -print -quit)

# Validate required files exist
for file in "$gff_file" "$ffn_file" "$faa_file"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: Required file $file not found. Skipping processing of Prokka annotation."
        exit 0
    fi
done

# Log R processing configuration
log "$LOGS_DIR" "$LOG_NAME" "Processing Prokka results with R"

# Run R script for processing
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$GENOME_ANNOTATION_DIR":/mnt/input_prokka \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input_fasta \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_prokka.R" \
    --output_dir "/mnt/output" \
    --gff "/mnt/input_prokka/$(basename "$gff_file")" \
    --ffn "/mnt/input_prokka/$(basename "$ffn_file")" \
    --faa "/mnt/input_prokka/$(basename "$faa_file")" \
    --fasta "/mnt/input_fasta/$(basename "$analysis_assembly_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed. Skipping processing of Prokka annotation."
    exit 0
}

log "$LOGS_DIR" "$LOG_NAME" "Prokka analysis completed successfully"