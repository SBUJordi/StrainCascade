#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_MicrobeAnnotator_annotation.sh
# Description: Performs MicrobeAnnotator functional annotation as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_annotation_main_abs> <functional_analysis_main_abs> 
          <results_integration_abs> <databases_dir> <version> <reproducibility_mode>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 14 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly OUTPUT_DIR="$6"
readonly SAMPLE_NAME="$7"
readonly INITIAL_THREADS="$8"
readonly GENOME_ANNOTATION_DIR="$9"
readonly FUNCTIONAL_ANALYSIS_DIR="${10}"
readonly RESULTS_INTEGRATION_DIR="${11}"
readonly DATABASES_DIR="${12}"
readonly VERSION="${13}"
readonly REPRODUCIBILITY_MODE="${14}"

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
readonly MICROBEANNOTATOR_OUTPUT_DIR="$OUTPUT_DIR/MicrobeAnnotator_annotation_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Create required directories
create_directory "$MICROBEANNOTATOR_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Find required container images
readonly straincascade_genome_annotation_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_annotation*.sif')
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find input protein file (.faa) with precedence order
faa_file_to_use=""
file_origin=""

# Try Bakta output first
if bakta_faa=$(find "$GENOME_ANNOTATION_DIR" -type f -name "*annotation_bakta.faa" -print -quit) && [[ -n "$bakta_faa" ]]; then
    faa_file_to_use="$bakta_faa"
    file_origin="bakta"
    log "$LOGS_DIR" "$LOG_NAME" "Using Bakta-generated .faa file: $bakta_faa"
# Try Prokka output second
elif prokka_faa=$(find "$GENOME_ANNOTATION_DIR" -type f -name "*annotation_prokka.faa" -print -quit) && [[ -n "$prokka_faa" ]]; then
    faa_file_to_use="$prokka_faa"
    file_origin="prokka"
    log "$LOGS_DIR" "$LOG_NAME" "Using Prokka-generated .faa file: $prokka_faa"
# Try any .faa file as last resort
elif any_faa=$(find "$GENOME_ANNOTATION_DIR" -type f -name "*.faa" -print -quit) && [[ -n "$any_faa" ]]; then
    faa_file_to_use="$any_faa"
    file_origin="non_StrainCascade_tool"
    log "$LOGS_DIR" "$LOG_NAME" "Using non-StrainCascade .faa file: $any_faa (Consider running Bakta/Prokka first)"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: No suitable .faa file found. Please run Bakta/Prokka first."
    exit 0
fi

# Prepare standardized input filename
readonly INPUT_FAA="${SAMPLE_NAME}_${file_origin}_annotation_microbeannotator.faa"

# Run MicrobeAnnotator annotation
log "$LOGS_DIR" "$LOG_NAME" "Running MicrobeAnnotator annotation"

# Execute MicrobeAnnotator in container
apptainer exec \
    --bind "$(dirname "$faa_file_to_use")":/mnt/input \
    --bind "$MICROBEANNOTATOR_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/microbeannotator_db":/mnt/microbeannotator_db \
    --bind "$ENTROPY_FILE":/dev/random \
    --bind "$ENTROPY_FILE":/dev/urandom \
    "$straincascade_genome_annotation_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate microbeannotator_env && \
             cp /mnt/input/$(basename "$faa_file_to_use") /mnt/input/$INPUT_FAA && \
             export TMPDIR=/tmp && \
             microbeannotator -i /mnt/input/$INPUT_FAA \
                             -d /mnt/microbeannotator_db \
                             -o /mnt/output \
                             -m blast \
                             -p 1 \
                             -t $THREADS \
                             --refine && \
             rm /mnt/input/$INPUT_FAA" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: MicrobeAnnotator annotation failed. Skipping MicrobeAnnotator annotation."
    exit 0
}

# Process output directories
find "$MICROBEANNOTATOR_OUTPUT_DIR" -mindepth 1 -type d | while read -r dir; do
    mv "$dir" "${dir}_microbeannotator"
done

# Clean up filenames
find "$MICROBEANNOTATOR_OUTPUT_DIR" -name '*__*' -exec bash -c 'mv "$1" "${1//__/_}"' _ {} \;

# Copy output files to appropriate directories
if output_files=$(find "$MICROBEANNOTATOR_OUTPUT_DIR" -mindepth 1 -name "*microbeannotator.faa.annot" -type f); then
    for file in $output_files; do
        cp "$file" "$GENOME_ANNOTATION_DIR/"
    done
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No annotation output files found"
fi

if metabolic_files=$(find "$MICROBEANNOTATOR_OUTPUT_DIR" -mindepth 1 -name "metabolic_summary*" -type f); then
    for file in $metabolic_files; do
        cp "$file" "$FUNCTIONAL_ANALYSIS_DIR/"
    done
else
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No metabolic summary files found"
fi

# Process results with R
log "$LOGS_DIR" "$LOG_NAME" "Processing MicrobeAnnotator results with R"

# Find required input files for R processing
annot_file=$(find "$GENOME_ANNOTATION_DIR" -name '*microbeannotator.faa.annot' -print -quit)
tab_file=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name '*module_completeness.tab' -print -quit)

# Process annotation results
if [[ -n "$annot_file" ]]; then
    apptainer exec \
        --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
        --bind "$GENOME_ANNOTATION_DIR":/mnt/input \
        --bind "$QS_FILES_DIR":/mnt/output \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_microbeannotator.R" \
        --output_dir "/mnt/output" \
        --annot "/mnt/input/$(basename "$annot_file")" \
        --version "$VERSION" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed (annotation). Skipping processing of MicrobeAnnotator annotation."
            exit 0
        }
fi

# Process KEGG module results 
if [[ -n "$tab_file" ]]; then
    apptainer exec \
        --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
        --bind "$FUNCTIONAL_ANALYSIS_DIR":/mnt/input \
        --bind "$QS_FILES_DIR":/mnt/output \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_microbeannotator_KEGG_modules.R" \
        --output_dir "/mnt/output" \
        --tab "/mnt/input/$(basename "$tab_file")" \
        --version "$VERSION" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed (KEGG modules). Skipping processing of MicrobeAnnotator pathway analysis."
            exit 0
        }
fi

log "$LOGS_DIR" "$LOG_NAME" "MicrobeAnnotator analysis completed successfully"