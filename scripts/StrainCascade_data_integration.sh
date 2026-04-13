#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_data_integration.sh
# Description: Integrates and processes genome annotation results from multiple tools
#             (Bakta, Prokka, MicrobeAnnotator) and generates comprehensive reports

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> 
          <genome_assembly_main_abs> <results_integration_abs> <version> <sample_name>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 9 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly GENOME_ASSEMBLY_DIR="$6"
readonly RESULTS_DIR="$7"
readonly VERSION="$8"
readonly SAMPLE_NAME="$9"

# Derived constants
readonly QS_FILES_DIR="$RESULTS_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$QS_FILES_DIR"

# Find required Apptainer image
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Log script initialization
log "$LOGS_DIR" "$LOG_NAME" "Starting StrainCascade data integration for sample: $SAMPLE_NAME"

# Retrieve and validate assembly file
readonly ANALYSIS_ASSEMBLY_FILE=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$ANALYSIS_ASSEMBLY_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping data integration."
    exit 0
fi
log "$LOGS_DIR" "$LOG_NAME" "Using assembly file: $ANALYSIS_ASSEMBLY_FILE"

# =============================================================================
# STEP 0: NORMALIZE FINAL ASSEMBLY (if not already done)
# =============================================================================
# Normalization is normally performed at the end of assembly polishing (SC12).
# This step is a FALLBACK for cases where:
# - Assembly-only input was provided (polishing was skipped)
# - Normalization failed during polishing
# - Pipeline was run with a custom module selection that skipped polishing
#
# If contig_name_mapping.tsv already exists, normalization was already performed.

if [[ -f "${QS_FILES_DIR}/contig_name_mapping.tsv" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Contig names already normalized (mapping file exists). Skipping normalization."
else
    log "$LOGS_DIR" "$LOG_NAME" "Contig name mapping not found - running normalization (fallback for assembly-only input)"
    
    # Check if Circlator results exist
    CIRCLATOR_QS=$(find "$QS_FILES_DIR" -name 'circlator_results.qs' -print -quit 2>/dev/null || true)
    circlator_arg=""
    if [[ -n "$CIRCLATOR_QS" && -f "$CIRCLATOR_QS" ]]; then
        circlator_arg="--circlator_qs /mnt/qs_files/$(basename "$CIRCLATOR_QS")"
        log "$LOGS_DIR" "$LOG_NAME" "Will update Circlator results with normalized names"
    fi

    apptainer exec \
        --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
        --bind "$(dirname "$ANALYSIS_ASSEMBLY_FILE")":/mnt/assembly \
        --bind "${QS_FILES_DIR}:/mnt/qs_files" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_normalize_final_assembly.R" \
        --fasta "/mnt/assembly/$(basename "$ANALYSIS_ASSEMBLY_FILE")" \
        --output_dir "/mnt/qs_files" \
        $circlator_arg || {
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Assembly normalization failed. Continuing with original names."
    }
    
    # Copy mapping file to genome assembly directory if created
    if [[ -f "${QS_FILES_DIR}/contig_name_mapping.tsv" ]]; then
        cp "${QS_FILES_DIR}/contig_name_mapping.tsv" "$(dirname "$ANALYSIS_ASSEMBLY_FILE")/contig_name_mapping.tsv"
        log "$LOGS_DIR" "$LOG_NAME" "Contig name mapping saved to genome assembly directory"
    fi
fi

# =============================================================================
# STEP 1: PROCESS SELECTED ASSEMBLY
# =============================================================================
# Process the (normalized) assembly to create selected_assembly_results.qs

log "$LOGS_DIR" "$LOG_NAME" "Processing selected assembly"

apptainer exec \
    --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
    --bind "$(dirname "$ANALYSIS_ASSEMBLY_FILE")":/mnt/assembly \
    --bind "${QS_FILES_DIR}:/mnt/output" \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_selected_assembly.R" \
    --fasta "/mnt/assembly/$(basename "$ANALYSIS_ASSEMBLY_FILE")" \
    --output_dir "/mnt/output" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Warning: Selected assembly processing failed"
}

# =============================================================================
# STEP 2: FIND AND AGGREGATE ANNOTATION RESULTS
# =============================================================================

# Find annotation result files
readonly BAKTA_QS=$(find "$QS_FILES_DIR" -name 'bakta_results.qs' -print -quit)
readonly PROKKA_QS=$(find "$QS_FILES_DIR" -name 'prokka_results.qs' -print -quit)
readonly MICROBEANNOTATOR_QS=$(find "$QS_FILES_DIR" -name 'microbeannotator_results.qs' -print -quit)
readonly DEEPFRI_QS=$(find "$QS_FILES_DIR" -name 'deepfri_results.qs' -print -quit)

# Build R script arguments
r_args=""
[[ -n "$BAKTA_QS" ]] && r_args+=" --bakta $(basename "$BAKTA_QS")"
[[ -n "$PROKKA_QS" ]] && r_args+=" --prokka $(basename "$PROKKA_QS")"
[[ -n "$MICROBEANNOTATOR_QS" ]] && r_args+=" --microbeannotator $(basename "$MICROBEANNOTATOR_QS")"
[[ -n "$DEEPFRI_QS" ]] && r_args+=" --deepfri $(basename "$DEEPFRI_QS")"

# Log annotation files found
log "$LOGS_DIR" "$LOG_NAME" "Found annotation files:
Bakta: ${BAKTA_QS:-not available}
Prokka: ${PROKKA_QS:-not available}
MicrobeAnnotator: ${MICROBEANNOTATOR_QS:-not available}
DeepFRI: ${DEEPFRI_QS:-not available}"

# Process annotation results if available
if [[ -n "$r_args" ]]; then
    # Run genome annotation aggregation
    log "$LOGS_DIR" "$LOG_NAME" "Running genome annotation aggregation"
    
    apptainer exec \
        --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
        --bind "${QS_FILES_DIR}:/mnt/input_output" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_genome_annotation_aggregation.R" \
        --output_dir "/mnt/input_output" \
        $r_args \
        --version "$VERSION"

    # Run genome annotation integration if aggregation successful
    if [[ -f "${QS_FILES_DIR}/annotation_results_aggregated.qs" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Running genome annotation integration"
        
        # Check if contig name mapping file exists
        mapping_arg=""
        if [[ -f "${QS_FILES_DIR}/contig_name_mapping.tsv" ]]; then
            mapping_arg="--mapping_file /mnt/input_output/contig_name_mapping.tsv"
            log "$LOGS_DIR" "$LOG_NAME" "Using contig name mapping for annotation integration"
        fi
        
        apptainer exec \
            --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
            --bind "${QS_FILES_DIR}:/mnt/input_output" \
            "$r_sif" \
            Rscript "/mnt/r_script_dir/R_genome_annotation_integration.R" \
            --input_file "/mnt/input_output/annotation_results_aggregated.qs" \
            --output_dir "/mnt/input_output" \
            $mapping_arg \
            --version "$VERSION"
    else
        log "$LOGS_DIR" "$LOG_NAME" "Error: Aggregation failed. Skipping integration step."
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "No input files available for aggregation and integration"
fi

# Clean up previous runs
rm -f "${RESULTS_DIR}/*.RData"

# Check if contig name mapping file exists for RData generation
qs2rdata_mapping_arg=""
if [[ -f "${QS_FILES_DIR}/contig_name_mapping.tsv" ]]; then
    qs2rdata_mapping_arg="--mapping_file /mnt/input/contig_name_mapping.tsv"
    log "$LOGS_DIR" "$LOG_NAME" "Using contig name mapping for RData generation"
fi

# Generate final RData file
log "$LOGS_DIR" "$LOG_NAME" "Generating final RData file"

apptainer exec \
    --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
    --bind "${QS_FILES_DIR}:/mnt/input" \
    --bind "${RESULTS_DIR}:/mnt/output" \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_qs2RData.R" \
    --sample_name "${SAMPLE_NAME}" \
    --input_dir "/mnt/input" \
    --output_dir "/mnt/output" \
    $qs2rdata_mapping_arg || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Final data integration failed"
    exit 0
}

# Locate and validate RData file
readonly RDATA_FILE=$(find "$RESULTS_DIR" -name "StrainCascade_Results_${SAMPLE_NAME}*.RData" -print -quit)
if [[ -z "$RDATA_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: RData file not found for sample ${SAMPLE_NAME}"
    exit 0
fi

# Generate final report
log "$LOGS_DIR" "$LOG_NAME" "Generating StrainCascade report"

apptainer exec \
    --bind "${R_SCRIPT_DIR}:/mnt/r_script_dir" \
    --bind "${RESULTS_DIR}:/mnt/input_output" \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/StrainCascade_report_generation.R" \
    --sample_name "${SAMPLE_NAME}" \
    --RData "/mnt/input_output/$(basename "$RDATA_FILE")" \
    --output_dir "/mnt/input_output" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: StrainCascade report generation failed"
    rm -f "$RESULTS_DIR"/StrainCascade_Analysis_Report_*.Rmd 2>/dev/null || true
    exit 0
}

# Remove any existing StrainCascade Analysis Report Rmd files
rm -f "$RESULTS_DIR"/StrainCascade_Analysis_Report_*.Rmd 2>/dev/null || true

log "$LOGS_DIR" "$LOG_NAME" "StrainCascade data integration completed successfully"