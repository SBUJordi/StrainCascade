#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_ResFinder_antimicrobial_resistance_identification.sh
# Description: Performs ResFinder antimicrobial resistance identification as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <sequencing_type> <genome_assembly_main_abs> 
          <functional_analysis_main_abs> <results_integration_abs> <databases_dir> <version>
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
readonly THREADS="$8"
readonly SEQUENCING_TYPE="$9"
readonly GENOME_ASSEMBLY_DIR="${10}"
readonly FUNCTIONAL_ANALYSIS_DIR="${11}"
readonly RESULTS_INTEGRATION_DIR="${12}"
readonly DATABASES_DIR="${13}"
readonly VERSION="${14}"

# Derived constants
readonly RESFINDER_OUTPUT_DIR="$OUTPUT_DIR/ResFinder_AMR_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly TEMP_INPUT_FILE="temp_input_assembly_ResFinder.fasta"

# Source utility functions
source "$UTILS_FILE"

# Determine sequencing type flag for ResFinder
readonly RESFINDER_SEQ_TYPE=$(case $SEQUENCING_TYPE in
    *nano*) echo "--nanopore" ;;
    *) echo "" ;;
esac)

# Find required Apptainer images
readonly straincascade_taxonomic_functional_analysis_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_taxonomic_functional_analysis*.sif')
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Create required directories
create_directory "$RESFINDER_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Find analysis assembly file
readonly ANALYSIS_ASSEMBLY_FILE=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -z "$ANALYSIS_ASSEMBLY_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found. Skipping ResFinder identification."
    exit 0
fi

# Build ResFinder command
readonly RESFINDER_BASE_CMD="python /opt/conda/envs/resfinder_env/bin/run_resfinder.py \
    -ifa /mnt/input/$TEMP_INPUT_FILE \
    -o /mnt/output \
    -s 'Other' \
    --acquired \
    --point \
    --ignore_missing_species \
    -l 0.6 \
    -t 0.8 \
    -l_p 0.6 \
    -t_p 0.8 \
    -db_res /mnt/resfinder_db \
    -db_point /mnt/pointfinder_db \
    -db_disinf /mnt/disinfinder_db \
    --output_aln \
    -acq \
    -d \
    -u"

# Run ResFinder analysis
log "$LOGS_DIR" "$LOG_NAME" "Running ResFinder AMR identification for $ANALYSIS_ASSEMBLY_FILE"

apptainer exec \
    --bind "$(dirname "$ANALYSIS_ASSEMBLY_FILE")":/mnt/input \
    --bind "$RESFINDER_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/":/mnt \
    "$straincascade_taxonomic_functional_analysis_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate resfinder_env && \
             cp /mnt/input/$(basename "$ANALYSIS_ASSEMBLY_FILE") /mnt/input/$TEMP_INPUT_FILE && \
             $RESFINDER_BASE_CMD $RESFINDER_SEQ_TYPE && \
             rm /mnt/input/$TEMP_INPUT_FILE" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: ResFinder analysis failed for $ANALYSIS_ASSEMBLY_FILE"
    exit 1
}

# Rename and copy output files
log "$LOGS_DIR" "$LOG_NAME" "Processing ResFinder output files"

find "$RESFINDER_OUTPUT_DIR" -type f -print0 | while IFS= read -r -d '' file; do
    dir=$(dirname "$file")
    base=$(basename "$file")
    ext="${base##*.}"
    filename="${base%.*}"
    new_base="${filename}_${SAMPLE_NAME}_AMR_identification_ResFinder.${ext}"
    mv "$file" "$dir/$new_base"
done

# Copy relevant files to functional analysis directory
if ! find "$RESFINDER_OUTPUT_DIR" -type f \( -name "*.txt" -o -name "*.fsa" \) -exec cp {} "$FUNCTIONAL_ANALYSIS_DIR" \;; then
    log "$LOGS_DIR" "$LOG_NAME" "Warning: No suitable output files found in $RESFINDER_OUTPUT_DIR"
fi

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Find required output files
readonly TAB_FILE=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ResFinder_results_tab_*' -print -quit)
readonly PHENO_FILE=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'pheno_table_*' -print -quit)
readonly HIT_FSA_FILE=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ResFinder_Hit_in_genome_seq*' -print -quit)
readonly REF_FSA_FILE=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ResFinder_Resistance_gene_seq*' -print -quit)

# Validate required files exist
for file in "$TAB_FILE" "$PHENO_FILE" "$HIT_FSA_FILE" "$REF_FSA_FILE"; do
    if [[ ! -f "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Required file $file not found"
        exit 1
    fi
done

# Log R processing configuration
log "$LOGS_DIR" "$LOG_NAME" "Processing ResFinder results with R:
Tab file: $TAB_FILE
Pheno file: $PHENO_FILE
Hit FSA file: $HIT_FSA_FILE
Ref FSA file: $REF_FSA_FILE
FASTA file: $ANALYSIS_ASSEMBLY_FILE"

# Run R script for processing
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$FUNCTIONAL_ANALYSIS_DIR":/mnt/input \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_resfinder.R" \
    --output_dir "/mnt/output" \
    --tab_file "/mnt/input/$(basename "$TAB_FILE")" \
    --pheno_file "/mnt/input/$(basename "$PHENO_FILE")" \
    --hit_fsa "/mnt/input/$(basename "$HIT_FSA_FILE")" \
    --ref_fsa "/mnt/input/$(basename "$REF_FSA_FILE")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "ResFinder analysis completed successfully"