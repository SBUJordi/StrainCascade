#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_ISEScan_IS_elements_identification.sh
# Description: Identifies insertion sequence (IS) elements using ISEScan as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <functional_analysis_main_abs> 
          <results_integration_abs> <version>
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
readonly THREADS="$8"
readonly GENOME_ASSEMBLY_DIR="$9"
readonly FUNCTIONAL_ANALYSIS_DIR="${10}"
readonly RESULTS_INTEGRATION_DIR="${11}"
readonly VERSION="${12}"

# Derived constants
readonly ISESCAN_OUTPUT_DIR="$OUTPUT_DIR/ISEScan_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"

# Source utility functions
source "$UTILS_FILE"

# Create required directories
create_directory "$ISESCAN_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Find required Apptainer images
readonly straincascade_crisprcas_phage_is_elements_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_crisprcas_phage_is_elements*.sif')
readonly r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
[[ -z "$analysis_assembly_file" ]] && {
    echo "Error: No assembly files found. Skipping ISEScan identification."
    exit 0
}

# Log analysis start
log "$LOGS_DIR" "$LOG_NAME" "Starting ISEScan analysis for $analysis_assembly_file"

# Build ISEScan command
isescan_cmd="isescan.py \
    --seqfile /mnt/output/temp/$(basename "$analysis_assembly_file") \
    --output /mnt/output \
    --nthread $THREADS"

# Run ISEScan identification
apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$ISESCAN_OUTPUT_DIR":/mnt/output \
    "$straincascade_crisprcas_phage_is_elements_sif" \
    bash -c "mkdir -p /mnt/output/temp && \
             cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/output/temp/$(basename "$analysis_assembly_file") && \
             source /opt/conda/etc/profile.d/conda.sh && \
             conda activate isescan_env && \
             $isescan_cmd" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: ISEScan analysis failed for $analysis_assembly_file"
    exit 1
}

# Copy output files to functional analysis directory
readonly OUTPUT_EXTENSIONS=(.gff .tsv .is.fna .orf.fna .orf.faa)

base_filename="$(basename "$analysis_assembly_file")"
base_filename="${base_filename%.*}"

for ext in "${OUTPUT_EXTENSIONS[@]}"; do
    if files=$(find "$ISESCAN_OUTPUT_DIR/temp" -type f -name "${base_filename}*${ext}"); then
        for file in $files; do
            cp "$file" "$FUNCTIONAL_ANALYSIS_DIR/ISEScan_${SAMPLE_NAME}_${ext}"
        done
    else
        echo "Warning: No *${ext} files found in $ISESCAN_OUTPUT_DIR/temp"
    fi
done

# Process results with R
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Find required output files
tsv_file=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ISEScan_*.tsv' -print -quit)
is_fna_file=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ISEScan_*.is.fna' -print -quit)
orf_fna_file=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ISEScan_*.orf.fna' -print -quit)
orf_faa_file=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name 'ISEScan_*.orf.faa' -print -quit)

# Validate required files exist
for file in "$tsv_file" "$is_fna_file" "$orf_fna_file" "$orf_faa_file"; do
    if [[ ! -f "$file" ]]; then
        echo "Warning: Required file $file not found. Skipping R processing."
        exit 0
    fi
done

# Log R processing configuration
cat << EOF
Processing ISEScan results with R using the following configuration:
TSV file: $tsv_file
IS DNA sequences: $is_fna_file
ORF DNA sequences: $orf_fna_file
ORF amino acid sequences: $orf_faa_file
Output directory: $QS_FILES_DIR
StrainCascade version v$VERSION
EOF

# Run R script for processing
apptainer exec \
    --bind "$R_SCRIPT_DIR":/mnt/r_script_dir \
    --bind "$FUNCTIONAL_ANALYSIS_DIR":/mnt/input \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_isescan.R" \
    --output_dir "/mnt/output" \
    --tsv "/mnt/input/$(basename "$tsv_file")" \
    --is_fna "/mnt/input/$(basename "$is_fna_file")" \
    --orf_fna "/mnt/input/$(basename "$orf_fna_file")" \
    --orf_faa "/mnt/input/$(basename "$orf_faa_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "ISEScan analysis completed successfully"