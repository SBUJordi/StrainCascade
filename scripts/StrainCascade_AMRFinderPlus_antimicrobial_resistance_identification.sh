#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh
# Description: Performs AMRFinderPlus antimicrobial resistance identification as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_dir> 
          <sample_name> <threads> <genome_assembly_main_abs> <functional_analysis_main_abs>
          <results_integration_abs> <taxonomic_classification_main_abs> <databases_dir> <version>
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
readonly GENOME_ASSEMBLY_DIR="$9"
readonly TAXONOMIC_CLASS_DIR="${10}"
readonly FUNCTIONAL_ANALYSIS_DIR="${11}"
readonly RESULTS_INTEGRATION_DIR="${12}"
readonly DATABASES_DIR="${13}"
readonly VERSION="${14}"

# Derived constants
readonly AMRFINDER_OUTPUT_DIR="$OUTPUT_DIR/AMRFinderPlus_AMR_identification_results"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly TEMP_DIR="${AMRFINDER_OUTPUT_DIR}/temp_${SAMPLE_NAME}_$$"
readonly OUTPUT_BASE="${SAMPLE_NAME}_identification_amrfinderplus"
readonly MUTATION_BASE="${SAMPLE_NAME}_mutations_amrfinderplus"

# Source utility functions
source "$UTILS_FILE"

# Supported AMRFinder organisms array
readonly SUPPORTED_ORGANISMS=(
    "Acinetobacter_baumannii" "Burkholderia_cepacia" "Burkholderia_mallei"
    "Burkholderia_pseudomallei" "Campylobacter" "Citrobacter_freundii"
    "Corynebacterium_diphtheriae" "Clostridioides_difficile" "Enterobacter_cloacae"
    "Enterobacter_asburiae" "Enterococcus_faecalis" "Enterococcus_faecium"
    "Escherichia" "Klebsiella_oxytoca" "Klebsiella_pneumoniae"
    "Neisseria_gonorrhoeae" "Neisseria_meningitidis" "Pseudomonas_aeruginosa"
    "Salmonella" "Serratia_marcescens" "Staphylococcus_aureus"
    "Staphylococcus_pseudintermedius" "Streptococcus_agalactiae"
    "Streptococcus_pneumoniae" "Streptococcus_pyogenes" "Vibrio_cholerae"
    "Vibrio_vulnificus" "Vibrio_parahaemolyticus"
)

# Cleanup function
cleanup() {
    log "$LOGS_DIR" "$LOG_NAME" "Cleaning up temporary files..."
    rm -rf "$TEMP_DIR"
    rm -f "$GENOME_ASSEMBLY_DIR/temp_input_assembly_AMRFinder.fasta"
    log "$LOGS_DIR" "$LOG_NAME" "Cleanup complete."
}
trap cleanup EXIT

# Create required directories
create_directory "$AMRFINDER_OUTPUT_DIR"
create_directory "$TEMP_DIR"
create_directory "$QS_FILES_DIR"

# Find required SIF files
straincascade_taxonomic_functional_analysis_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_taxonomic_functional_analysis*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Determine organism name from taxonomic classification
get_organism_name() {
    local organism_name=""
    local organism_name_file
    organism_name_file=$(find "$TAXONOMIC_CLASS_DIR" -name '*gtdbtk_organism_name.txt' -print -quit)
    
    if [[ -s "$organism_name_file" ]] && grep -q "species" "${TAXONOMIC_CLASS_DIR}/level_organism_name.txt" 2>/dev/null; then
        organism_name=$(sed -E 's/[gs]__([^;]*).*/\1/' "$organism_name_file" | sed 's/ /_/g')
        
        # Check if organism is supported
        for supported in "${SUPPORTED_ORGANISMS[@]}"; do
            if [[ "$organism_name" == "$supported"* ]]; then
                echo "$supported"
                return 0
            fi
        done
    fi
    echo ""
}

# Get organism name
organism_name=$(get_organism_name)

# Build AMRFinderPlus command
amrfinder_cmd="export TMPDIR=/tmp && mkdir -p \$TMPDIR && \
    amrfinder \
    -n /mnt/input/temp_input_assembly_AMRFinder.fasta \
    -o /mnt/output/${OUTPUT_BASE}.tsv \
    --database $DATABASES_DIR/amrfinderplus_db/latest \
    --threads $THREADS \
    --plus \
    --mutation_all /mnt/output/${MUTATION_BASE}.tsv"

[[ -n "$organism_name" ]] && amrfinder_cmd+=" --organism $organism_name"

# Log AMRFinderPlus execution
log "$LOGS_DIR" "$LOG_NAME" "Running AMRFinderPlus AMR identification for $analysis_assembly_file"
log "$LOGS_DIR" "$LOG_NAME" "Organism name: ${organism_name:-'Not specified'}"

# Run AMRFinderPlus
apptainer exec \
    --bind "$TEMP_DIR":/mnt/input \
    --bind "$TEMP_DIR":/tmp \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input2 \
    --bind "$AMRFINDER_OUTPUT_DIR":/mnt/output \
    --bind "$DATABASES_DIR/amrfinderplus_db":/opt/conda/envs/amrfinderplus_env/share/amrfinderplus/data \
    "$straincascade_taxonomic_functional_analysis_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate amrfinderplus_env && \
             cp /mnt/input2/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_AMRFinder.fasta && \
             $amrfinder_cmd" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: AMRFinderPlus analysis failed"
    exit 1
}

# Copy output files to functional analysis directory
if ! find "$AMRFINDER_OUTPUT_DIR" -name "${SAMPLE_NAME}_*.tsv" -exec cp {} "$FUNCTIONAL_ANALYSIS_DIR" \;; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No output files found to copy"
    exit 1
fi

# Validate required files for R processing
tsv_file1=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name '*identification_amrfinderplus.tsv' -print -quit)
tsv_file2=$(find "$FUNCTIONAL_ANALYSIS_DIR" -name '*mutations_amrfinderplus.tsv' -print -quit)

for file in "$tsv_file1" "$tsv_file2"; do
    if [[ ! -f "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Required file $file not found"
        exit 1
    fi
done

# Log R processing configuration
log "$LOGS_DIR" "$LOG_NAME" "Processing AMRFinderPlus results with R"

# Process results with R
apptainer exec \
    --bind "$SCRIPT_DIR/R_scripts":/mnt/r_script_dir \
    --bind "$FUNCTIONAL_ANALYSIS_DIR":/mnt/input \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input_fasta \
    --bind "$QS_FILES_DIR":/mnt/output \
    "$r_sif" \
    Rscript "/mnt/r_script_dir/R_process_amrfinderplus.R" \
    --output_dir "/mnt/output" \
    --tsv_file1 "/mnt/input/$(basename "$tsv_file1")" \
    --tsv_file2 "/mnt/input/$(basename "$tsv_file2")" \
    --fasta "/mnt/input_fasta/$(basename "$analysis_assembly_file")" \
    --version "$VERSION" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: R processing failed"
    exit 1
}

log "$LOGS_DIR" "$LOG_NAME" "AMRFinderPlus analysis completed successfully"