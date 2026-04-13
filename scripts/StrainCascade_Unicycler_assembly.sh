#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_Unicycler_assembly.sh
# Description: Performs Unicycler genome assembly (long-read-only or hybrid) as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <input_file> 
          <output_dir> <sample_name> <sequencing_type> <threads> <genome_assembly_main_abs> 
          <reproducibility_mode> <short_reads_r1> <short_reads_r2>
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
readonly INPUT_FILE="$6"
readonly OUTPUT_DIR="$7"
readonly SAMPLE_NAME="$8"
readonly SEQUENCING_TYPE="$9"
readonly INITIAL_THREADS="${10}"
readonly GENOME_ASSEMBLY_DIR="${11}"
readonly REPRODUCIBILITY_MODE="${12}"
readonly SHORT_READS_R1="${13}"
readonly SHORT_READS_R2="${14}"

# Set threads based on reproducibility mode
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    readonly THREADS=1
else
    readonly THREADS="${INITIAL_THREADS}"
fi

# Derived constants
readonly UNICYCLER_OUTPUT_DIR="$OUTPUT_DIR/Unicycler_assembly_results"
readonly MIN_VALID_GENOME_SIZE=1300000
readonly MAX_PLAUSIBLE_GENOME_SIZE=10000000
readonly ASSEMBLY_PREFIX="${SAMPLE_NAME}_assembly_uc"
readonly FINAL_ASSEMBLY_NAME="${ASSEMBLY_PREFIX}.fasta"

# Source utility functions
source "$UTILS_FILE"

log "$LOGS_DIR" "$LOG_NAME" "Threads set to ${THREADS} based on reproducibility mode ${REPRODUCIBILITY_MODE}."

# Clean up existing Unicycler output directory to prevent reuse issues
if [[ -d "$UNICYCLER_OUTPUT_DIR" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Removing existing Unicycler output directory to ensure clean run"
    rm -rf "$UNICYCLER_OUTPUT_DIR"
fi

# Create required output directory
create_directory "$UNICYCLER_OUTPUT_DIR"

# Clean up empty FASTA files
log "$LOGS_DIR" "$LOG_NAME" "Checking for empty FASTA files in $GENOME_ASSEMBLY_DIR"
while IFS= read -r -d '' file; do
    if [[ -f "$file" && ! -s "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Removing empty file: $(basename "$file")"
        rm "$file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

# Find required Apptainer image
readonly straincascade_genome_assembly_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_assembly*.sif')

# Determine assembly mode
if [[ -n "$SHORT_READS_R1" && -n "$SHORT_READS_R2" \
   && "$SHORT_READS_R1" != "not_provided" && "$SHORT_READS_R2" != "not_provided" \
   && "$SHORT_READS_R1" != "not_available" && "$SHORT_READS_R2" != "not_available" ]]; then
    readonly ASSEMBLY_MODE="hybrid"
    log "$LOGS_DIR" "$LOG_NAME" "Running Unicycler in hybrid mode (long reads + short reads)"
    log "$LOGS_DIR" "$LOG_NAME" "Long reads: $(basename "$INPUT_FILE")"
    log "$LOGS_DIR" "$LOG_NAME" "Short reads R1: $(basename "$SHORT_READS_R1")"
    log "$LOGS_DIR" "$LOG_NAME" "Short reads R2: $(basename "$SHORT_READS_R2")"
    
    # Validate short read files exist
    if [[ ! -f "$SHORT_READS_R1" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Short reads R1 file not found: $SHORT_READS_R1. Skipping Unicycler assembly."
        exit 0
    fi
    if [[ ! -f "$SHORT_READS_R2" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Short reads R2 file not found: $SHORT_READS_R2. Skipping Unicycler assembly."
        exit 0
    fi
else
    readonly ASSEMBLY_MODE="long"
    log "$LOGS_DIR" "$LOG_NAME" "Running Unicycler in long-read-only mode"
    log "$LOGS_DIR" "$LOG_NAME" "Long reads: $(basename "$INPUT_FILE")"
fi

# Genome size estimation file location
readonly GENOME_SIZE_FILE="$GENOME_ASSEMBLY_DIR/informed_genome_size_estimation.txt"

# Entropy source binding for deterministic reproducibility only
# Using /dev/zero (infinite stream) instead of a finite file to prevent
# exhaustion-related hangs during long-running assembly stages
ENTROPY_ARGS=()
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    ENTROPY_ARGS=(--bind /dev/zero:/dev/random --bind /dev/zero:/dev/urandom)
    log "$LOGS_DIR" "$LOG_NAME" "Deterministic mode: binding /dev/zero as entropy source"
fi

# Build Unicycler command
unicycler_cmd="unicycler -o /mnt/output -t $THREADS"

# Add deterministic flags if in deterministic mode
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    unicycler_cmd="$unicycler_cmd --no_rotate --no_pilon"
    log "$LOGS_DIR" "$LOG_NAME" "Deterministic mode enabled: using --no_rotate --no_pilon"
fi

# Add input files based on assembly mode
if [[ "$ASSEMBLY_MODE" == "hybrid" ]]; then
    unicycler_cmd="$unicycler_cmd -1 /mnt/short_r1/$(basename "$SHORT_READS_R1") -2 /mnt/short_r2/$(basename "$SHORT_READS_R2") -l /mnt/long/$(basename "$INPUT_FILE")"
else
    unicycler_cmd="$unicycler_cmd -l /mnt/long/$(basename "$INPUT_FILE")"
fi

# Run Unicycler
log "$LOGS_DIR" "$LOG_NAME" "Running Unicycler assembler"
log "$LOGS_DIR" "$LOG_NAME" "Output directory: $UNICYCLER_OUTPUT_DIR"
log "$LOGS_DIR" "$LOG_NAME" "Using $THREADS threads"

if [[ "$ASSEMBLY_MODE" == "hybrid" ]]; then
    apptainer exec \
        --bind "$(dirname "$INPUT_FILE")":/mnt/long \
        --bind "$(dirname "$SHORT_READS_R1")":/mnt/short_r1 \
        --bind "$(dirname "$SHORT_READS_R2")":/mnt/short_r2 \
        --bind "$UNICYCLER_OUTPUT_DIR":/mnt/output \
        ${ENTROPY_ARGS[@]+"${ENTROPY_ARGS[@]}"} \
        "$straincascade_genome_assembly_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate genome_assembly_env && \
                 $unicycler_cmd" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Unicycler hybrid assembly failed for $INPUT_FILE. Skipping Unicycler assembly."
        exit 0
    }
else
    apptainer exec \
        --bind "$(dirname "$INPUT_FILE")":/mnt/long \
        --bind "$UNICYCLER_OUTPUT_DIR":/mnt/output \
        ${ENTROPY_ARGS[@]+"${ENTROPY_ARGS[@]}"} \
        "$straincascade_genome_assembly_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate genome_assembly_env && \
                 $unicycler_cmd" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Unicycler long-read assembly failed for $INPUT_FILE. Skipping Unicycler assembly."
        exit 0
    }
fi

# Process assembly output
if [[ -f "$UNICYCLER_OUTPUT_DIR/assembly.fasta" ]]; then
    # Check if file is empty
    if [[ ! -s "$UNICYCLER_OUTPUT_DIR/assembly.fasta" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Unicycler produced an empty assembly file. Skipping Unicycler assembly."
        exit 0
    fi

    mv "$UNICYCLER_OUTPUT_DIR/assembly.fasta" "$UNICYCLER_OUTPUT_DIR/$FINAL_ASSEMBLY_NAME"
    log "$LOGS_DIR" "$LOG_NAME" "Renamed assembly.fasta to $FINAL_ASSEMBLY_NAME"
    
    cp "$UNICYCLER_OUTPUT_DIR/$FINAL_ASSEMBLY_NAME" "$GENOME_ASSEMBLY_DIR/"
    log "$LOGS_DIR" "$LOG_NAME" "Copied assembly to $GENOME_ASSEMBLY_DIR"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Assembly file not found in $UNICYCLER_OUTPUT_DIR. Skipping Unicycler assembly."
    exit 0
fi

# Perform genome size estimation if needed
if [[ ! -f "$GENOME_SIZE_FILE" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Performing genome size estimation"
    
    # Select first assembly file for estimation
    readarray -t fasta_files < <(find "$GENOME_ASSEMBLY_DIR" -name "*.fasta")
    if [[ ${#fasta_files[@]} -eq 0 ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: No FASTA files found for genome size estimation. Skipping Unicycler assembly."
        exit 0
    fi
    selected_fasta="${fasta_files[0]}"
    
    # Calculate genome size with error handling
    if ! genome_size=$(grep -v ">" "$selected_fasta" | tr -d '\n' | tr -cd 'ACGTNacgtn' | wc -c); then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to calculate genome size from assembly. Skipping Unicycler assembly."
        exit 0
    fi
    log "$LOGS_DIR" "$LOG_NAME" "Estimated genome size: $genome_size bp"
    
    # Validate genome size
    if (( genome_size > MAX_PLAUSIBLE_GENOME_SIZE )); then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Genome size $genome_size bp is implausibly large for gut bacteria \
(https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a). Skipping Unicycler assembly."
        exit 0
    elif (( genome_size <= MIN_VALID_GENOME_SIZE )); then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Genome size $genome_size bp may be too small for independent cell replication \
(https://www.science.org/doi/10.1126/science.1114057, Fig. 1). Skipping Unicycler assembly."
        exit 0
    else
        log "$LOGS_DIR" "$LOG_NAME" "Genome size $genome_size bp is within expected range"
        echo "$genome_size" > "$GENOME_SIZE_FILE"
    fi
fi

# Clean up empty FASTA files
log "$LOGS_DIR" "$LOG_NAME" "Checking for empty FASTA files in $GENOME_ASSEMBLY_DIR"
while IFS= read -r -d '' file; do
    if [[ -f "$file" && ! -s "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Removing empty file: $(basename "$file")"
        rm "$file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

log "$LOGS_DIR" "$LOG_NAME" "Unicycler assembly completed successfully"
