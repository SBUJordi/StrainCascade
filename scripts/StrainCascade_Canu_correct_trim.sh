#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_Canu_correct_trim.sh
# Description: Performs Canu read correction and trimming as part of StrainCascade pipeline

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> 
          <input_file> <output_dir> <sample_name> <sequencing_type> 
           <threads> <genome_assembly_main_abs> <reproducibility_mode>
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
readonly INPUT_FILE="$6"
readonly OUTPUT_DIR="$7"
readonly SAMPLE_NAME="$8"
readonly SEQUENCING_TYPE="$9"
readonly INITIAL_THREADS="${10}"
readonly GENOME_ASSEMBLY_DIR="${11}"
readonly REPRODUCIBILITY_MODE="${12}"

# Set threads based on algorithm type
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    readonly THREADS=1
else
    readonly THREADS="${INITIAL_THREADS}"
fi

# Derived constants
readonly CANU_OUTPUT_DIR="$OUTPUT_DIR/Canu_correct_trim_results"
readonly GENOME_SIZE_FILE="$GENOME_ASSEMBLY_DIR/informed_genome_size_estimation.txt"
readonly DEFAULT_GENOME_SIZE="4.5m"

# Source utility functions
source "$UTILS_FILE"

log "$LOGS_DIR" "$LOG_NAME" "Threads set to ${THREADS} based on algorithm type ${REPRODUCIBILITY_MODE}."

# Create output directory
create_directory "$CANU_OUTPUT_DIR"

# Find required Apptainer image
readonly canu_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_assembly*.sif')

# Determine sequencing configuration
readonly SEQ_CONFIG=$(case "$SEQUENCING_TYPE" in
    pacbio-raw)
        echo "-pacbio true true"
        ;;
    pacbio-corr)
        echo "-pacbio -corrected false true"
        ;;
    pacbio-hifi)
        echo "-pacbio-hifi false false"
        ;;
    nano-raw)
        echo "-nanopore true true"
        ;;
    nano-corr|nano-hq)
        echo "-nanopore -corrected false true"
        ;;
    *)
        log "$LOGS_DIR" "$LOG_NAME" "Error: Invalid sequencing type: $SEQUENCING_TYPE"
        log "$LOGS_DIR" "$LOG_NAME" "Valid types: pacbio-raw|pacbio-corr|pacbio-hifi|nano-raw|nano-corr|nano-hq"
        log "$LOGS_DIR" "$LOG_NAME" "Skipping Canu correction and trimming."
        exit 0
        ;;
esac)

# Parse sequencing configuration
read -r -a CONFIG_ARRAY <<< "$SEQ_CONFIG"
if [[ "${CONFIG_ARRAY[1]}" == "-corrected" ]]; then
    CANU_SEQ_TYPE="${CONFIG_ARRAY[0]} ${CONFIG_ARRAY[1]}"
    DO_CORRECTION="${CONFIG_ARRAY[2]}"
    DO_TRIMMING="${CONFIG_ARRAY[3]}"
else
    CANU_SEQ_TYPE="${CONFIG_ARRAY[0]}"
    DO_CORRECTION="${CONFIG_ARRAY[1]}"
    DO_TRIMMING="${CONFIG_ARRAY[2]}"
fi

# Early exit for HiFi data
if [[ "$DO_CORRECTION" == "false" && "$DO_TRIMMING" == "false" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Note: PacBio HiFi data is pre-corrected and trimmed. Skipping processing."
    exit 0
fi

# Determine genome size
readonly GENOME_SIZE=$([[ -f "$GENOME_SIZE_FILE" ]] && 
    awk '{printf "%.1fm\n", $1/1000000}' "$GENOME_SIZE_FILE" || 
    echo "$DEFAULT_GENOME_SIZE") || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Could not determine genome size. Using default."
        echo "$DEFAULT_GENOME_SIZE"
    }

# Log processing parameters
log "$LOGS_DIR" "$LOG_NAME" "Starting Canu processing for $INPUT_FILE"
log "$LOGS_DIR" "$LOG_NAME" "Sequencing type: $SEQUENCING_TYPE"
log "$LOGS_DIR" "$LOG_NAME" "Using genome size: $GENOME_SIZE"
log "$LOGS_DIR" "$LOG_NAME" "Using $THREADS threads"

# Initialize current input tracking
CURRENT_INPUT="$INPUT_FILE"
PROCESSED_OUTPUT="" # Track the most recently processed output

# Create deterministic entropy source
readonly ENTROPY_FILE="$CANU_OUTPUT_DIR/deterministic_entropy_file"
if [[ ! -f "$ENTROPY_FILE" ]]; then
    dd if=/dev/zero bs=1024 count=100 > "$ENTROPY_FILE"
    log "$LOGS_DIR" "$LOG_NAME" "Deterministic entropy file created at $ENTROPY_FILE"
fi

# Correction step if needed
if [[ "$DO_CORRECTION" == "true" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Running correction step"
    
    apptainer exec \
        --bind "$(dirname "$CURRENT_INPUT")":/mnt/input \
        --bind "$CANU_OUTPUT_DIR":/mnt/output \
        --bind "$ENTROPY_FILE":/dev/random \
        --bind "$ENTROPY_FILE":/dev/urandom \
        "$canu_sif" canu \
        -correct \
        -p "${SAMPLE_NAME}_correction" \
        -d /mnt/output \
        genomeSize="$GENOME_SIZE" \
        $CANU_SEQ_TYPE "/mnt/input/$(basename "$CURRENT_INPUT")" \
        useGrid=false \
        cnsThreads="$THREADS" \
        corThreads="$THREADS" \
        redThreads="$THREADS" \
        oeaThreads="$THREADS" \
        batThreads="$THREADS" \
        ovlThreads="$THREADS" \
        mhapThreads="$THREADS" \
        mmapThreads="$THREADS" \
        ovbThreads="$THREADS" \
        ovsThreads="$THREADS" \
        cnsConcurrency="$THREADS" \
        corConcurrency="$THREADS" \
        redConcurrency="$THREADS" \
        oeaConcurrency="$THREADS" \
        batConcurrency="$THREADS" \
        ovlConcurrency="$THREADS" \
        mhapConcurrency="$THREADS" \
        mmapConcurrency="$THREADS" \
        ovbConcurrency="$THREADS" \
        ovsConcurrency="$THREADS" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Canu correction step failed."
        if [[ "$DO_TRIMMING" == "false" ]]; then
            exit 0
        fi
    }
    
    # Process correction output
    readonly CORRECTED_OUTPUT="$(find "$CANU_OUTPUT_DIR" -name "${SAMPLE_NAME}_correction.correctedReads.fasta.gz")"
    if [[ ! -f "$CORRECTED_OUTPUT" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Correction output file not found, proceeding with original input"
    else
        gunzip -c "$CORRECTED_OUTPUT" > "$CANU_OUTPUT_DIR/${SAMPLE_NAME}_corrected.fasta"
        log "$LOGS_DIR" "$LOG_NAME" "Created ${SAMPLE_NAME}_corrected.fasta"
        
        CURRENT_INPUT="$CANU_OUTPUT_DIR/${SAMPLE_NAME}_corrected.fasta"
        PROCESSED_OUTPUT="$CURRENT_INPUT"
        CANU_SEQ_TYPE="$CANU_SEQ_TYPE -corrected"
    fi
fi

# Trimming step if needed
if [[ "$DO_TRIMMING" == "true" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Running trimming step"
    
    apptainer exec \
        --bind "$(dirname "$CURRENT_INPUT")":/mnt/input \
        --bind "$CANU_OUTPUT_DIR":/mnt/output \
        --bind "$ENTROPY_FILE":/dev/random \
        --bind "$ENTROPY_FILE":/dev/urandom \
        "$canu_sif" canu \
        -trim \
        -p "${SAMPLE_NAME}_trimming" \
        -d /mnt/output \
        genomeSize="$GENOME_SIZE" \
        $CANU_SEQ_TYPE "/mnt/input/$(basename "$CURRENT_INPUT")" \
        useGrid=false \
        cnsThreads="$THREADS" \
        corThreads="$THREADS" \
        redThreads="$THREADS" \
        oeaThreads="$THREADS" \
        batThreads="$THREADS" \
        ovlThreads="$THREADS" \
        mhapThreads="$THREADS" \
        mmapThreads="$THREADS" \
        ovbThreads="$THREADS" \
        ovsThreads="$THREADS" \
        cnsConcurrency="$THREADS" \
        corConcurrency="$THREADS" \
        redConcurrency="$THREADS" \
        oeaConcurrency="$THREADS" \
        batConcurrency="$THREADS" \
        ovlConcurrency="$THREADS" \
        mhapConcurrency="$THREADS" \
        mmapConcurrency="$THREADS" \
        ovbConcurrency="$THREADS" \
        ovsConcurrency="$THREADS" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Canu trimming step failed."
        if [[ -n "$PROCESSED_OUTPUT" ]]; then
            log "$LOGS_DIR" "$LOG_NAME" "Using corrected output as final result"
            mv "$PROCESSED_OUTPUT" "$INPUT_FILE"
        fi
        exit 0
    }
    
    # Process trimming output
    readonly TRIMMED_OUTPUT="$(find "$CANU_OUTPUT_DIR" -name "${SAMPLE_NAME}_trimming.trimmedReads.fasta.gz")"
    if [[ ! -f "$TRIMMED_OUTPUT" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Trimming output file not found."
        if [[ -n "$PROCESSED_OUTPUT" ]]; then
            log "$LOGS_DIR" "$LOG_NAME" "Using corrected output as final result"
            mv "$PROCESSED_OUTPUT" "$INPUT_FILE"
        fi
        exit 0
    fi
    
    gunzip -c "$TRIMMED_OUTPUT" > "$CANU_OUTPUT_DIR/${SAMPLE_NAME}_processed.fasta"
    log "$LOGS_DIR" "$LOG_NAME" "Created ${SAMPLE_NAME}_processed.fasta"
    PROCESSED_OUTPUT="$CANU_OUTPUT_DIR/${SAMPLE_NAME}_processed.fasta"
fi

# Final output handling
if [[ -n "$PROCESSED_OUTPUT" ]]; then
    # Validate FASTA format using samtools in container
    if apptainer exec \
        "$canu_sif" samtools faidx "$PROCESSED_OUTPUT" 2>/dev/null; then
        mv "$PROCESSED_OUTPUT" "$INPUT_FILE"
        log "$LOGS_DIR" "$LOG_NAME" "Replaced $(basename "$INPUT_FILE") with processed reads"
        log "$LOGS_DIR" "$LOG_NAME" "Canu processing completed successfully"
    else
        log "$LOGS_DIR" "$LOG_NAME" "Error: Invalid FASTA format detected by samtools in $PROCESSED_OUTPUT - will not be used"
        exit 0
    fi
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: No processing was completed successfully"
    exit 0
fi