#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_minimap2_BBMap_coverage.sh
# Description: Performs read mapping and coverage analysis using minimap2 (long reads) or BBMap (short reads)

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <input_file> 
          <output_dir> <sample_name> <sequencing_type> <threads> <genome_assembly_main_abs>
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
readonly INPUT_FILE="$6"
readonly OUTPUT_DIR="$7"
readonly SAMPLE_NAME="$8"
readonly SEQUENCING_TYPE="$9"
readonly THREADS="${10}"
readonly GENOME_ASSEMBLY_DIR="${11}"

# Derived constants
readonly COVERAGE_OUTPUT_DIR="$OUTPUT_DIR/minimap2_BBMap_coverage_results"
readonly READ_LENGTH_THRESHOLD=600
readonly READ_LENGTH_BUFFER=100

# Source utility functions
source "$UTILS_FILE"

# Validate sequencing type and set minimap2 preset
case "$SEQUENCING_TYPE" in
    *"pacbio-hifi"*) readonly MINIMAP2_PRESET="map-hifi" ;;
    *"pacbio"*)      readonly MINIMAP2_PRESET="map-pb" ;;
    *"nano"*)        readonly MINIMAP2_PRESET="map-ont" ;;
    *)
        log "$LOGS_DIR" "$LOG_NAME" "Skipping: Unsupported sequencing type: $SEQUENCING_TYPE. Requires 'pacbio' or 'nanopore' data."
        exit 0
        ;;
esac

# Create output directory
create_directory "$COVERAGE_OUTPUT_DIR"

# Find analysis assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
[[ -z "$analysis_assembly_file" ]] && {
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found. Skipping coverage analysis."
    exit 0
}

# Validate input file
[[ ! -f "$INPUT_FILE" ]] && {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Input file $INPUT_FILE does not exist."
    exit 0
}

# Calculate maximum read length (sample first 1000 reads for efficiency)
log "$LOGS_DIR" "$LOG_NAME" "Calculating maximum read length from $INPUT_FILE (sampling first 1000 reads, assuming FASTA format)"
max_read_length=$(awk '/^>/ {n++; if (n>1000) exit; if (seqlen>max) max=seqlen; seqlen=0; next} {seqlen+=length($0)} END {if (seqlen>max) max=seqlen; print max}' "$INPUT_FILE") || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to calculate maximum read length from file."
    exit 1
}

# Check if max_read_length was calculated correctly
if [[ -z "$max_read_length" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Maximum read length is empty."
else
    log "$LOGS_DIR" "$LOG_NAME" "Maximum read length: $max_read_length"

    # Add buffer to max read length
    max_read_length=$((max_read_length + READ_LENGTH_BUFFER))
    log "$LOGS_DIR" "$LOG_NAME" "Maximum read length (with ${READ_LENGTH_BUFFER}bp buffer): $max_read_length"
fi

# Find container image
straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')

# Prepare mapping command based on read length
if [[ "$max_read_length" -le $READ_LENGTH_THRESHOLD ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Running BBMap alignment (short reads detected)"
    USE_BBMAP=true
    mapping_cmd="bbmap.sh \
        in=/mnt/sequencing_file_dir/$(basename "$INPUT_FILE") \
        ref=/mnt/assembly_file_dir/$(basename "$analysis_assembly_file") \
        out=/mnt/output/${SAMPLE_NAME}_mapped.sam \
        fastareadlen=$max_read_length \
        overwrite=true \
        path=/mnt/output"
else
    log "$LOGS_DIR" "$LOG_NAME" "Running minimap2 alignment with preset: $MINIMAP2_PRESET"
    USE_BBMAP=false
fi

# Run mapping
if [[ "$USE_BBMAP" == "true" ]]; then
    apptainer exec \
        --bind "$(dirname "$INPUT_FILE")":/mnt/sequencing_file_dir \
        --bind "$(dirname "$analysis_assembly_file")":/mnt/assembly_file_dir \
        --bind "$COVERAGE_OUTPUT_DIR":/mnt/output \
        "$straincascade_assembly_qc_refinement_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate bbmap_env && \
                 $mapping_cmd" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: BBMap mapping failed"
        exit 0
    }
else
    # Use minimap2 for long-read mapping
    apptainer exec \
        --bind "$(dirname "$INPUT_FILE")":/mnt/sequencing_file_dir \
        --bind "$(dirname "$analysis_assembly_file")":/mnt/assembly_file_dir \
        --bind "$COVERAGE_OUTPUT_DIR":/mnt/output \
        "$straincascade_assembly_qc_refinement_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate racon_env && \
                 minimap2 -ax $MINIMAP2_PRESET -t $THREADS \
                     /mnt/assembly_file_dir/$(basename "$analysis_assembly_file") \
                     /mnt/sequencing_file_dir/$(basename "$INPUT_FILE") \
                     > /mnt/output/${SAMPLE_NAME}_mapped.sam" || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: minimap2 mapping failed"
        exit 0
    }
fi

# Generate coverage statistics
log "$LOGS_DIR" "$LOG_NAME" "Generating coverage statistics"
apptainer exec \
    --bind "$COVERAGE_OUTPUT_DIR":/mnt/input_output \
    "$straincascade_assembly_qc_refinement_sif" \
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             conda activate bbmap_env && \
             pileup.sh \
             in=/mnt/input_output/${SAMPLE_NAME}_mapped.sam \
             out=/mnt/input_output/$(basename "${analysis_assembly_file%.*}_coverage.txt") \
             overwrite=true" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Coverage statistics generation failed"
    exit 1
}

# Copy coverage file to assembly directory
coverage_file=$(find "$COVERAGE_OUTPUT_DIR" -type f -name "*_coverage.txt" -print -quit)
if [[ -n "$coverage_file" ]]; then
    cp "$coverage_file" "$GENOME_ASSEMBLY_DIR"
    log "$LOGS_DIR" "$LOG_NAME" "Coverage statistics copied to assembly directory"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Coverage statistics file not found"
    exit 1
fi

# Cleanup
rm -f "$GENOME_ASSEMBLY_DIR"/*.mmi "$GENOME_ASSEMBLY_DIR"/*.fai
log "$LOGS_DIR" "$LOG_NAME" "Temporary index files cleaned up"