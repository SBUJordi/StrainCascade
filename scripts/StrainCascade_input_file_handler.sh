#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_input_file_handler.sh
# Description: Handles input file processing and format conversion for StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <apptainer_images_dir> <input_file> <output_dir> 
          <sequencing_reads_main_abs> <genome_assembly_main_abs> <input_type>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 9 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly APPTAINER_DIR="$4"
readonly INPUT_FILE="$5"
readonly OUTPUT_DIR="$6"
readonly SEQUENCING_READS_DIR="$7"
readonly GENOME_ASSEMBLY_DIR="$8"
readonly INPUT_TYPE="$9"

# Source utility functions
readonly UTILS_FILE="$SCRIPT_DIR/utils.sh"
[[ -f "$UTILS_FILE" ]] || {
    echo "Error: utils.sh not found in $SCRIPT_DIR"
    exit 1
}
source "$UTILS_FILE"

# Find appropriate Apptainer image
log "$LOGS_DIR" "$LOG_NAME" "Finding Apptainer image for processing"
straincascade_genome_assembly=$(find "$APPTAINER_DIR" -name 'straincascade_genome_assembly*.sif' -print -quit)

[[ -f "$straincascade_genome_assembly" ]] || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: No matching .sif file found in $APPTAINER_DIR"
    exit 1
}

# Initialize working variables
processed_input="$INPUT_FILE"
bam_file="not_available"

# Process input file
[[ -n "$processed_input" ]] || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: No input file provided"
    exit 1
}

# Extract file information
dir=$(dirname "$processed_input")
base=$(basename "$processed_input")
filename=$(basename "$base" .gz)
extension="${filename##*.}"
filename="${filename%.*}"

# Handle double extensions (e.g., fastq.gz)
[[ "$base" == *.fastq.gz ]] && extension="fastq.gz"

# Process assembly input
if [[ "$INPUT_TYPE" == "assembly" ]]; then
    case "$extension" in
        fasta|fa|fna)
            log "$LOGS_DIR" "$LOG_NAME" "Processing assembly file: $base"
            cp "$processed_input" "$GENOME_ASSEMBLY_DIR/${filename}.fasta"
            processed_input="$GENOME_ASSEMBLY_DIR/${filename}.fasta"
            ;;
        *)
            log "$LOGS_DIR" "$LOG_NAME" "Error: Unsupported assembly extension. Expected .fasta, .fa, or .fna"
            exit 1
            ;;
    esac
else
    # Process reads input
    case "$extension" in
        fasta|fa|fna)
            log "$LOGS_DIR" "$LOG_NAME" "Processing FASTA reads: $base"
            cp "$processed_input" "$SEQUENCING_READS_DIR/${filename}.fasta"
            processed_input="$SEQUENCING_READS_DIR/${filename}.fasta"
            ;;
        fastq|fastq.gz|bam)
            log "$LOGS_DIR" "$LOG_NAME" "Processing $extension reads: $base"
            
            # Handle BAM files
            if [[ "$extension" == "bam" ]]; then
                bam_file="$processed_input"
                cp "$processed_input" "$SEQUENCING_READS_DIR/$base"
            else
                cp "$processed_input" "$SEQUENCING_READS_DIR/$base"
            fi
            
            # Add deterministic entropy source
            readonly ENTROPY_FILE="$SEQUENCING_READS_DIR/deterministic_entropy_file"
            if [[ ! -f "$ENTROPY_FILE" ]]; then
                dd if=/dev/zero bs=1024 count=100 > "$ENTROPY_FILE"
                log "$LOGS_DIR" "$LOG_NAME" "Deterministic entropy file created"
            fi

            # Convert to FASTA using Apptainer
            log "$LOGS_DIR" "$LOG_NAME" "Converting to FASTA format"
            apptainer exec \
                --bind "$SEQUENCING_READS_DIR":/mnt/input \
                --bind "$ENTROPY_FILE":/dev/random \
                --bind "$ENTROPY_FILE":/dev/urandom \
                "$straincascade_genome_assembly" \
                bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                        conda activate tools_env && \
                        samtools fasta /mnt/input/$base > /mnt/input/${filename}.fasta"
            
            processed_input="$SEQUENCING_READS_DIR/${filename}.fasta"
            
            # Extract BAM metadata if applicable
            if [[ "$extension" == "bam" ]]; then
                log "$LOGS_DIR" "$LOG_NAME" "Extracting BAM metadata"
                apptainer exec \
                    --bind "$SEQUENCING_READS_DIR":/mnt/input \
                    --bind "$ENTROPY_FILE":/dev/random \
                    --bind "$ENTROPY_FILE":/dev/urandom \
                    "$straincascade_genome_assembly" \
                    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate tools_env && \
                            samtools view -H /mnt/input/$base | awk '
                            BEGIN {
                                print \"BAM File: '$bam_file'\"
                                print \"--------------------\"
                            }
                            {
                                if (\$1 ~ /^@(RG|SQ|PG|HD)/) {
                                    for (i=1; i<=NF; i++) {
                                        tag = substr(\$i,1,3)
                                        value = substr(\$i, 4)
                                        if (tag ~ /^(ID|SM|LB|PL|PU|SN|LN|PN|CL|VN|SO|DS):/) {
                                            print tag \" \" value
                                        }
                                    }
                                }
                            }' > /mnt/input/${filename}_bam_info.txt || \
                            echo \"Error: Failed to extract BAM metadata\" > /mnt/input/${filename}_bam_info.txt"
            fi
            ;;
        *)
            log "$LOGS_DIR" "$LOG_NAME" "Error: Unsupported file extension"
            exit 1
            ;;
    esac
fi

# Output results
echo -e "$processed_input\t$bam_file"