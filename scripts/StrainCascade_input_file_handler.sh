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
          <sequencing_reads_main_abs> <genome_assembly_main_abs> <input_type> <short_reads_r1> <short_reads_r2> <threads>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 12 ]] || show_usage

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
readonly SHORT_READS_R1="${10}"
readonly SHORT_READS_R2="${11}"
readonly THREADS="${12}"

# Source utility functions
readonly UTILS_FILE="$SCRIPT_DIR/utils.sh"
[[ -f "$UTILS_FILE" ]] || {
    echo "Error: utils.sh not found in $SCRIPT_DIR"
    exit 1
}
source "$UTILS_FILE"

# Find appropriate Apptainer image
log "$LOGS_DIR" "$LOG_NAME" "Finding Apptainer image for processing"
straincascade_genome_assembly=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_genome_assembly*.sif')

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
            target="$GENOME_ASSEMBLY_DIR/${filename}.fasta"
            if [[ "$(realpath "$processed_input")" != "$(realpath "$target" 2>/dev/null)" ]]; then
                cp "$processed_input" "$target"
            fi
            processed_input="$target"
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
            target="$SEQUENCING_READS_DIR/${filename}.fasta"
            if [[ "$(realpath "$processed_input")" != "$(realpath "$target" 2>/dev/null)" ]]; then
                cp "$processed_input" "$target"
            fi
            processed_input="$target"
            ;;
        fastq|fastq.gz|bam)
            log "$LOGS_DIR" "$LOG_NAME" "Processing $extension reads: $base"
            
            # Handle BAM files
            if [[ "$extension" == "bam" ]]; then
                bam_file="$processed_input"
            fi
            if [[ "$(realpath "$processed_input")" != "$(realpath "$SEQUENCING_READS_DIR/$base" 2>/dev/null)" ]]; then
                cp "$processed_input" "$SEQUENCING_READS_DIR/$base"
            fi
            
            # Convert to FASTA using Apptainer
            log "$LOGS_DIR" "$LOG_NAME" "Converting to FASTA format using $THREADS threads"
            
            # Create local temp directory within output to avoid /scratch/local issues
            local_tmpdir="$SEQUENCING_READS_DIR/tmp_conversion"
            mkdir -p "$local_tmpdir"
            
            apptainer exec \
                --bind "$SEQUENCING_READS_DIR":/mnt/input \
                --bind "$local_tmpdir":/tmp \
                "$straincascade_genome_assembly" \
                bash -c "export TMPDIR=/tmp && \
                        source /opt/conda/etc/profile.d/conda.sh && \
                        conda activate tools_env && \
                        samtools fasta -@ $THREADS /mnt/input/$base > /mnt/input/${filename}.fasta" || {
                log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to convert $extension to FASTA format"
                rm -rf "$local_tmpdir"
                exit 0
            }
            
            # Cleanup temp directory
            rm -rf "$local_tmpdir"
            
            processed_input="$SEQUENCING_READS_DIR/${filename}.fasta"
            
            # Verify output file is non-empty
            if [[ ! -s "$processed_input" ]]; then
                log "$LOGS_DIR" "$LOG_NAME" "Error: Conversion produced empty FASTA file"
                exit 0
            fi
            
            # Extract BAM metadata if applicable
            if [[ "$extension" == "bam" ]]; then
                log "$LOGS_DIR" "$LOG_NAME" "Extracting BAM metadata"
                local_tmpdir="$SEQUENCING_READS_DIR/tmp_bam_meta"
                mkdir -p "$local_tmpdir"
                apptainer exec \
                    --bind "$SEQUENCING_READS_DIR":/mnt/input \
                    --bind "$local_tmpdir":/tmp \
                    "$straincascade_genome_assembly" \
                    bash -c "export TMPDIR=/tmp && \
                            source /opt/conda/etc/profile.d/conda.sh && \
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
                rm -rf "$local_tmpdir"
            fi
            ;;
        *)
            log "$LOGS_DIR" "$LOG_NAME" "Error: Unsupported file extension"
            exit 1
            ;;
    esac
fi

# Process short read files if provided for hybrid assembly
processed_short_r1="not_available"
processed_short_r2="not_available"

if [[ "$SHORT_READS_R1" != "not_provided" && "$SHORT_READS_R2" != "not_provided" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Processing short read files for hybrid assembly"
    
    # Validate and process R1
    if [[ ! -f "$SHORT_READS_R1" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Short reads R1 file not found: $SHORT_READS_R1"
        echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"
        exit 1
    fi
    
    if [[ ! -r "$SHORT_READS_R1" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Short reads R1 file not readable: $SHORT_READS_R1"
        echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"
        exit 1
    fi
    
    # Extract R1 file information
    sr1_base=$(basename "$SHORT_READS_R1")
    sr1_filename=$(basename "$sr1_base" .gz)
    sr1_extension="${sr1_filename##*.}"
    
    # Validate R1 extension
    case "$sr1_extension" in
        fastq|fq|fastq.gz|fq.gz)
            log "$LOGS_DIR" "$LOG_NAME" "Processing short reads R1: $sr1_base (format: $sr1_extension)"
            cp "$SHORT_READS_R1" "$SEQUENCING_READS_DIR/$sr1_base"
            processed_short_r1="$SEQUENCING_READS_DIR/$sr1_base"
            ;;
        *)
            log "$LOGS_DIR" "$LOG_NAME" "Error: Unsupported short reads R1 extension. Expected .fastq, .fq, .fastq.gz, or .fq.gz"
            echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"
            exit 1
            ;;
    esac
    
    # Validate and process R2
    if [[ ! -f "$SHORT_READS_R2" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Short reads R2 file not found: $SHORT_READS_R2"
        echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"
        exit 1
    fi
    
    if [[ ! -r "$SHORT_READS_R2" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Short reads R2 file not readable: $SHORT_READS_R2"
        echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"
        exit 1
    fi
    
    # Extract R2 file information
    sr2_base=$(basename "$SHORT_READS_R2")
    sr2_filename=$(basename "$sr2_base" .gz)
    sr2_extension="${sr2_filename##*.}"
    
    # Validate R2 extension
    case "$sr2_extension" in
        fastq|fq|fastq.gz|fq.gz)
            log "$LOGS_DIR" "$LOG_NAME" "Processing short reads R2: $sr2_base (format: $sr2_extension)"
            cp "$SHORT_READS_R2" "$SEQUENCING_READS_DIR/$sr2_base"
            processed_short_r2="$SEQUENCING_READS_DIR/$sr2_base"
            ;;
        *)
            log "$LOGS_DIR" "$LOG_NAME" "Error: Unsupported short reads R2 extension. Expected .fastq, .fq, .fastq.gz, or .fq.gz"
            echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"
            exit 1
            ;;
    esac
    
    log "$LOGS_DIR" "$LOG_NAME" "Short read files processed successfully"
fi

# Output results
echo -e "$processed_input\t$bam_file\t$processed_short_r1\t$processed_short_r2"