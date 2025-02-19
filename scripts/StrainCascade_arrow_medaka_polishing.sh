#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_arrow_medaka_polishing.sh
# Description: Performs Arrow (PacBio) or Medaka (Nanopore) polishing as part of StrainCascade

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <input_file>
          <bam_file> <output_dir> <sample_name> <sequencing_type> <threads>
          <genome_assembly_dir> <reproducibility_mode>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 13 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly INPUT_FILE="$6"
readonly BAM_FILE="$7"
readonly OUTPUT_DIR="$8"
readonly SAMPLE_NAME="$9"
readonly SEQUENCING_TYPE="${10}"
readonly INITIAL_THREADS="${11}"
readonly GENOME_ASSEMBLY_DIR="${12}"
readonly REPRODUCIBILITY_MODE="${13}"

# Set threads based on algorithm type
if [[ "${REPRODUCIBILITY_MODE}" == "deterministic" ]]; then
    readonly THREADS=1
else
    readonly THREADS="${INITIAL_THREADS}"
fi

# Derived constants
readonly POLISHING_ROUNDS=2
readonly POLISHING_TOOL=$(
    if [[ "$SEQUENCING_TYPE" == *"pacbio"* && -n "$BAM_FILE" && "$BAM_FILE" != "not_available" ]]; then
        echo "arrow"
    elif [[ "$SEQUENCING_TYPE" == *"nano"* ]]; then
        echo "medaka"
    else
        echo "no"
    fi
)
readonly POLISHING_OUTPUT_DIR="$OUTPUT_DIR/${POLISHING_TOOL}_polishing_results"
readonly FIXED_BAM_FILE="$POLISHING_OUTPUT_DIR/fixed_$(basename "$BAM_FILE")"

# Source utility functions
source "$UTILS_FILE"

log "$LOGS_DIR" "$LOG_NAME" "Threads set to ${THREADS} based on algorithm type ${REPRODUCIBILITY_MODE}."

# Create required directories
create_directory "$POLISHING_OUTPUT_DIR"

# Find required SIF file and assembly file
straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")

# Validate prerequisites
if [[ -z "$POLISHING_TOOL" || "$POLISHING_TOOL" == "no" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No valid sequencing type or necessary BAM file identified. Skipping polishing."
    exit 0
fi

if [[ -z "$analysis_assembly_file" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: No assembly files found. Skipping polishing."
    exit 0
fi

for file in "$INPUT_FILE" $([ "$POLISHING_TOOL" = "arrow" ] && echo "$BAM_FILE"); do
    if [[ -n "$file" && ! -f "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Required file $file not found. Skipping polishing."
        exit 0
    fi
done

# Create deterministic entropy source
readonly ENTROPY_FILE="$POLISHING_OUTPUT_DIR/deterministic_entropy_file"
if [[ ! -f "$ENTROPY_FILE" ]]; then
    dd if=/dev/zero bs=1024 count=100 > "$ENTROPY_FILE"
    log "$LOGS_DIR" "$LOG_NAME" "Deterministic entropy file created at $ENTROPY_FILE"
fi

log "$LOGS_DIR" "$LOG_NAME" "Starting $POLISHING_TOOL polishing with $POLISHING_ROUNDS rounds"

if [[ "$POLISHING_TOOL" == "arrow" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Processing BAM file for Arrow polishing"
    
    # Pre-calculate basename to avoid escaping issues
    bam_basename=$(basename "$BAM_FILE")
    
    # Run the standardization process inside the Apptainer
    output=$(apptainer exec \
        --bind "$(dirname "$BAM_FILE")":/mnt/input \
        --bind "$POLISHING_OUTPUT_DIR":/mnt/output \
        --bind "$ENTROPY_FILE":/dev/random \
        --bind "$ENTROPY_FILE":/dev/urandom \
        "$straincascade_assembly_qc_refinement_sif" \
        /bin/bash -c "
        new_rg_id=\$(printf '%08x' \$RANDOM) && \
        source /opt/conda/etc/profile.d/conda.sh && \
        conda activate tools_env && \
        samtools addreplacerg \
        -r \"ID:\${new_rg_id}\" \
        -r \"PL:PACBIO\" \
        -r \"SM:$SAMPLE_NAME\" \
        -o '/mnt/output/fixed_${bam_basename}' \
        '/mnt/input/${bam_basename}' && \
        samtools index '/mnt/output/fixed_${bam_basename}' && \
        samtools view -H '/mnt/output/fixed_${bam_basename}'
        " 2>&1)

    # Check if the output contains error messages about PacBio BAM compatibility
    if echo "$output" | grep -q "This does not appear to be a valid PacBio BAM file"; then
        log "$LOGS_DIR" "$LOG_NAME" "The BAM file is not compatible with Arrow polishing (not a valid PacBio BAM file). Skipping this module."
        exit 0
    elif [ $? -ne 0 ]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to standardize RG IDs and add necessary tags in the BAM file. Skipping this module."
        exit 0
    fi
fi

# Perform polishing rounds
current_assembly="$analysis_assembly_file"
for ((round=1; round<=POLISHING_ROUNDS; round++)); do
    round_dir="$POLISHING_OUTPUT_DIR/round_${round}"
    consensus_file="$round_dir/consensus.fasta"
    create_directory "$round_dir"
    
    log "$LOGS_DIR" "$LOG_NAME" "Starting polishing round $round"
    
    # Prepare polishing command based on tool
    if [[ "$POLISHING_TOOL" == "medaka" ]]; then
    
    # First check if bacterial model is compatible
    log "$LOGS_DIR" "$LOG_NAME" "Checking bacterial model compatibility"
    
    model_check=$(apptainer exec \
        --bind "$(dirname "$INPUT_FILE")":/mnt/input_file \
        "$straincascade_assembly_qc_refinement_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                conda activate medaka_env && \
                medaka tools list_models 2>/dev/null | grep -q 'r1041_e82_400bps_bacterial_methylation' && \
                medaka tools check_basecalls -i /mnt/input_file/$(basename "$INPUT_FILE") 2>/dev/null") || {
        log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to check bacterial model compatibility. Using default model."
        bacteria_param=""
    }
    
        # Set bacterial flag if compatible
        bacteria_param=""
        if [[ $? -eq 0 ]]; then
            bacteria_param="--model r1041_e82_400bps_bacterial_methylation"
            log "$LOGS_DIR" "$LOG_NAME" "Using bacterial-specific model"
        else
            log "$LOGS_DIR" "$LOG_NAME" "Using default model (bacterial model not compatible)"
        fi
        
        for batch_size in "" "100" "50"; do
            batch_param=""
            [[ -n "$batch_size" ]] && batch_param="-b $batch_size"
            
            log "$LOGS_DIR" "$LOG_NAME" "Attempting Medaka polishing${batch_size:+ with batch size $batch_size}${bacteria_param:+ using bacterial model}"
            
            if apptainer exec \
                --bind "$(dirname "$INPUT_FILE")":/mnt/input_file \
                --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
                --bind "$round_dir":/mnt/output \
                --bind "$ENTROPY_FILE":/dev/random \
                --bind "$ENTROPY_FILE":/dev/urandom \
                "$straincascade_assembly_qc_refinement_sif" \
                bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                        conda activate medaka_env && \
                        export TF_FORCE_GPU_ALLOW_GROWTH=true && \
                        medaka_consensus \
                            -i /mnt/input_file/$(basename "$INPUT_FILE") \
                            -d /mnt/input_assembly/$(basename "$current_assembly") \
                            -o /mnt/output \
                            $batch_param \
                            $bacteria_param \
                            -t $THREADS && \
                        [[ -f /mnt/output/consensus.fasta ]]"; then
                break
            elif [[ "$batch_size" == "50" ]]; then
                log "$LOGS_DIR" "$LOG_NAME" "Error: Medaka polishing round $round failed with all batch sizes. Skipping this module."
                exit 0
            else
                log "$LOGS_DIR" "$LOG_NAME" "Retrying with smaller batch size"
                rm -f "$round_dir"/* # Clean failed attempt
            fi
        done
    else
        apptainer exec \
            --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
            --bind "$round_dir":/mnt/output \
            --bind "$(dirname "$FIXED_BAM_FILE")":/mnt/bam_file_dir \
            --bind "$ENTROPY_FILE":/dev/random \
            --bind "$ENTROPY_FILE":/dev/urandom \
            "$straincascade_assembly_qc_refinement_sif" \
            bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                    conda activate genomicconsensus_env && \
                    variantCaller --algorithm=arrow \
                        -j $THREADS \
                        -r /mnt/input_assembly/$(basename "$current_assembly") \
                        -o /mnt/output/consensus.fasta \
                        /mnt/bam_file_dir/$(basename "$FIXED_BAM_FILE") && \
                    [[ -f /mnt/output/consensus.fasta ]]" || {
            log "$LOGS_DIR" "$LOG_NAME" "Error: Arrow polishing round $round failed. Skipping this module."
            exit 0
        }
    fi
    
    # Verify final output exists for both tools
    if [[ ! -f "$consensus_file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Arrow Consensus file not produced in round $round. Skipping this module."
        exit 0
    fi
    
    current_assembly="$consensus_file"
    log "$LOGS_DIR" "$LOG_NAME" "Completed polishing round $round"
done

# Update final assembly
if [[ -f "$current_assembly" ]]; then
    cp "$current_assembly" "$analysis_assembly_file"
    log "$LOGS_DIR" "$LOG_NAME" "Successfully updated assembly with polished version"
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Final polished assembly not found"
    exit 1
fi

# Cleanup temporary files
if [[ "$POLISHING_TOOL" == "arrow" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Cleaning up temporary files"
    rm -f "$FIXED_BAM_FILE" "${FIXED_BAM_FILE}.bai" "${FIXED_BAM_FILE}.pbi"
fi

log "$LOGS_DIR" "$LOG_NAME" "$POLISHING_TOOL polishing completed successfully"