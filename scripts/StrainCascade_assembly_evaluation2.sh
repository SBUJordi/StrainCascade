#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_assembly_evaluation2.sh
# Description: Evaluates bacterial genome assemblies using QUAST and selects the best assembly based on specified selection algorithm as part of StrainCascade (Second evaluation phase)

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_directory> 
          <sample_name> <threads> <genome_assembly_main_abs> <selection_algorithm>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 10 ]] || show_usage

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
readonly SELECTION_ALGORITHM="${10}"

# Derived constants
readonly EVAL_OUTPUT_DIR="$OUTPUT_DIR/StrainCascade_assembly_evaluation"
readonly QUAST_REPORT_FILE="${SAMPLE_NAME}_quast_assembly_evaluation.tsv"

# Source utility functions
source "$UTILS_FILE"

# Find required Apptainer images
straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
python_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'python_3.12.4*.sif')

# Find assembly files, excluding previously processed ones
mapfile -t assemblies < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" \
    ! -name "*_best_ev2*" ! -name "*_best_ev3*" ! -name "*_circularised.fasta")

# Check if any assemblies were found
if [[ ${#assemblies[@]} -eq 0 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found in '$GENOME_ASSEMBLY_DIR'. Skipping reevaluation step."
    exit 0
fi

# Clean up old evaluation results to handle re-runs properly
log "$LOGS_DIR" "$LOG_NAME" "Removing existing evaluation directory for clean re-run"
rm -rf "$EVAL_OUTPUT_DIR"

# Remove old best_ev2 files from genome assembly directory
while IFS= read -r -d '' file; do
    log "$LOGS_DIR" "$LOG_NAME" "Removing old evaluation file: $(basename "$file")"
    rm -f "$file"
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*_best_ev2.fasta" -print0)

create_directory "$EVAL_OUTPUT_DIR"

# Process each assembly with QUAST (parallel execution, each using 1 thread)
quast_output_dirs=()
quast_files=()
running_jobs=0

for assembly in "${assemblies[@]}"; do
    file=$(basename "$assembly")
    pattern_name=$(basename "$file" .fasta)
    quast_output_dir="${EVAL_OUTPUT_DIR}/${SAMPLE_NAME}_${pattern_name}_quast_output"
    
    create_directory "$quast_output_dir"
    log "$LOGS_DIR" "$LOG_NAME" "Launching QUAST analysis for $file"
    
    apptainer exec \
        --bind "$GENOME_ASSEMBLY_DIR:/data" \
        --bind "$quast_output_dir:/output" \
        "$straincascade_assembly_qc_refinement_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate quast_env && \
                 quast.py /data/$file -o /output -t 1" &
    
    quast_output_dirs+=("$quast_output_dir")
    quast_files+=("$file")
    running_jobs=$((running_jobs + 1))
    
    # Limit concurrent QUAST jobs to THREADS
    if (( running_jobs >= THREADS )); then
        wait -n 2>/dev/null || true
        running_jobs=$((running_jobs - 1))
    fi
done

# Wait for all remaining QUAST jobs to complete
wait 2>/dev/null || true

# Collect QUAST reports sequentially
for i in "${!quast_output_dirs[@]}"; do
    quast_output_dir="${quast_output_dirs[$i]}"
    file="${quast_files[$i]}"
    
    # Append QUAST report if available (skip header if report already has content)
    if [[ -f "$quast_output_dir/transposed_report.tsv" ]]; then
        if [[ -s "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE" ]]; then
            tail -n +2 "$quast_output_dir/transposed_report.tsv" >> "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE"
        else
            cat "$quast_output_dir/transposed_report.tsv" >> "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE"
        fi
    else
        log "$LOGS_DIR" "$LOG_NAME" "Warning: QUAST failed to produce output for $file"
    fi
done

# Run assembly selection
log "$LOGS_DIR" "$LOG_NAME" "Running assembly selection algorithm: $SELECTION_ALGORITHM"

best_assembly=$(apptainer exec \
    --bind "$EVAL_OUTPUT_DIR:/data" \
    --bind "$SCRIPT_DIR:/scripts" \
    "$python_sif" \
    python "/scripts/StrainCascade_assembly_selection_${SELECTION_ALGORITHM}.py" \
    "/data/$QUAST_REPORT_FILE" "/data")

# Process best assembly result
if [[ -n "$best_assembly" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Best Assembly (ev2): $best_assembly"
    
    # Find and copy the best assembly
    while IFS= read -r assembly; do
        if [[ -f "$assembly" ]]; then
            prefix=$(basename "${assembly%.*}")
            cp "$assembly" "${EVAL_OUTPUT_DIR}/${prefix}_best_ev2.fasta"
            cp "$assembly" "${GENOME_ASSEMBLY_DIR}/${prefix}_best_ev2.fasta"
            log "$LOGS_DIR" "$LOG_NAME" "Best assembly copied: ${prefix}_best_ev2.fasta"
            echo "Reevaluation (II) complete. Best reevaluated assembly: ${prefix}_best_ev2.fasta"
            exit 0
        fi
    done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*${best_assembly}*")
    
    log "$LOGS_DIR" "$LOG_NAME" "Error: Could not find best assembly file"
    exit 1
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Assembly selection failed to return valid assembly name"
    exit 1
fi