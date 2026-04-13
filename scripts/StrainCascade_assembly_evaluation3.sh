#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_assembly_evaluation3.sh - Version 2.0.0
# Description: Final evaluation of bacterial genome assemblies using QUAST, selects best assembly, and processes results with R as part of StrainCascade (Third evaluation phase)

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> <output_directory> 
          <sample_name> <threads> <genome_assembly_main_abs> <results_integration_abs> 
          <version> <selection_algorithm>
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
readonly RESULTS_INTEGRATION_DIR="${10}"
readonly VERSION="${11}"
readonly SELECTION_ALGORITHM="${12}"

# Derived constants
readonly EVAL_OUTPUT_DIR="$OUTPUT_DIR/StrainCascade_assembly_evaluation"
readonly QUAST_REPORT_FILE="${SAMPLE_NAME}_quast_assembly_evaluation.tsv"
readonly QS_FILES_DIR="$RESULTS_INTEGRATION_DIR/qs_files"
readonly R_SCRIPT_DIR="$SCRIPT_DIR/R_scripts"

# Source utility functions
source "$UTILS_FILE"

# Clean up old evaluation results to handle re-runs properly
if [[ -d "$EVAL_OUTPUT_DIR" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Removing existing evaluation directory for clean re-run"
    rm -rf "$EVAL_OUTPUT_DIR"
fi

# Remove old best_ev3 files from genome assembly directory
while IFS= read -r -d '' file; do
    log "$LOGS_DIR" "$LOG_NAME" "Removing old evaluation file: $(basename "$file")"
    rm -f "$file"
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*_best_ev3.fasta" -print0)

# Find required Apptainer images
straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')
python_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'python_3.12.4*.sif')
r_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'r_4.4.1*.sif')

# Find assembly files, excluding previously processed ones
mapfile -t assemblies < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" ! -name "*_best_ev3*")

# Check if any assemblies were found
if [[ ${#assemblies[@]} -eq 0 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No assembly files found in '$GENOME_ASSEMBLY_DIR'. Skipping reevaluation step."
    exit 0
fi

# Create output directories
create_directory "$EVAL_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

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
    log "$LOGS_DIR" "$LOG_NAME" "Best Assembly (ev3): $best_assembly"
    
    # Find and process the best assembly
    while IFS= read -r assembly; do
        if [[ -f "$assembly" ]]; then
            file=$(basename "$assembly")
            
            # Extract the base assembler name (e.g., "flye" from "SRR123_assembly_flye_best_ev1_best_ev2.fasta")
            # Remove sample prefix, _assembly_ prefix, and all _best_ev* suffixes
            assembler_name=$(echo "$file" | sed -E 's/^[^_]+_assembly_//' | sed -E 's/_best_ev[0-9]+//g' | sed -E 's/_mac2\.0_merged_assembly//' | sed -E 's/_circularised//' | sed 's/\.fasta$//')
            
            # Create clean final assembly name: {sample}_assembly_{assembler}_final.fasta
            readonly CLEAN_FINAL_NAME="${SAMPLE_NAME}_assembly_${assembler_name}_final.fasta"
            
            # Copy to eval output and genome assembly directory with clean name
            cp "$assembly" "${EVAL_OUTPUT_DIR}/${CLEAN_FINAL_NAME}"
            cp "$assembly" "${GENOME_ASSEMBLY_DIR}/${CLEAN_FINAL_NAME}"
            
            # Document assembly selection
            selection_message="$assembly was chosen as best assembly and renamed to \"${CLEAN_FINAL_NAME}\" for all further downstream analysis"
            echo -e "\n$selection_message" >> "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE"
            echo -e "\n$selection_message" >> "${EVAL_OUTPUT_DIR}/chosen_assembly.tsv"
            
            # Copy evaluation files
            cp "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE" "$GENOME_ASSEMBLY_DIR/"
            cp "${EVAL_OUTPUT_DIR}/chosen_assembly.tsv" "$GENOME_ASSEMBLY_DIR/"
            
            # Copy detailed QUAST output (use original file prefix for matching)
            prefix=${file%.*}
            quast_dir=$(find "$EVAL_OUTPUT_DIR" -type d -name "*${prefix}_quast_output*" -print -quit)
            if [[ -n "$quast_dir" ]]; then
                cp -r "$quast_dir" "$GENOME_ASSEMBLY_DIR/detailed_quast_output_of_chosen_assembly"
            fi
            
            log "$LOGS_DIR" "$LOG_NAME" "Reevaluation (III) complete. Best assembly saved as: ${CLEAN_FINAL_NAME}"
            break
        fi
    done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*${best_assembly}*")
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Assembly selection failed to return valid assembly name"
    exit 1
fi

# Clean up genome assembly directory: move intermediate assemblies to archive
# The final assembly (*_final.fasta) remains in the main directory
readonly INTERMEDIATE_ARCHIVE_DIR="$GENOME_ASSEMBLY_DIR/assembly_intermediates"
create_directory "$INTERMEDIATE_ARCHIVE_DIR"

while IFS= read -r -d '' file; do
    filename=$(basename "$file")
    
    # Keep files that end with _final.fasta
    if [[ "$filename" == *"_final.fasta" ]]; then
        continue
    fi
    
    # Move all other .fasta files to archive
    log "$LOGS_DIR" "$LOG_NAME" "Archiving intermediate: $filename"
    mv "$file" "$INTERMEDIATE_ARCHIVE_DIR/"
done < <(find "$GENOME_ASSEMBLY_DIR" -maxdepth 1 -type f -name "*.fasta" -print0)

# Log cleanup summary
final_count=$(find "$GENOME_ASSEMBLY_DIR" -maxdepth 1 -type f -name "*_final.fasta" | wc -l)
archived_count=$(find "$INTERMEDIATE_ARCHIVE_DIR" -type f -name "*.fasta" 2>/dev/null | wc -l)
log "$LOGS_DIR" "$LOG_NAME" "Final assemblies: $final_count, Archived intermediates: $archived_count"

# Process results with R
log "$LOGS_DIR" "$LOG_NAME" "Processing QUAST results with R"

# Process QUAST results
tsv_file=$(find "$GENOME_ASSEMBLY_DIR" -name 'chosen_assembly.tsv' -print -quit)
if [[ -n "$tsv_file" ]]; then
    apptainer exec \
        --bind "$R_SCRIPT_DIR:/mnt/r_script_dir" \
        --bind "$GENOME_ASSEMBLY_DIR:/mnt/input" \
        --bind "$QS_FILES_DIR:/mnt/output" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_quast.R" \
        --output_dir "/mnt/output" \
        --tsv "/mnt/input/$(basename "$tsv_file")" \
        --version "$VERSION"
fi

# Process selected assembly
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
if [[ -n "$analysis_assembly_file" ]]; then
    apptainer exec \
        --bind "$R_SCRIPT_DIR:/mnt/r_script_dir" \
        --bind "$(dirname "$analysis_assembly_file"):/mnt/input" \
        --bind "$QS_FILES_DIR:/mnt/output" \
        "$r_sif" \
        Rscript "/mnt/r_script_dir/R_process_selected_assembly.R" \
        --output_dir "/mnt/output" \
        --fasta "/mnt/input/$(basename "$analysis_assembly_file")" \
        --version "$VERSION"
fi