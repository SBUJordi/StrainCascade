#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_assembly_evaluation3.sh - Version 1.0.0
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

# Remove and recreate output directory
rm -rf "$EVAL_OUTPUT_DIR"
create_directory "$EVAL_OUTPUT_DIR"
create_directory "$QS_FILES_DIR"

# Process each assembly with QUAST
for assembly in "${assemblies[@]}"; do
    file=$(basename "$assembly")
    pattern_name=$(basename "$file" .fasta)
    quast_output_dir="${EVAL_OUTPUT_DIR}/${SAMPLE_NAME}_${pattern_name}_quast_output"
    
    create_directory "$quast_output_dir"
    log "$LOGS_DIR" "$LOG_NAME" "Running QUAST analysis for $file"
    
    apptainer exec \
        --bind "$GENOME_ASSEMBLY_DIR:/data" \
        --bind "$quast_output_dir:/output" \
        "$straincascade_assembly_qc_refinement_sif" \
        bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                 conda activate quast_env && \
                 quast.py /data/$file -o /output -t $THREADS"
    
    # Append QUAST report if available
    if [[ -f "$quast_output_dir/transposed_report.tsv" ]]; then
        cat "$quast_output_dir/transposed_report.tsv" >> "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE"
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
            prefix=${file%.*}
            
            # Copy and rename assembly files
            cp "$assembly" "${EVAL_OUTPUT_DIR}/${prefix}_best_ev3.fasta"
            cp "$assembly" "${GENOME_ASSEMBLY_DIR}/${prefix}_best_ev3.fasta"
            
            # Document assembly selection
            selection_message="$assembly was chosen as best assembly and renamed to \"${prefix}_best_ev3.fasta\" for all further downstream analysis"
            echo -e "\n$selection_message" >> "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE"
            echo -e "\n$selection_message" >> "${EVAL_OUTPUT_DIR}/chosen_assembly.tsv"
            
            # Copy evaluation files
            cp "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE" "$GENOME_ASSEMBLY_DIR/"
            cp "${EVAL_OUTPUT_DIR}/chosen_assembly.tsv" "$GENOME_ASSEMBLY_DIR/"
            
            # Copy detailed QUAST output
            if quast_dir=$(find "$EVAL_OUTPUT_DIR" -type d -name "*${prefix}_quast_output*" -print -quit); then
                cp -r "$quast_dir" "$GENOME_ASSEMBLY_DIR/detailed_quast_output_of_chosen_assembly"
            fi
            
            log "$LOGS_DIR" "$LOG_NAME" "Reevaluation (III) complete. Best assembly copied: ${prefix}_best_ev3.fasta"
            break
        fi
    done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*${best_assembly}*")
else
    log "$LOGS_DIR" "$LOG_NAME" "Error: Assembly selection failed to return valid assembly name"
    exit 1
fi

# Clean up genome assembly directory
declare -A file_categories
while IFS= read -r -d '' file; do
    if [[ $file =~ best_ev3 ]]; then
        file_categories["ev3"]+=" $file"
    elif [[ $file =~ best_ev2 ]]; then
        file_categories["ev2"]+=" $file"
    elif [[ $file =~ best_ev1 ]]; then
        file_categories["ev1"]+=" $file"
    else
        file_categories["none"]+=" $file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

# Determine files to keep
keep_files="${file_categories[none]}"
if [[ -n "${file_categories[ev3]:-}" ]]; then
    keep_files+=" ${file_categories[ev3]}"
elif [[ -n "${file_categories[ev2]:-}" ]]; then
    keep_files+=" ${file_categories[ev2]}"
elif [[ -n "${file_categories[ev1]:-}" ]]; then
    keep_files+=" ${file_categories[ev1]}"
fi

# Remove files not in keep list
while IFS= read -r -d '' file; do
    if [[ ! " $keep_files " =~ " $file " ]]; then
        rm "$file"
        log "$LOGS_DIR" "$LOG_NAME" "Removed: $file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

# Select final assembly
readonly PRIORITY_PATTERNS=('best_ev3' '_circularised\.fasta' 'best_ev2' '_mac2\.0_merged_assembly\.fasta' 'best_ev1')
selected_file=""

for pattern in "${PRIORITY_PATTERNS[@]}"; do
    if selected_file=$(find "$GENOME_ASSEMBLY_DIR" -type f -name "*${pattern}*" -print -quit); then
        break
    fi
done

# Fall back to first available file if no priority match
[[ -z "$selected_file" ]] && selected_file=$(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print -quit)

if [[ -n "$selected_file" ]]; then
    # Rename to final version
    new_file="${selected_file%.*}_final.fasta"
    mv "$selected_file" "$new_file"
    log "$LOGS_DIR" "$LOG_NAME" "Final assembly selected: $(basename "$new_file")"
    
    # Document final selection
    echo -e "\n$(basename "$selected_file") was chosen for all further downstream analysis and renamed to $(basename "$new_file")" >> \
        "${EVAL_OUTPUT_DIR}/$QUAST_REPORT_FILE" "${EVAL_OUTPUT_DIR}/chosen_assembly.tsv"
fi

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