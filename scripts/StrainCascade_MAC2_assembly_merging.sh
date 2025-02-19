#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_MAC2_assembly_merging.sh
# Description: Performs MAC2 optimization and merging of genome assemblies as part of StrainCascade

set -euo pipefail

# Function to display usage information
show_usage() {
    cat << EOF
Usage: $0 <script_dir> <logs_dir> <log_name> <utils_file> <apptainer_images_dir> 
          <output_directory> <sample_name> <genome_assembly_main_abs>
EOF
    exit 1
}

# Validate input parameters
[[ $# -eq 8 ]] || show_usage

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_DIR="$5"
readonly OUTPUT_DIR="$6"
readonly SAMPLE_NAME="$7"
readonly GENOME_ASSEMBLY_DIR="$8"

# Derived constants
readonly MAC2_OUTPUT_DIR="${OUTPUT_DIR}/MAC2_optimization_results"
readonly MAC2_WORKDIR="${MAC2_OUTPUT_DIR}/workdir"
readonly MAX_RUNTIME_SECONDS=720000  # 200 hours -> adjust if you want real time-constraints
readonly FILTER_THRESHOLD_SECONDS=86400  # 24 hours
readonly MIN_ASSEMBLIES=5

# Source utility functions
source "$UTILS_FILE"

# Define additional utility functions
count_contigs() {
    local file="$1"
    
    if [[ -z "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: No file for contig counting provided"
        echo "0"
        return 0
    fi
    if [[ ! -f "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: File for contig counting not found: $file"
        echo "0"
        return 0
    fi
    
    # Count headers based on file type
    if [[ "$file" =~ \.gz$ ]]; then
        zcat "$file" | awk '/^>/ { count++ } END { print count+0 }'
    else
        awk '/^>/ { count++ } END { print count+0 }'  "$file"
    fi
}

# Create necessary directory structure
create_directory "$MAC2_OUTPUT_DIR"
create_directory "${MAC2_WORKDIR}"
create_directory "${MAC2_WORKDIR}/input"
create_directory "${MAC2_WORKDIR}/output"
create_directory "${MAC2_WORKDIR}/temp"

# Clean up empty FASTA files
log "$LOGS_DIR" "$LOG_NAME" "Checking for empty FASTA files in $GENOME_ASSEMBLY_DIR"
while IFS= read -r -d '' file; do
    if [[ -f "$file" && ! -s "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Removing empty file: $(basename "$file")"
        rm "$file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

# Find MAC2 container
readonly straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_DIR" 'straincascade_assembly_qc_refinement*.sif')

# Try to find best_ev1 assembly first
initial_best_assembly=$(find "$GENOME_ASSEMBLY_DIR" -type f -name "*best_ev1*.fasta" | head -n 1)

if [[ -z "$initial_best_assembly" ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No best_ev1 assembly found, searching for assembly with least contigs..."
    
    readarray -t fasta_files < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta")
    
    if [[ ${#fasta_files[@]} -eq 0 ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "No assembly found in $GENOME_ASSEMBLY_DIR that could serve as query. Skipping MAC2.0 assembly merging."
        exit 0
    fi
    
    lowest_count=999999
    lowest_file=""
    
    for file in "${fasta_files[@]}"; do
        if ! count=$(count_contigs "$file"); then
            log "$LOGS_DIR" "$LOG_NAME" "Warning: Failed to count contigs in $file, skipping"
            continue
        fi
        
        if [[ $count -lt $lowest_count ]]; then
            lowest_count=$count
            lowest_file=$file
        fi
    done
    
    if [[ -z "$lowest_file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Could not find valid assembly with countable contigs"
        exit 0
    fi
    
    initial_best_assembly=$lowest_file
    log "$LOGS_DIR" "$LOG_NAME" "Selected available assembly with least contigs ($lowest_count): $(basename "$lowest_file")"
fi

# Check if initial best assembly has only one contig
if ! initial_contig_count=$(count_contigs "$initial_best_assembly"); then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Could not count contigs in best assembly"
    exit 0
fi

if [[ $initial_contig_count -eq 1 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Initial best assembly already has only one contig. No merging needed."
    exit 0
fi

# Copy initial best assembly as query
cp "$initial_best_assembly" "${MAC2_WORKDIR}/input/query.fa"

# Find remaining assemblies excluding best assembly
readarray -t assembly_files < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" ! -name "$(basename "$initial_best_assembly")")

if [[ ${#assembly_files[@]} -eq 0 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No additional assemblies found for merging. Skipping MAC2.0 assembly merging."
    exit 0
fi

# Create array of tuples (file,contig_count) for sorting
declare -a assembly_tuples
for file in "${assembly_files[@]}"; do
    if ! count=$(count_contigs "$file"); then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Failed to count contigs in $file, skipping"
        continue
    fi
    assembly_tuples+=("$count:$file")
done

if [[ ${#assembly_tuples[@]} -eq 0 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "No valid additional assemblies found after contig counting. Skipping MAC2.0 assembly merging."
    exit 0
fi

# Sort by contig count
readarray -t sorted_assemblies < <(printf '%s\n' "${assembly_tuples[@]}" | sort -t':' -n | cut -d':' -f2-)

log "$LOGS_DIR" "$LOG_NAME" "Found ${#sorted_assemblies[@]} valid additional assemblies, ordered by contig count"

# Initialize array for ordered basenames
declare -a ordered_reference_assemblies

# Copy sorted assemblies to work directory with sequential naming as references
for index in "${!sorted_assemblies[@]}"; do
    src_file="${sorted_assemblies[$index]}"
    basename_file=$(basename "$src_file" .fasta).fa
    ordered_reference_assemblies+=("$basename_file")
    
    # Copy file as reference instead of query
    if ! cp "$src_file" "${MAC2_WORKDIR}/input/reference_${index}.fa"; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to copy $basename_file to work directory. Trying next assembly..."
        continue
    fi
done

log "$LOGS_DIR" "$LOG_NAME" "Copied ${#ordered_reference_assemblies[@]} reference assemblies to work directory"
log "$LOGS_DIR" "$LOG_NAME" "List of copied files: $(IFS=,; echo "${ordered_reference_assemblies[*]}")"

# Initialize filtering variables
filtering_enabled=false
max_contig_count=0

# Start time tracking
readonly SCRIPT_START_TIME=$(date +%s)

# Iterate through numbered reference assemblies
for index in "${!ordered_reference_assemblies[@]}"; do
    current_time=$(date +%s)
    total_runtime=$((current_time - SCRIPT_START_TIME))
    
    # Check if we should switch to filtered mode
    if [[ $total_runtime -gt $FILTER_THRESHOLD_SECONDS && $filtering_enabled == false ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Switching to filtered mode after ${FILTER_THRESHOLD_SECONDS} seconds"
        filtering_enabled=true
        
        # Calculate filtering threshold based on current contigs
        declare -a current_counts
        for tuple in "${assembly_tuples[@]}"; do
            count="${tuple%%:*}"
            current_counts+=("$count")
        done
        
        if [[ ${#current_counts[@]} -ge $MIN_ASSEMBLIES ]]; then
            sorted_counts=($(printf '%s\n' "${current_counts[@]}" | sort -n))
            mid=$(( ${#sorted_counts[@]} / 2 ))
            median=$(( (${sorted_counts[mid]} + ${sorted_counts[mid - 1]}) / 2 ))
            max_contig_count=$((median * 3))
        else
            sorted_counts=($(printf '%s\n' "${current_counts[@]}" | sort -n))
            mid=$(( ${#sorted_counts[@]} / 2 ))
            max_contig_count=$((${sorted_counts[mid]} * 3))
        fi
        
        log "$LOGS_DIR" "$LOG_NAME" "Set maximum contig count threshold to $max_contig_count"
    fi
    
    # Exit if we've exceeded maximum runtime
    if [[ $total_runtime -gt $MAX_RUNTIME_SECONDS ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Maximum runtime exceeded. Finalizing MAC2.0 assembly merging."
        break
    fi
    
    # Use current reference file and assembly name
    reference_file="reference_${index}.fa"
    assembly_name="${ordered_reference_assemblies[$index]}"
    
    # Skip if filtering is enabled and count is too high
    if [[ $filtering_enabled == true ]]; then
        count="${assembly_tuples[$index]%%:*}"
        if [[ $count -gt $max_contig_count ]]; then
            log "$LOGS_DIR" "$LOG_NAME" "Skipping ${assembly_name} - too many contigs ($count > $max_contig_count)"
            continue
        fi
    fi
    
    log "$LOGS_DIR" "$LOG_NAME" "Starting MAC2.0 with reference: ${assembly_name}"
    
    # Run MAC2.0 with proper paths
    if ! apptainer exec \
        --bind "${MAC2_WORKDIR}":/mnt/mac2 \
        "$straincascade_assembly_qc_refinement_sif" \
        /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                      conda activate MAC2.0_env && \
                      cd /mnt/mac2 && \
                      MAC2.0 query.fa ${reference_file}"; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: MAC2.0 process failed for $assembly_name. Trying next assembly..."
        continue
    fi
    
    # Check if merge produced output
    if [[ ! -f "${MAC2_WORKDIR}/output/scaffold.fasta" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Warning: No merged output produced for $assembly_name. Trying next assembly..."
        continue
    fi
    
    # Update query for next iteration
    mv "${MAC2_WORKDIR}/output/scaffold.fasta" "${MAC2_WORKDIR}/input/query.fa"
done

# Save final merged assembly from query file
if [[ -f "${MAC2_WORKDIR}/input/query.fa" ]]; then
    mv "${MAC2_WORKDIR}/input/query.fa" "${MAC2_OUTPUT_DIR}/${SAMPLE_NAME}_mac2.0_merged_assembly.fasta"
    cp "${MAC2_OUTPUT_DIR}/${SAMPLE_NAME}_mac2.0_merged_assembly.fasta" "${GENOME_ASSEMBLY_DIR}/${SAMPLE_NAME}_mac2.0_merged_assembly.fasta"
    log "$LOGS_DIR" "$LOG_NAME" "MAC2.0 assembly merging completed successfully"
else
    log "$LOGS_DIR" "$LOG_NAME" "MAC2.0 assembly merging produced no final merged assembly"
fi

# Clean up temporary files
rm -rf "${MAC2_WORKDIR}"

# Clean up empty FASTA files
log "$LOGS_DIR" "$LOG_NAME" "Checking for empty FASTA files in $GENOME_ASSEMBLY_DIR"
while IFS= read -r -d '' file; do
    if [[ -f "$file" && ! -s "$file" ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Removing empty file: $(basename "$file")"
        rm "$file"
    fi
done < <(find "$GENOME_ASSEMBLY_DIR" -type f -name "*.fasta" -print0)

exit 0