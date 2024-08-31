#!/bin/bash

# StrainCascade_MAC2_assembly_merging.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <output_directory> <sample_name> <genome_assembly_main_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
output_dir=$5
sample_name=$6
genome_assembly_main_abs=$7

source "$utils_file"

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'straincascade_assembly_qc_refinement*.sif' -print)

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_assembly_qc_refinement=${matching_files[0]}

# Define paths and variables
mac2_output_dir="${output_dir}/MAC2_optimization_results"
create_directory "$mac2_output_dir"
create_directory "${mac2_output_dir}/input"
create_directory "${mac2_output_dir}/output"

# Function to count contigs in a FASTA file
count_contigs() {
    grep -c "^>" "$1"
}

# Find the best original assembly
reference_assembly=$(find "$genome_assembly_main_abs" -type f -name "*_best_ev1.fasta")

# If no best original assembly is found, choose the one with the lowest contig count
if [ ! -f "$reference_assembly" ]; then
    log "$logs_dir" "MAC2_optimization.log" "Best original assembly ("*_best_ev1.fasta") not found. Choosing assembly with lowest contig count as reference."
    reference_assembly=$(find "$genome_assembly_main_abs" -type f -name "*.fasta" | 
                             sort -n -k2 -t: <(while read -r file; do 
                                 echo "$file:$(count_contigs "$file")"; 
                             done < <(find "$genome_assembly_main_abs" -type f -name "*.fasta")) | 
                             head -n1 | cut -d: -f1)
fi

if [ ! -f "$reference_assembly" ]; then
    log "$logs_dir" "MAC2_optimization.log" "No suitable reference assembly found. Skipping MAC2.0 assembly merging."
    exit 0
fi

# Count contigs in the reference assembly
contig_count=$(count_contigs "$reference_assembly")

if [ "$contig_count" -gt 1 ]; then
    log "$logs_dir" "MAC2_optimization.log" "Reference assembly has multiple contigs. Proceeding with MAC2.0 iterative assembly merging."
    
    # Find all assembly files except the reference assembly
    assembly_files=($(find "$genome_assembly_main_abs" -type f -name "*.fasta" ! -name "$(basename "$reference_assembly")"))
    
    # Sort assemblies by contig count
    sorted_assemblies=($(for file in "${assembly_files[@]}"; do
        echo "$file:$(count_contigs "$file")"
    done | sort -n -k2 -t: | cut -d: -f1))
    
    # Copy reference assembly to input directory
    cp "$reference_assembly" "${mac2_output_dir}/input/reference.fa"
    
    # Iterative merging
    merged_assembly="${mac2_output_dir}/input/reference.fa"
    for assembly in "${sorted_assemblies[@]}"; do
        assembly_name=$(basename "$assembly" .fasta)
        
        # Copy query assembly to input directory
        cp "$assembly" "${mac2_output_dir}/input/query.fa"
        
        log "$logs_dir" "MAC2_optimization.log" "Merging $assembly_name with reference"
        
        apptainer exec \
            --bind "${mac2_output_dir}":/mnt/mac2 \
            "$straincascade_assembly_qc_refinement" \
            /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                          conda activate mac2_env && \
                          cd /mnt/mac2 && \
                          MAC2.0 input/query.fa input/$(basename "$merged_assembly")"
        
        # Update merged assembly for next iteration
        merged_assembly="${mac2_output_dir}/output/merged_contigs.fa"
        mv "$merged_assembly" "${mac2_output_dir}/input/reference.fa"
    done
    
    # Rename the final merged assembly
    mv "${mac2_output_dir}/input/reference.fa" "${mac2_output_dir}/${sample_name}_mac2.0_merged_assembly.fasta"
    cp "${mac2_output_dir}/${sample_name}_mac2.0_merged_assembly.fasta" "${genome_assembly_main_abs}/${sample_name}_mac2.0_merged_assembly.fasta"
else
    log "$logs_dir" "MAC2_optimization.log" "Reference assembly has only one contig (= already good enough). Skipping MAC2.0 assembly merging."
fi