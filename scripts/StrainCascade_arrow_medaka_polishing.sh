#!/bin/bash

# StrainCascade_arrow_medaka_polishing.sh - Version 1.0.4
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 11 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <input_file> <bam_file> <output_dir> <sample_name> <sequencing_type> <threads> <genome_assembly_main_abs>"
    exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
input_file=$5
bam_file=$6
output_dir=$7
sample_name=$8
sequencing_type=$9
threads=${10}
genome_assembly_main_abs=${11}

source "$utils_file"

# Determine polishing tool based on sequencing type
if [[ "$sequencing_type" == *"pacbio"* && -n "$bam_file" && "$bam_file" != "not_available" ]]; then
    polishing_tool="arrow"
elif [[ "$sequencing_type" == *"nano"* ]]; then
    polishing_tool="medaka"
else
    echo "No valid sequencing type or necessary BAM file identified. Skipping this module and continuing with the next script in the pipeline."
    exit 0
fi

# Find matching .sif files
matching_files=($(ls "$apptainer_images_dir"/straincascade_assembly_qc_refinement*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Skipping this module and continuing with the next script in the pipeline."
    exit 0
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_assembly_qc_refinement=${matching_files[0]}

# Create output directory
polishing_output_dir="$output_dir/${polishing_tool}_polishing_results"
if ! create_directory "$polishing_output_dir"; then
    echo "Error: Failed to create output directory. Skipping this module and continuing with the next script in the pipeline."
    exit 0
fi

# Retrieve the analysis assembly file from genome_assembly_main_abs using the new function
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file assembly for analysis."
fi

# Check if input files exist
if [ ! -f "$input_file" ]; then
    echo "Error: Input file $input_file does not exist. Skipping this module and continuing with the next script in the pipeline."
    exit 0
fi
if [ "$polishing_tool" = "arrow" ] && [ ! -f "$bam_file" ]; then
    echo "Error: BAM file $bam_file does not exist. Skipping this module and continuing with the next script in the pipeline."
    exit 0
fi

# Standardize RG IDs and add necessary tags in the BAM file
echo "Standardizing RG IDs and adding necessary tags in the BAM file..."
log "$logs_dir" "${polishing_tool}_polishing.log" "Standardizing RG IDs and adding necessary tags in the BAM file..."

fixed_bam_file="${polishing_output_dir}/fixed_$(basename "$bam_file")"

# Run the standardization process inside the Apptainer
if ! apptainer exec \
    --bind "$(dirname "$bam_file")":/mnt/input \
    --bind "$polishing_output_dir":/mnt/output \
    "$straincascade_assembly_qc_refinement" \
    /bin/bash -c "
    set -e && \
    new_rg_id=$(printf "%08x" $RANDOM) && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate tools_env && \
    samtools addreplacerg \
    -r \"ID:\${new_rg_id}\" \
    -r \"PL:PACBIO\" \
    -r \"PU:unknown\" \
    -r \"LB:unknown\" \
    -r \"SM:sample\" \
    -o /mnt/output/fixed_\$(basename \"$bam_file\") \
    /mnt/input/\$(basename \"$bam_file\") && \
    samtools index /mnt/output/fixed_\$(basename \"$bam_file\") && \
    samtools view -H /mnt/output/fixed_\$(basename \"$bam_file\")
    "
then
    echo "Error: Failed to standardize RG IDs and add necessary tags in the BAM file. Skipping this module and continuing with the next script in the pipeline."
    exit 0
fi

# Update the bam_file variable to use the fixed version
bam_file="$fixed_bam_file"
echo "BAM file with standardized RG IDs and necessary tags: $bam_file"
log "$logs_dir" "${polishing_tool}_polishing.log" "BAM file with standardized RG IDs and necessary tags: $bam_file"

# Create PacBio index file (.pbi)
echo "Creating PacBio index file..."
log "$logs_dir" "${polishing_tool}_polishing.log" "Creating PacBio index file..."

if ! apptainer exec \
    --bind "$(dirname "$bam_file")":/mnt/bam_file_dir \
    "$straincascade_assembly_qc_refinement" \
    /bin/bash -c "
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate tools_env && \
    pbindex /mnt/bam_file_dir/$(basename "$bam_file")
    "
then
    echo "Error: Failed to create PacBio index file. Skipping this module and continuing with the next script in the pipeline."
    exit 0
fi

echo "PacBio index file created successfully."
log "$logs_dir" "${polishing_tool}_polishing.log" "PacBio index file created successfully."

# Polishing rounds
polishing_rounds=2
current_assembly="$analysis_assembly_file"

for ((round=1; round<=polishing_rounds; round++)); do
    echo "Starting polishing round $round"
    log "$logs_dir" "${polishing_tool}_polishing.log" "Starting polishing round $round"

    round_output_dir="${polishing_output_dir}/round_${round}"
    if ! create_directory "$round_output_dir"; then
        echo "Error: Failed to create round output directory. Skipping this module and continuing with the next script in the pipeline."
        exit 0
    fi

    if [ "$polishing_tool" = "medaka" ]; then
        polishing_cmd="conda activate medaka_env && \
        medaka_consensus \
                -i /mnt/input_file/$(basename "$input_file") \
                -d /mnt/input_assembly/$(basename "$current_assembly") \
                -o /mnt/output -t $threads"
    elif [ "$polishing_tool" = "arrow" ]; then
        polishing_cmd="
        conda activate genomicconsensus_env && \
        variantCaller --algorithm=arrow \
        -j $threads \
        -r /mnt/input_assembly/$(basename "$current_assembly") \
        -o /mnt/output/consensus.fasta \
        /mnt/bam_file_dir/$(basename "$bam_file")"
    fi

    # Run the polishing command
    if ! apptainer exec \
        --bind "$(dirname "$input_file")":/mnt/input_file \
        --bind "$(dirname "$current_assembly")":/mnt/input_assembly \
        --bind "$round_output_dir":/mnt/output \
        --bind "$(dirname "$bam_file")":/mnt/bam_file_dir \
        "$straincascade_assembly_qc_refinement" \
        /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
        $polishing_cmd" 2>&1
    then
        echo "Error: Round $round of ${polishing_tool} polishing failed. Skipping this module and continuing with the next script in the pipeline."
        exit 0
    fi

    echo "Round $round of ${polishing_tool} polishing completed successfully."
    log "$logs_dir" "${polishing_tool}_polishing.log" "Round $round of ${polishing_tool} polishing completed successfully."
    
    current_assembly="$round_output_dir/consensus.fasta"
done

echo "Final polished assembly: $current_assembly"
log "$logs_dir" "${polishing_tool}_polishing.log" "Final polished assembly: $current_assembly"

# Overwrite the initial assembly file with the polished assembly
if [ -f "$current_assembly" ]; then
    cp "$current_assembly" "$analysis_assembly_file"
    echo "Initial assembly file overwritten with polished assembly."
    log "$logs_dir" "${polishing_tool}_polishing.log" "Initial assembly file overwritten with polished assembly."
else
    echo "Error: Final polished assembly not found. Initial assembly file not overwritten."
    log "$logs_dir" "${polishing_tool}_polishing.log" "Error: Final polished assembly not found. Initial assembly file not overwritten."
fi

# Cleanup
if [ "$polishing_tool" = "arrow" ]; then
    echo "Cleaning up temporary BAM file..."
    log "$logs_dir" "${polishing_tool}_polishing.log" "Cleaning up temporary BAM file..."
    rm -f "$fixed_bam_file" "${fixed_bam_file}.bai"
fi

echo "Polishing module completed successfully."
log "$logs_dir" "${polishing_tool}_polishing.log" "Polishing module completed successfully."

exit 0