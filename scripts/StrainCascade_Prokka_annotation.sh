#!/bin/bash

# StrainCascade_Prokka_annotation.sh - Version 1.3.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 11 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_assembly_main_abs> <genome_annotation_main_abs> <databases_dir> <locus_tag> <force_overwrite>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
threads=$6
genome_assembly_main_abs=$7
genome_annotation_main_abs=$8
databases_dir=$9
locus_tag=${10}
force_overwrite=${11}

# Load utils from the script directory
utils_file="${script_dir}/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

## Define paths and variables for this script ##
# List all matching .sif files and store them in an array
matching_files=($(ls "$apptainer_images_dir"/straincascade_genome_annotation*.sif 2> /dev/null))

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_genome_annotation=${matching_files[0]}

# Create output directory
prokka_output_dir="$output_dir/Prokka_annotation_results"
create_directory "$prokka_output_dir"  # Ensure the output directory is created

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (Prokka genome annotation) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for Prokka."
fi

# Run Prokka annotation
log "$logs_dir" "Prokka_annotation.log" "Running Prokka annotation for $analysis_assembly_file in $prokka_output_dir"

# Prepare Prokka command (--compliant mode sets minimal contig length to 200 bp which is NCBI-compliant: "Contigs should be >199nt, unless they are part of multi-component scaffolds in an AGP file"; https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/)
# INSDC Locus tag prefixes must contain between 3 and 12 alphanumeric uppercase characters and start with a letter.
# Set locus tag
if [ "$locus_tag" == "automatic" ]; then
    current_date=$(date +%y%m%d)
    random_digits=$(shuf -i 0-9 -n 3 | tr -d '\n')
    locus_tag="SC${random_digits}P${current_date}"
fi

prokka_cmd="prokka --outdir /mnt/output \
                   --prefix ${sample_name}_annotation_prokka \
                   --cpus $threads \
                   --compliant \
                   --addgenes \
                   --rfam \
                   --mincontiglen 200 \
                   --locustag $locus_tag"

# Add --force option if force_overwrite is yes
if [ "$force_overwrite" = "yes" ]; then
    prokka_cmd="${prokka_cmd} --force"
fi

# Add the input file path
prokka_cmd="$prokka_cmd /mnt/input/temp_input_assembly_Prokka.fasta"

# Run Prokka using Apptainer
apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$prokka_output_dir":/mnt/output \
    "$straincascade_genome_annotation" \
    /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate prokka_env && \
                  
                  # Copy and rename the input file (too long names may cause problems; this way the name of the input files is standardised)
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_Prokka.fasta && \
                  
                  # Set TMPDIR environment variable
                  export TMPDIR=/tmp && \

                  # Run Prokka command
                  $prokka_cmd && \
                  
                  # Remove the temporary file
                  rm /mnt/input/temp_input_assembly_Prokka.fasta" 2>&1

# Check exit status of Prokka
if [ $? -ne 0 ]; then
    log "$logs_dir" "Prokka_annotation.log" "Error: Prokka annotation failed for $analysis_assembly_file"
    exit 1
fi

# Copy the output files to genome_annotation_main_abs
for ext in _annotation_prokka.gff _annotation_prokka.gbk _annotation_prokka.fna _annotation_prokka.faa _annotation_prokka.ffn _annotation_prokka.sqn _annotation_prokka.fsa _annotation_prokka.tbl _annotation_prokka.err _annotation_prokka.log _annotation_prokka.txt _annotation_prokka.tsv; do
    files=$(find "$prokka_output_dir" -type f -name "*.$ext")
    if [ -n "$files" ]; then
        cp $files "$genome_annotation_main_abs"
    else
        echo "Warning: No .$ext files found in $prokka_output_dir"
    fi
done