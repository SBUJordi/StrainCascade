#!/bin/bash

# StrainCascade_Bakta_annotation.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 12 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_dir> <sample_name> <threads> <genome_assembly_main_abs> <taxonomic_classification_main_abs> <genome_annotation_main_abs> <databases_dir> <locus_tag> <force_overwrite>"
    exit 1
fi

script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
threads=$6
genome_assembly_main_abs=$7
taxonomic_classification_main_abs=$8
genome_annotation_main_abs=$9
databases_dir=${10}
locus_tag=${11}
force_overwrite=${12}

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
bakta_output_dir="$output_dir/Bakta_annotation_results"
create_directory "$bakta_output_dir"  

# Retrieve the analysis assembly file from genome_assembly_main_abs
analysis_assembly_file=$(find_analysis_assembly_file "$genome_assembly_main_abs")
if [ -z "$analysis_assembly_file" ]; then
    echo "Error: No assembly files found. Skipping this module (Bakta genome annotation) and continuing with the next script in the pipeline."
    exit 0
else
    echo "Using $analysis_assembly_file for Bakta."
fi

# Run Bakta annotation
log "$logs_dir" "Bakta_annotation.log" "Running Bakta annotation for $analysis_assembly_file in $bakta_output_dir"

# Check the taxonomic level
taxonomic_level_file=$(ls $taxonomic_classification_main_abs/*level_organism_name.txt)
taxonomic_level=$(cat "$taxonomic_level_file")

# Check the organism name
organism_name_file=$(ls $taxonomic_classification_main_abs/*gtdbtk_organism_name.txt)
organism_name=$(cat "$organism_name_file")

# Prepare Bakta command based on taxonomic level (--compliant mode sets minimal contig length to 200 bp which is NCBI-compliant: "Contigs should be >199nt, unless they are part of multi-component scaffolds in an AGP file"; https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/)
# INSDC Locus tag prefixes must contain between 3 and 12 alphanumeric uppercase characters and start with a letter.
# Set locus tag
if [ "$locus_tag" == "automatic" ]; then
    current_date=$(date +%y%m%d)
    random_digits=$(shuf -i 0-9 -n 3 | tr -d '\n')
    locus_tag="SC${random_digits}B${current_date}"
fi

# Prepare Bakta command
bakta_cmd="bakta --db /mnt/bakta_db/db \
                 --output /mnt/output \
                 --prefix ${sample_name}_annotation_bakta \
                 --compliant \
                 --threads $threads \
                 --locus-tag $locus_tag \
                 --tmp-dir /tmp" 

# Add appropriate taxonomic info option to Bakta command if available
if [ "$taxonomic_level" == "genus" ]; then
    bakta_cmd="$bakta_cmd --genus \"$organism_name\""
elif [ "$taxonomic_level" == "species" ]; then
    bakta_cmd="$bakta_cmd --species \"$organism_name\""
fi

# Add --force option if force_overwrite is yes
if [ "$force_overwrite" = "yes" ]; then
    bakta_cmd="$bakta_cmd --force"
fi

# Add the input file path
bakta_cmd="$bakta_cmd /mnt/input/temp_input_assembly_Bakta.fasta"

apptainer exec \
    --bind "$(dirname "$analysis_assembly_file")":/mnt/input \
    --bind "$bakta_output_dir":/mnt/output \
    --bind "$databases_dir/bakta_db/db":/mnt/bakta_db/db \
    "$straincascade_genome_annotation" \
    /bin/bash -c "
                  source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate bakta_env && \
                  
                  # Copy and rename the input file (too long names may cause problems; this way the name of the input files is standardised)
                  cp /mnt/input/$(basename "$analysis_assembly_file") /mnt/input/temp_input_assembly_Bakta.fasta && \
                  
                  # Set TMPDIR environment variable
                  export TMPDIR=/tmp && \
                  
                  # Update AMRFinderPlus's internal database (otherwise error)
                  amrfinder_update --force_update --database /mnt/bakta_db/db/amrfinderplus-db && \
                  
                  # Run Bakta command
                  $bakta_cmd && \
                  
                  # Remove the temporary file
                  rm /mnt/input/temp_input_assembly_Bakta.fasta" 2>&1

# Check exit status of Bakta
if [ $? -ne 0 ]; then
    log "$logs_dir" "Bakta_annotation.log" "Error: Bakta annotation failed for $analysis_assembly_file"
    exit 1
fi

# Copy the output files to genome_annotation_main_abs
for ext in annotation_bakta.tsv annotation_bakta.ffn annotation_bakta.faa annotation_bakta.gbff annotation_bakta.png annotation_bakta.svg; do
    files=$(find "$bakta_output_dir" -type f -name "*.$ext")
    if [ -n "$files" ]; then
        cp $files "$genome_annotation_main_abs"
    else
        echo "Warning: No .$ext files found in $bakta_output_dir"
    fi
done