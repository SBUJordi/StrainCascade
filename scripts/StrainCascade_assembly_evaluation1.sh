#!/bin/bash

# StrainCascade_assembly_evaluation1.sh - Version 1.1.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <script_dir> <logs_dir> <apptainer_images_dir> <output_directory> <sample_name> <threads> <genome_assembly_main_abs>"
    exit 1
fi
script_dir=$1
logs_dir=$2
apptainer_images_dir=$3
output_dir=$4
sample_name=$5
threads=$6
genome_assembly_main_abs=$7

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

# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'python_3.12.4*.sif' -print)

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_python=${matching_files[0]}

# Create output directory and copy assembly files
assembly_evaluation_output_dir="$output_dir/StrainCascade_assembly_evaluation"
create_directory "$assembly_evaluation_output_dir"

# Find only the original assembly generated directly by assemblers in genome_assembly_main_abs
assemblies=($(find "$genome_assembly_main_abs" -type f -name "*.fasta" ! -name "*_best_ev*" ! -name "*_circularised.fasta" ! -name "*_mac2.0_merged_assembly.fasta"))

# Check if any matching files were found
if [ ${#assemblies[@]} -eq 0 ]; then
    error_message="No matching files (=assemblies created by StrainCascade) found in '$genome_assembly_main_abs'. Skipping this evaluation step. Possible issue: No assembler module was run or any suitable assembler module should be run first, respectively."
    echo "$error_message"
    log "$logs_dir" "StrainCascade_assembly_evaluation.log" "$error_message"
    exit 0  # Exit gracefully, allowing the pipeline to continue
fi

# Process each assembly file
for assembly in "${assemblies[@]}"
do
  # Determine the file name and pattern name
  file=$(basename "$assembly")
  pattern_name=$(basename "$file" .fasta)
  quast_output_dir="${assembly_evaluation_output_dir}/${sample_name}_${pattern_name}_quast_output"
  create_directory "$quast_output_dir"

  log "$logs_dir" "StrainCascade_assembly_evaluation.log" "Running QUAST in Apptainer for $file"

  # Run QUAST using Apptainer, binding both the original assembly directory and the output directory
  apptainer exec --bind "$genome_assembly_main_abs:/data" \
                --bind "$quast_output_dir:/output" \
                "$straincascade_assembly_qc_refinement" \
                /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                              conda activate quast_env && \
                              quast.py /data/$file -o /output -t \"$threads\""

  # Append QUAST report to the assembly evaluation file
  if [ -f "$quast_output_dir/transposed_report.tsv" ]; then
    cat "$quast_output_dir/transposed_report.tsv" >> "${assembly_evaluation_output_dir}/${sample_name}_quast_assembly_evaluation.tsv"
  else
    log "$logs_dir" "StrainCascade_assembly_evaluation.log" "Warning: QUAST did not produce output for $file"
  fi
done

# Run the assembly selection script using Apptainer
log "$logs_dir" "StrainCascade_assembly_evaluation.log" "Running assembly selection script in Apptainer"
best_assembly=$(apptainer exec --bind "$assembly_evaluation_output_dir:/data" \
                                 --bind "$script_dir:/scripts" \
                                 "$straincascade_python" \
                                 /bin/bash -c "python /scripts/StrainCascade_assembly_selection.py /data/${sample_name}_quast_assembly_evaluation.tsv /data")

if [ -n "$best_assembly" ]; then
  echo "Best Assembly (ev1): $best_assembly"
  log "$logs_dir" "StrainCascade_assembly_evaluation.log" "Best Assembly (ev1): $best_assembly"
else
  echo "Error: Assembly selection script did not return a valid assembly name"
  log "$logs_dir" "StrainCascade_assembly_evaluation.log" "Error: Assembly selection script did not return a valid assembly name"
  exit 1
fi

# Find the original assembly file that matches the pattern
for assembly in $(find "$genome_assembly_main_abs" -type f -name "*${best_assembly}*")
do
  if [ -f "$assembly" ]; then
    # Copy the assembly file to the output directory with a new name
    file=$(basename "$assembly")
    prefix=${file%.*}
    cp "$assembly" "${assembly_evaluation_output_dir}/${prefix}_best_ev1.fasta"
    
    # Copy the assembly file to the genome_assembly_main_abs with a new name
    cp "$assembly" "${genome_assembly_main_abs}/${prefix}_best_ev1.fasta"

    break  # Break the loop as soon as the file is found and copied
  fi
done

echo "Evaluation (I) complete. Best evaluated assembly: ${prefix}_best_ev1.fasta"