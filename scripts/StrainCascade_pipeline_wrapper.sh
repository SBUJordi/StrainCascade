#!/bin/bash

# Wrapper script:
# StrainCascade_WG_seq_pipeline.sh - Version 1.1.0
# Author: Sebastian Bruno Ulrich Jordi

input_file=$1
output_directory=$2
sequencing_type=$3
threads=$4
selected_modules=$5
result_type=$6
force_overwrite=$7
databases_dir=$8
locus_tag=$9

## Define variables and create directories ##
# Convert the selected_modules string to an array
IFS=' ' read -r -a selected_modules <<< "$selected_modules"

# Get the directory of this script (don't move the pipeline scripts from the original location)
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Load utils from the script directory
utils_file="$script_dir/utils.sh"
if [ -f "$utils_file" ]; then
  source "$utils_file"
else
  echo "Error: utils.sh not found in $script_dir"
  exit 1
fi

# Define the pipeline directory (= StrainCascade installation directory)
pipeline_dir="$(dirname "$script_dir")"

# Define the apptainer images directory
apptainer_images_dir="$pipeline_dir/apptainer_images"

# Define the sample name variable
sample_name="$input_file"
for suffix in .fasta .fa .fastq .fastq.gz .fna .bam; do
  sample_name=$(basename "$sample_name" "$suffix")
done

# Define base_dir=$output_directory (to be used by the pipeline scripts)
base_dir="$output_directory"

# Create the output_dir (the output_dir is derived from the argument output_directory and will be the actual folder to contain results)
output_dir="${base_dir}/output_${sample_name}"
create_directory "$output_dir"

# Define the main results directory
main_results_dir_abs="$(get_absolute_path "${output_dir}/${sample_name}_StrainCascade_main_results")"
create_directory "$main_results_dir_abs"
# Create subdirectories in the main results directory
subdirs=("00_sequencing_reads" "01_genome_assembly" "02_taxonomic_classification" "03_genome_annotation" "04_functional_analysis" "05_results_integration")
for subdir in "${subdirs[@]}"; do
  create_directory "${main_results_dir_abs}/${subdir}"
done

# Define absolute path variables for the subdirectories
sequencing_reads_main_abs="$(get_absolute_path "${main_results_dir_abs}/00_sequencing_reads")"
genome_assembly_main_abs="$(get_absolute_path "${main_results_dir_abs}/01_genome_assembly")"
taxonomic_classification_main_abs="$(get_absolute_path "${main_results_dir_abs}/02_taxonomic_classification")"
genome_annotation_main_abs="$(get_absolute_path "${main_results_dir_abs}/03_genome_annotation")"
functional_analysis_main_abs="$(get_absolute_path "${main_results_dir_abs}/04_functional_analysis")"
results_integration_abs="$(get_absolute_path "${main_results_dir_abs}/05_results_integration")"

# Define the logs directory
logs_dir="$output_dir/pipeline_logs"
create_directory "$logs_dir"

## Handle input file types
# List all matching .sif files and store them in an array
readarray -t matching_files < <(find "$apptainer_images_dir" -name 'straincascade_genome_assembly*.sif' -print)

# Check the number of matching files
if [ ${#matching_files[@]} -eq 0 ]; then
    echo "No matching .sif files found in $apptainer_images_dir. Continuing with the next script in the pipeline."
    exit 0  # Exit gracefully, allowing the pipeline to continue
elif [ ${#matching_files[@]} -gt 1 ]; then
    echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
fi

# Proceed with the first match
straincascade_genome_assembly=${matching_files[0]}

if [ -n "$input_file" ]; then
  # Extract the file name, and extension
  dir=$(dirname "$input_file")
  base=$(basename "$input_file")
  extension="${base##*.}"
  filename="${base%.*}"
  
  # Check if the extension is one of the accepted formats
  if [[ "$extension" =~ ^(fasta|fa|fastq|fastq.gz|fna|bam)$ ]]; then
    if [ "$extension" = "bam" ]; then
      
      bam_file="$input_file"
      
      apptainer exec \
        --bind "$dir":/mnt/input \
        "$straincascade_genome_assembly" \
        /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                  conda activate tools_env && \
                  samtools fasta /mnt/input/$base > /mnt/input/${filename}.fasta" 2>&1
      
      input_file="${dir}/${filename}.fasta"

    else 
      bam_file="not_available"
    fi
  else
    echo "Error: Unsupported file extension"
  fi
else
  echo "Error: No input file identified"
fi


## Start processing ##
# Log the attempt to change directory
log "$logs_dir" "pipeline.log" "Changing directory to: $base_dir"
# Change directory to base_dir, log error if it fails
cd "$base_dir" || { log "$logs_dir" "pipeline.log" "Error changing directory to $base_dir."; exit 1; }

# Copy the input file to the main results directory
cp "$input_file" "$sequencing_reads_main_abs"

# Iterate through selected modules
for module_script in "${selected_modules[@]}"; do
  module_name="${module_script%%.*}"

  case "$module_name" in
    "StrainCascade_LJA_assembly")
      # Run La Jolla Assembler module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_SPAdes_assembly")
      # Run SPAdes module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_Canu_assembly")
      # Run Canu module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_Flye_assembly")
      # Run Flye module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_assembly_evaluation1")
      # Run first StrainCascade_assembly_evaluation module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_MAC2_assembly_merging")
      # Run StrainCascade_MAC2.0_assembly_merging module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_assembly_evaluation2")
      # Run second StrainCascade_assembly_evaluation module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_Circlator_circularisation")
      # Run Circlator module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_assembly_evaluation3")
      # Run run final StrainCascade_assembly_evaluation module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_arrow_medaka_polishing")
      # Run BBMap_coverage module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$bam_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_NGMLR_BBMap_coverage")
      # Run BBMap_coverage module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_CheckM2_QC")
      # Run CheckM2 QC module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$databases_dir"
      ;;
    "StrainCascade_GTDB-Tk_taxonomy")
      # Run GTDB-Tk taxonomy classification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$databases_dir"
      ;;
    "StrainCascade_GTDB-Tk_de_novo_tree")
      # Run GTDB-Tk taxonomy classification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$databases_dir"
      ;;
    "StrainCascade_Bakta_annotation")
      # Run Bakta annotation module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$genome_annotation_main_abs" "$databases_dir" "$locus_tag"  "$force_overwrite"
      ;;
    "StrainCascade_Prokka_annotation")
      # Run Prokka annotation module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$genome_annotation_main_abs" "$databases_dir" "$locus_tag" "$force_overwrite"
      ;;
    "StrainCascade_MicrobeAnnotator_annotation")
      # Run MicrobeAnnotator annotation module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_annotation_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_PlasmidFinder_identification")
      # Run PlasmidFinder identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$genome_assembly_main_abs" "$databases_dir"
      ;;
    "StrainCascade_RGI_antimicrobial_resistance_identification")
      # Run RGI antimicrobial resistance identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_ResFinder_antimicrobial_resistance_identification")
      # Run ResFinder antimicrobial resistance identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$sequencing_type" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_dbCAN3_CAZymes_identification")
      # Run dbCAN CAZymes identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_IslandPath_genomic_islands_identification")
      # Run IslandPath genomic island identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$genome_annotation_main_abs" "$functional_analysis_main_abs"
      ;;
    "StrainCascade_VirSorter2_phage_identification")
      # Run VirSorter2 phage identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_DeepVirFinder_phage_identification")
      # Run VirSorter2 phage identification module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs"
      ;;
    "StrainCascade_results_summary")
      # Run results summary module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$output_dir" "$sample_name" "$miniconda_path" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$genome_annotation_main_abs"
      ;;
    "StrainCascade_result_processing")
      # Run StrainCascade_result_processing module
      execute_module "$module_script" "$module_name" "$script_dir" "$logs_dir" "$output_dir" "$sample_name" "$miniconda_path" "$main_results_dir_abs" "$results_integration_abs"
      ;;
    *)
      # Handle other modules if needed
      log "$logs_dir" "pipeline.log" "Warning: Unrecognized module $module_name"
      ;;
  esac
done

# Arrange the results based on the result type
case "$result_type" in
  "main")
    # Check if force overwrite is enabled
    if [ "$force_overwrite" = "yes" ]; then
      # Remove everything from output_dir except main_results_dir_abs
      find "$output_dir" -mindepth 1 -maxdepth 1 ! -name "$(basename "$main_results_dir_abs")" -exec rm -rf {} +
      log "$logs_dir" "pipeline.log" "Removed all contents from $output_dir except $(basename "$main_results_dir_abs")"
    else
      # Remove everything from output_dir except main_results_dir_abs, without force
      find "$output_dir" -mindepth 1 -maxdepth 1 ! -name "$(basename "$main_results_dir_abs")" -exec rm -rf {} +
      log "$logs_dir" "pipeline.log" "Removed all contents from $output_dir except $(basename "$main_results_dir_abs")"
    fi
    ;;
  "R")
    echo "Warning: R result type is still in development. Treating it as 'main' for now."
    # Check if force overwrite is enabled
    if [ "$force_overwrite" = "yes" ]; then
      # Remove everything from output_dir except main_results_dir_abs
      find "$output_dir" -mindepth 1 -maxdepth 1 ! -name "$(basename "$main_results_dir_abs")" -exec rm -rf {} +
      log "$logs_dir" "pipeline.log" "Removed all contents from $output_dir except $(basename "$main_results_dir_abs")"
    else
      # Remove everything from output_dir except main_results_dir_abs, without force
      find "$output_dir" -mindepth 1 -maxdepth 1 ! -name "$(basename "$main_results_dir_abs")" -exec rm -rf {} +
      log "$logs_dir" "pipeline.log" "Removed all contents from $output_dir except $(basename "$main_results_dir_abs")"
    fi
    ;;
  "all")
    # Do nothing
    ;;
  *)
    echo "Error: Invalid result type. Valid options are: all, main, R" >&2
    exit 1
    ;;
esac