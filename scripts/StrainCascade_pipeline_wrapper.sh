#!/bin/bash

# StrainCascade_pipeline_wrapper.sh - Version 1.1.0
# Author: Sebastian Bruno Ulrich Jordi

# Wrapper script for StrainCascade pipeline

script_dir=$1
utils_file=$2
input_file=$3
output_directory=$4
databases_dir=$5
apptainer_images_dir=$6
sequencing_type=$7
input_type=$8
threads=$9
selected_modules=${10}
result_type=${11}
force_overwrite=${12}
locus_tag=${13}
version=${14}

## Define variables and create directories ##
selected_modules_string="$selected_modules"

# Convert the selected_modules string to an array
IFS=' ' read -r -a selected_modules <<< "$selected_modules"

source "$utils_file"

# Define the sample name variable
sample_name="$input_file"
for suffix in .fasta .fa .fastq .fastq.gz .fna .bam; do
  sample_name=$(basename "$sample_name" "$suffix")
done

# Define base_dir=$output_directory (to be used by the pipeline scripts)
base_dir="$output_directory"

# Create the output_dir (the output_dir is derived from the argument output_directory and will be the actual folder to contain results)
output_dir="${base_dir}/StrainCascade_output_${sample_name}"
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

## Handle inputfile
# Use execute_module to run the input file handler
input_handler_output=$(execute_module "StrainCascade_input_file_handler.sh" "$apptainer_images_dir" "$input_file" "$output_dir" "$sequencing_reads_main_abs" "$sequencing_reads_main_abs" "$input_type")

if [ $? -ne 0 ]; then
    log "$logs_dir" "pipeline.log" "Error: Input file handling failed."
    exit 1
fi

# Parse the output to get both input_file and bam_file | input_file is redefined as the original or processed input file copied to sequencing_reads_main_abs
IFS=$'\t' read -r original_input input_file bam_file <<< "$input_handler_output"

# Log the results
log "$logs_dir" "pipeline.log" "Original input: $original_input"
log "$logs_dir" "pipeline.log" "Updated input file: $input_file"
log "$logs_dir" "pipeline.log" "BAM file (if applicable): $bam_file"
log "$logs_dir" "pipeline.log" "Input type: $input_type"

## Start processing ##
# Log the attempt to change directory
log "$logs_dir" "pipeline.log" "Changing directory to: $base_dir"

# Change directory to base_dir, log error if it fails
cd "$base_dir" || { log "$logs_dir" "pipeline.log" "Error changing directory to $base_dir."; exit 1; }

# Iterate through selected modules (execute_module has passing of script_dir and logs_dir as first variables hard-coded)
for module_script in "${selected_modules[@]}"; do
  case "$module_script" in
    "StrainCascade_Canu_correct_trim.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run Canu correction and trimming module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_LJA_assembly.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run La Jolla Assembler module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_SPAdes_assembly.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run SPAdes module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_Canu_assembly.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run Canu module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_Flye_assembly.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run Flye module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_assembly_evaluation1.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run first StrainCascade_assembly_evaluation module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_MAC2_assembly_merging.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run StrainCascade_MAC2.0_assembly_merging module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_assembly_evaluation2.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run second StrainCascade_assembly_evaluation module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_Circlator_circularisation.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run Circlator module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_assembly_evaluation3.sh")
      # Run run final StrainCascade_assembly_evaluation module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_arrow_medaka_polishing.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run arrow_medaka_polishing module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$bam_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_NGMLR_BBMap_coverage.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run NGMLR_BBMap_coverage module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$sequencing_type" "$threads" "$genome_assembly_main_abs"
      ;;
    "StrainCascade_CheckM2_QC.sh")
      # Run CheckM2 QC module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$databases_dir"
      ;;
    "StrainCascade_GTDB-Tk_taxonomy.sh")
      # Run GTDB-Tk taxonomy classification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$databases_dir"
      ;;
    "StrainCascade_GTDB-Tk_de_novo_tree.sh")
      # Run GTDB-Tk de novo tree classification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$databases_dir"
      ;;
    "StrainCascade_Bakta_annotation.sh")
      # Run Bakta annotation module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$genome_annotation_main_abs" "$databases_dir" "$locus_tag"  "$force_overwrite"
      ;;
    "StrainCascade_Prokka_annotation.sh")
      # Run Prokka annotation module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$genome_annotation_main_abs" "$databases_dir" "$locus_tag" "$force_overwrite"
      ;;
    "StrainCascade_MicrobeAnnotator_annotation.sh")
      # Run MicrobeAnnotator annotation module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_annotation_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_PlasmidFinder_identification.sh")
      # Run PlasmidFinder identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$genome_assembly_main_abs" "$databases_dir"
      ;;
    "StrainCascade_RGI_antimicrobial_resistance_identification.sh")
      # Run RGI antimicrobial resistance identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_ResFinder_antimicrobial_resistance_identification.sh")
      # Run ResFinder antimicrobial resistance identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$sequencing_type" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_dbCAN3_CAZymes_identification.sh")
      # Run dbCAN CAZymes identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_IslandPath_genomic_islands_identification.sh")
      # Run IslandPath genomic island identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$genome_annotation_main_abs" "$functional_analysis_main_abs"
      ;;
    "StrainCascade_VirSorter2_phage_identification.sh")
      # Run VirSorter2 phage identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs" "$databases_dir"
      ;;
    "StrainCascade_DeepVirFinder_phage_identification.sh")
      if [ "$input_type" = "assembly" ]; then
        log "$logs_dir" "pipeline.log" "Skipping $module_script because an assembly is given as input when reads are needed"
        continue
      fi
      # Run DeepVirFinder phage identification module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$input_file" "$output_dir" "$sample_name" "$threads" "$genome_assembly_main_abs" "$functional_analysis_main_abs"
      ;;
    "StrainCascade_data_integration.sh")
      # Run StrainCascade data integration module
      execute_module "$module_script" "$utils_file" "$apptainer_images_dir" "$genome_assembly_main_abs" "$taxonomic_classification_main_abs" "$genome_annotation_main_abs" "$functional_analysis_main_abs" "$results_integration_abs"
      ;;
    *)
      # Handle other modules if needed
      log "$logs_dir" "pipeline.log" "Warning: Unrecognized module $module_script"
      ;;
  esac
done

# Create a run summary
execute_module "StrainCascade_run_summary.sh" "$utils_file" "$apptainer_images_dir" "$input_file" "$sample_name" "$sequencing_type" "$threads" "$databases_dir" "$main_results_dir_abs" "$selected_modules_string" "$version"

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