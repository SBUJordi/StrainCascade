#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_pipeline_wrapper.sh
# Description: Optimized wrapper script for StrainCascade

set -euo pipefail
export LC_ALL=C

# Input validation and usage
if [[ $# -ne 17 ]]; then
    cat << EOF
Usage: $0 <script_dir> <utils_file> <input_file> <external_assembly_dir> <output_directory> <databases_dir> 
          <apptainer_images_dir> <sequencing_type> <input_type> <threads> 
          <selected_modules> <result_type> <selection_algorithm> <force_overwrite> 
          <locus_tag> <reproducibility_mode> <version>
EOF
    exit 1
fi

# Constants from command line arguments
readonly SCRIPT_DIR="$1"
readonly UTILS_FILE="$2"
readonly INPUT_FILE="$3"
readonly EXTERNAL_ASSEMBLY_DIR="$4"
readonly OUTPUT_DIRECTORY="$5"
readonly DATABASES_DIR="$6"
readonly APPTAINER_IMAGES_DIR="$7"
readonly SEQUENCING_TYPE="$8"
readonly INPUT_TYPE="$9"
readonly THREADS="${10}"
readonly SELECTED_MODULES="${11}"
readonly RESULT_TYPE="${12}"
readonly SELECTION_ALGORITHM="${13}"
readonly FORCE_OVERWRITE="${14}"
readonly LOCUS_TAG="${15}"
readonly REPRODUCIBILITY_MODE="${16}"
readonly VERSION="${17}"

# Get the current timestamp
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# Source utility functions
source "$UTILS_FILE"

# Extract sample name from input file, removing common bioinformatics extensions
sample_name=$(basename "$INPUT_FILE")
for ext in .fasta .fa .fastq .fastq.gz .fna .bam; do
    sample_name="${sample_name%$ext}"
done

# Define directory structure
readonly BASE_DIR="$OUTPUT_DIRECTORY"
readonly OUTPUT_DIR="${BASE_DIR}/StrainCascade_output_${sample_name}"
readonly MAIN_RESULTS_DIR="${OUTPUT_DIR}/${sample_name}_StrainCascade_main_results"
readonly LOG_NAME="StrainCascade.log"

# Define subdirectories array for cleaner iteration
readonly SUBDIRS=(
    "00_sequencing_reads"
    "01_genome_assembly"
    "02_taxonomic_classification"
    "03_genome_annotation"
    "04_functional_analysis"
    "05_results_integration"
)

# Create directory structure
create_directory "$OUTPUT_DIR"
create_directory "$MAIN_RESULTS_DIR"

# Create all subdirectories and store absolute paths
declare -A dir_paths
for subdir in "${SUBDIRS[@]}"; do
    create_directory "${MAIN_RESULTS_DIR}/${subdir}"
    case "${subdir}" in
        "00_sequencing_reads") key="sequencing_reads" ;;
        "01_genome_assembly") key="genome_assembly" ;;
        "02_taxonomic_classification") key="taxonomic_classification" ;;
        "03_genome_annotation") key="genome_annotation" ;;
        "04_functional_analysis") key="functional_analysis" ;;
        "05_results_integration") key="results_integration" ;;
    esac
    dir_paths["${key}"]="$(get_absolute_path "${MAIN_RESULTS_DIR}/${subdir}")"
done

# Create logs directory
readonly LOGS_DIR="${MAIN_RESULTS_DIR}/pipeline_logs"
create_directory "$LOGS_DIR"

# Function to execute a module
execute_module() {
    local module_script="$1"
    shift 1  # Shift the positional parameters to the left, so $@ contains the remaining arguments
    log "$LOGS_DIR" "$LOG_NAME" "Running $module_script module..."  # Log the start of the module

    bash "${SCRIPT_DIR}/${module_script}" "$SCRIPT_DIR" "$LOGS_DIR" "$@"  # Always passes LOGS_DIR and SCRIPT_DIR to the module script
    if [ $? -ne 0 ]; then
        log "$LOGS_DIR" "$LOG_NAME" "Error: $module_script module failed."  # Log an error if the module fails
        exit 1  # Exit the script with an error code
    fi
}

# Handle input file
log "$LOGS_DIR" "$LOG_NAME" "Processing input file: $INPUT_FILE"
input_handler_output=$(execute_module "StrainCascade_input_file_handler.sh" \
    "$LOG_NAME" \
    "$APPTAINER_IMAGES_DIR" \
    "$INPUT_FILE" \
    "$OUTPUT_DIR" \
    "${dir_paths[sequencing_reads]}" \
    "${dir_paths[genome_assembly]}" \
    "$INPUT_TYPE")

if [[ $? -ne 0 ]]; then
    log "$LOGS_DIR" "$LOG_NAME" "Error: Input file handling failed"
    exit 1
fi

# Parse input handler output
IFS=$'\t' read -r processed_input bam_file <<< "$input_handler_output"

# Log input processing results
log "$LOGS_DIR" "$LOG_NAME" "Original input: $INPUT_FILE"
log "$LOGS_DIR" "$LOG_NAME" "Processed input: $processed_input"
log "$LOGS_DIR" "$LOG_NAME" "BAM file: ${bam_file:-None}"

# Change to base directory
cd "$BASE_DIR" || {
    log "$LOGS_DIR" "$LOG_NAME" "Error: Failed to change to directory $BASE_DIR"
    exit 1
}

# Module execution mapping
declare -A module_params=(
    ["StrainCascade_Canu_correct_trim.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]} $REPRODUCIBILITY_MODE"
    ["StrainCascade_LJA_assembly.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]} $REPRODUCIBILITY_MODE"
    ["StrainCascade_Flye_assembly.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]} $REPRODUCIBILITY_MODE"
    ["StrainCascade_Canu_assembly.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]} $REPRODUCIBILITY_MODE"
    ["StrainCascade_SPAdes_assembly.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]} $REPRODUCIBILITY_MODE"
    ["StrainCascade_assembly_evaluation1.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} $SELECTION_ALGORITHM"
    ["StrainCascade_MAC2_assembly_merging.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name ${dir_paths[genome_assembly]}"
    ["StrainCascade_assembly_evaluation2.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} $SELECTION_ALGORITHM"
    ["StrainCascade_Circlator_circularisation.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name ${dir_paths[genome_assembly]} ${dir_paths[results_integration]} $VERSION"
    ["StrainCascade_assembly_evaluation3.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[results_integration]} $VERSION $SELECTION_ALGORITHM"
    ["StrainCascade_arrow_medaka_polishing.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $bam_file $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]} $REPRODUCIBILITY_MODE"
    ["StrainCascade_NGMLR_BBMap_coverage.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $SEQUENCING_TYPE $THREADS ${dir_paths[genome_assembly]}"
    ["StrainCascade_CheckM2_QC.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_GTDB-Tk_taxonomy.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[taxonomic_classification]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_GTDB-Tk_de_novo_tree.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $EXTERNAL_ASSEMBLY_DIR $THREADS ${dir_paths[genome_assembly]} ${dir_paths[taxonomic_classification]} $DATABASES_DIR"
    ["StrainCascade_Bakta_annotation.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[taxonomic_classification]} ${dir_paths[genome_annotation]} ${dir_paths[results_integration]} $DATABASES_DIR $LOCUS_TAG $FORCE_OVERWRITE $VERSION $REPRODUCIBILITY_MODE"
    ["StrainCascade_Prokka_annotation.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[genome_annotation]} ${dir_paths[results_integration]} $LOCUS_TAG $FORCE_OVERWRITE $VERSION $REPRODUCIBILITY_MODE"
    ["StrainCascade_MicrobeAnnotator_annotation.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_annotation]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION $REPRODUCIBILITY_MODE"
    ["StrainCascade_PlasmidFinder_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name ${dir_paths[genome_assembly]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[taxonomic_classification]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_ResFinder_antimicrobial_resistance_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS $SEQUENCING_TYPE ${dir_paths[genome_assembly]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_dbCAN3_CAZymes_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_IslandPath_genomic_islands_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name ${dir_paths[genome_assembly]} ${dir_paths[genome_annotation]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $VERSION"
    ["StrainCascade_VirSorter2_phage_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $DATABASES_DIR $VERSION"
    ["StrainCascade_DeepVirFinder_phage_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $processed_input $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[functional_analysis]}"
    ["StrainCascade_CRISPRCasFinder_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $VERSION"
    ["StrainCascade_ISEScan_IS_elements_identification.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR $OUTPUT_DIR $sample_name $THREADS ${dir_paths[genome_assembly]} ${dir_paths[functional_analysis]} ${dir_paths[results_integration]} $VERSION"
    ["StrainCascade_data_integration.sh"]="$LOG_NAME $UTILS_FILE $APPTAINER_IMAGES_DIR ${dir_paths[genome_assembly]} ${dir_paths[results_integration]} $VERSION $sample_name"
)

# Process selected modules
IFS=' ' read -r -a modules <<< "$SELECTED_MODULES"
for module in "${modules[@]}"; do
    # Skip assembly-dependent modules if input is assembly
    if [[ "$INPUT_TYPE" == "assembly" ]] && [[ "$module" =~ (Canu|LJA|SPAdes|Flye|evaluation[12]|MAC2|Circlator|arrow_medaka|NGMLR|DeepVirFinder) ]]; then
        log "$LOGS_DIR" "$LOG_NAME" "Skipping $module: Assembly input provided"
        continue
    fi
    
    if [[ -n "${module_params[$module]:-}" ]]; then
        execute_module "$module" ${module_params[$module]}
    else
        log "$LOGS_DIR" "$LOG_NAME" "Warning: Unrecognized module $module"
    fi
done

# Get the current timestamp as the stop time
STOP_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# Generate run summary
execute_module "StrainCascade_run_summary.sh" \
    "$LOG_NAME" \
    "$UTILS_FILE" \
    "$APPTAINER_IMAGES_DIR" \
    "$INPUT_FILE" \
    "$INPUT_TYPE" \
    "$sample_name" \
    "$SEQUENCING_TYPE" \
    "$THREADS" \
    "$DATABASES_DIR" \
    "$MAIN_RESULTS_DIR" \
    "${dir_paths[genome_assembly]}" \
    "${dir_paths[results_integration]}" \
    "$SELECTED_MODULES" \
    "$VERSION" \
    "$SELECTION_ALGORITHM" \
    "$REPRODUCIBILITY_MODE" \
    "$START_TIME" \
    "$STOP_TIME"

# Handle result type cleanup
if [[ "$RESULT_TYPE" != "all" ]]; then
    if [[ "$FORCE_OVERWRITE" == "yes" || "$RESULT_TYPE" == "main" || "$RESULT_TYPE" == "R" ]]; then
        find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 ! -name "$(basename "$MAIN_RESULTS_DIR")" -exec rm -rf {} +
        log "$LOGS_DIR" "$LOG_NAME" "Cleaned up output directory, keeping only main results"
    fi
fi

exit 0