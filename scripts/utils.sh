#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# utils.sh

## straincascade info functions
usage() {
    echo "Usage: straincascade [OPTIONS]"
    echo
    echo "straincascade options:"
    echo "  -i, --input FILE/DIR          Input reads file or directory (mutually exclusive with -a)"
    echo "  -a, --assembly_input FILE/DIR  Input assembly file or directory (mutually exclusive with -i)"
    echo "  -ea, --external_assembly_dir DIR  External assembly directory"
    echo "  -o, --output_dir DIR          Output directory (default: current directory)"
    echo "  -s, --sequencing_type TYPE    Sequencing type (default: pacbio-hifi)"
    echo "                               Options: pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw | nano-corr | nano-hq"
    echo "  -t, --threads INT             Number of threads (default: 32)"
    echo "  -e, --execution_mode MODE     Execution mode (mutually exclusive with -b, default: standard)"
    echo "                               Options: minimal | efficient | standard | comprehensive | custom:modules"
    echo "  -b, --bundle NAME             Bundle name (mutually exclusive with -e)"
    echo "                               Options: assembly | annotation | functional | phage"
    echo "  -r, --result_type TYPE        Result type (default: main)"
    echo "                               Options: all | main | R"
    echo "  -sa, --selection_algorithm TYPE Assembly selection algorithm (default: contig)"
    echo "                               Options: contig | continuity"
    echo "  -l, --locus-tag STRING        Locus tag (default: automatic)"
    echo "  --deterministic | --heuristic Choose deterministic to optimise for reproducibility"
    echo "                               at cost of computation time (default: heuristic)"
    echo
    echo "Update options:"
    echo "  -us, --update-software        Update StrainCascade software"
    echo "  -uai, --update-images         Update all Apptainer images" 
    echo "  -udb, --update-databases      Update all databases"
    echo
    echo "General options:"
    echo "  -h, --help                    Show detailed help message"
    echo "  -v, --version                 Show version information"
    exit 1
}

help() {
    echo "Usage: straincascade [OPTIONS]"
    echo
    echo "Pipeline Options:"
    echo "  -i, --input FILE/DIR          Input reads file or directory (mutually exclusive with -a)"
    echo "                               Valid extensions: .fasta, .fna, .fa, .fastq, .fastq.gz, .bam"
    echo
    echo "  -a, --assembly_input FILE/DIR  Input assembly file or directory (mutually exclusive with -i)"
    echo "                               Valid extensions: .fasta, .fna, .fa"
    echo
    echo "  -ea, --external_assembly_dir DIR  External assembly directory"
    echo
    echo "  -o, --output_dir DIR          Output directory (default: current directory)"
    echo
    echo "  -s, --sequencing_type TYPE    Sequencing type (default: pacbio-hifi)"
    echo "                               Options: pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw | nano-corr | nano-hq"
    echo
    echo "  -t, --threads INT             Number of threads to use (default: 32)"
    echo
    echo "  -e, --execution_mode MODE     Execution mode (mutually exclusive with -b, default: standard)"
    echo "                               Options: minimal | efficient | standard | comprehensive | custom:modules"
    echo "                               For custom mode, specify modules as comma-separated list: custom:SC1,SC2,SC3"
    echo
    echo "                               Execution Modes:"
    echo "                               - minimal: Runs essential modules for quick analysis"
    echo "                                         Includes: SC3, SC6, SC13, SC14, SC16, SC27"
    echo "                               - efficient: Balanced analysis for thorough results"
    echo "                                         Includes: SC1, SC2, SC3, SC6, SC7, SC8, SC11, SC13, SC14, SC16, SC19, SC27"
    echo "                               - standard: Comprehensive analysis with most modules"
    echo "                                         Includes: SC1, SC2, SC3, SC4, SC5, SC6, SC7, SC8, SC9, SC10, SC11, SC12,"
    echo "                                                   SC13, SC14, SC16, SC17, SC18, SC19, SC20, SC21, SC22, SC23, SC24,"
    echo "                                                   SC26, SC27, SC28"
    echo "                               - comprehensive: Runs all available modules"
    echo
    echo "  -b, --bundle NAME             Bundle of modules to run (mutually exclusive with -e)"
    echo "                               Options: assembly | annotation | functional | phage"
    echo
    echo "                               Bundles:"
    echo "                               - assembly: Genome assembly-related modules"
    echo "                                         Includes: SC1, SC2, SC3, SC4, SC5, SC6, SC7, SC8, SC9, SC10, SC11, SC12, SC13"
    echo "                               - annotation: Genome annotation-related modules"
    echo "                                         Includes: SC16, SC17, SC18"
    echo "                               - functional: Functional analysis-related modules"
    echo "                                         Includes: SC14, SC20, SC21, SC22, SC23"
    echo "                               - phage: Phage identification-related modules"
    echo "                                         Includes: SC24, SC25, SC26"
    echo
    echo "  -r  TYPE                Result type to generate (default: main)"
    echo "                          Options: all | main | R"
    echo
    echo "  -sa ALGORITHM           Assembly selection algorithm (contig or continuity). Default: contig"
#    echo
#    echo "  -f  yes/no              Force overwrite existing results (default: yes)"
    echo
    echo "  -l  STRING              Locus tag for genome annotation (default: automatic)"
    echo "                          Must contain 3-12 alphanumeric uppercase characters and start with a letter"
    echo "  --deterministic | --heuristic    Reproducibility control:

                                    deterministic:
                                      - Focuses on reproducibility by controlling entropy sources and reducing variability.
                                      - Stabilizes entropy sources in containerized environments.
                                      - Sets fixed seeds in supported tools and limits multithreading.
                                      - May increase computation time significantly.

                                    heuristic (default):
                                      - Favoring performance over reproducibility." 
    echo
    echo "Update options:"
    echo "  -us, --update-software  Update StrainCascade software to the latest version"
    echo "  -uai, --update-images   Update all Apptainer images"
    echo "  -udb, --update-databases Update all databases"
    echo
    echo "General options:"
    echo "  -h, --help              Show this help message"
    echo "  -v, --version           Show version information"
    echo
    echo "Available modules:"
    echo "  SC1  Canu Correction and Trimming            - Corrects and trims long reads using Canu"
    echo "  SC2  LJA Assembly                            - Assembles genomes using LJA"
    echo "  SC3  SPAdes Assembly                         - Assembles genomes using SPAdes"
    echo "  SC4  Canu Assembly                           - Assembles genomes using Canu"
    echo "  SC5  Flye Assembly                           - Assembles genomes using Flye"
    echo "  SC6  Assembly Evaluation 1                   - Evaluates assembly quality (round 1)"
    echo "  SC7  MAC2 Assembly Merging                   - Merges multiple assemblies using MAC2"
    echo "  SC8  Assembly Evaluation 2                   - Evaluates assembly quality (round 2)"
    echo "  SC9  Circlator Circularisation               - Circularizes contigs using Circlator"
    echo "  SC10 Assembly Evaluation 3                   - Evaluates assembly quality (round 3)"
    echo "  SC11 Arrow Medaka Polishing                  - Polishes assemblies using Arrow/Medaka"
    echo "  SC12 NGMLR BBMap Coverage                    - Maps reads and assesses coverage"
    echo "  SC13 CheckM2 QC                              - Performs quality control using CheckM2"
    echo "  SC14 GTDB-Tk Taxonomy                        - Assigns taxonomy using GTDB-Tk"
    echo "  SC15 GTDB-Tk De Novo Tree                    - Constructs phylogenetic tree using GTDB-Tk"
    echo "  SC16 Bakta Annotation                        - Annotates genomes using Bakta"
    echo "  SC17 Prokka Annotation                       - Annotates genomes using Prokka"
    echo "  SC18 MicrobeAnnotator Annotation             - Annotates genomes using MicrobeAnnotator"
    echo "  SC19 PlasmidFinder ID                        - Identifies plasmids using PlasmidFinder"
    echo "  SC20 AMRFinderPlus Resistance ID             - Identifies AMR genes using AMRFinderPlus"
    echo "  SC21 ResFinder Resistance ID                 - Identifies AMR genes using ResFinder"
    echo "  SC22 dbCAN3 CAZyme ID                        - Identifies CAZymes using dbCAN3"
    echo "  SC23 IslandPath Genomic Islands ID           - Identifies genomic islands using IslandPath"
    echo "  SC24 VirSorter2 Phage ID                     - Identifies phages using VirSorter2"
    echo "  SC25 DeepVirFinder Phage ID                  - Identifies phages using DeepVirFinder"
    echo "  SC26 CRISPRCasFinder ID                      - Identifies CRISPRCas systems"
    echo "  SC27 ISEScan IS Elements ID                  - Identifies IS elements using ISEScan"
    echo "  SC28 Data Integration                        - Integrates results from various modules"
    echo
    echo "CAVEAT: Some modules may be skipped automatically depending on specific parameters,"
    echo "        such as the chosen sequencing type or input file format. StrainCascade will detect"
    echo "        and skip any modules that are incompatible with your configuration."
    exit 1
}


## straincascade function update functions
# Function to handle software update
handle_software_update() {
    # Store the current directory
    local current_dir=$(pwd)
    local parent_dir="$(dirname "$script_dir")"

    echo "Updating StrainCascade software..."
    # Change to script directory
    cd "$script_dir" || {
        echo "Error: Failed to change to script directory" >&2
        return 1
    }
    
    # Call the function to update scripts
    update_scripts
    
    # Remove git history
    if [ -n "${parent_dir}" ]; then
        if [ -d "${parent_dir}/.git" ]; then
            rm -rf "${parent_dir}/.git" 2>/dev/null || {
                echo "Warning: Failed to remove .git directory in ${parent_dir}" >&2
                true
            }
        fi
    else
        echo "Warning: parent_dir is not defined" >&2
    fi

    # Return to original directory
    cd "$current_dir" || {
        echo "Error: Failed to return to original directory" >&2
        return 1
    }
}

# Function to handle Apptainer images update
handle_images_update() {
    # Store the current directory
    local current_dir=$(pwd)
    
    # Change to script directory
    cd "$script_dir" || {
        echo "Error: Failed to change to script directory" >&2
        return 1
    }
    
    echo "Deleting all current Apptainer images..."
    # Remove all files in the Apptainer images directory
    rm -rf "$apptainer_images_dir"/*
    
    echo "Updating all Apptainer images..."
    # Call the function to pull new Apptainer images
    pull_apptainer_images
    
    # Return to original directory
    cd "$current_dir" || {
        echo "Error: Failed to return to original directory" >&2
        return 1
    }
}

# Function to handle database update
handle_database_update() {
    # Store the current directory
    local current_dir=$(pwd)
    
    # Change to script directory
    cd "$script_dir" || {
        echo "Error: Failed to change to script directory" >&2
        return 1
    }
    
    echo "Deleting all current databases..."
    # Remove all files in the databases directory
    rm -rf "$databases_directory"/*

    echo "Updating all databases..."
    # Define the location of the databases
    local db_location="$databases_directory"
    # Call the function to install all databases
    install_all_databases "$db_location"
    
    # Return to original directory
    cd "$current_dir" || {
        echo "Error: Failed to return to original directory" >&2
        return 1
    }
}

## Various utilities
# Function to create an output directory if it doesn't exist
create_directory() {
  if [[ ! -d "$1" ]]; then
    echo "Creating output directory: $1"
    mkdir -p "$1" || { echo "Error creating output directory."; exit 1; }
  fi
}

# Function to write progress logs
log() {
  local logs_dir="$1"
  local log_file="$2"
  local message="$3"
  echo "$(date): $message" >> "$logs_dir/$log_file"
}

# Function to execute a module
execute_module() {
    local module_script="$1"
    shift 1  # Shift the positional parameters to the left, so $@ contains the remaining arguments
    log "$LOGS_DIR" "pipeline.log" "Running $module_script module..."  # Log the start of the module

    bash "${SCRIPT_DIR}/${module_script}" "$SCRIPT_DIR" "$LOGS_DIR" "$@"  # Always passes LOGS_DIR and SCRIPT_DIR to the module script
    if [ $? -ne 0 ]; then
        log "$LOGS_DIR" "pipeline.log" "Error: $module_script module failed."  # Log an error if the module fails
        exit 1  # Exit the script with an error code
    fi
}


# Function to retrieve the absolute path of a file or directory (more portable than readlink -f or realpath)
get_absolute_path() {
    local path="$1"
    
    # First, try the simple method
    if [[ "$path" = /* ]]; then
        echo "$path"
        return 0
    elif [[ "$path" = ./* ]]; then
        echo "$PWD/${path#./}"
        return 0
    elif [[ "$path" != *..* ]]; then
        echo "$PWD/$path"
        return 0
    fi

    # If the simple method didn't work (e.g., path contains '..'), try realpath
    if command -v realpath >/dev/null 2>&1; then
        realpath "$path" 2>/dev/null && return 0
    fi

    # If realpath failed or isn't available, try readlink -f
    if command -v readlink >/dev/null 2>&1; then
        readlink -f "$path" 2>/dev/null && return 0
    fi

    # If all methods failed, return an error
    echo "Error: Unable to determine absolute path for '$path'" >&2
    return 1
}

# Function to find the analysis assembly file based on specific patterns
find_analysis_assembly_file() {
  local genome_assembly_main_abs="$1"
  local analysis_assembly_file=""
  local patterns=("_final" "best_ev3" "_circularised" "best_ev2" "_mac2.0_merged_assembly" "best_ev1")

  # Check for files matching each pattern in order
  for pattern in "${patterns[@]}"; do
    analysis_assembly_file=$(find "$genome_assembly_main_abs" -type f -name "*${pattern}*.fasta" | head -n 1)
    if [[ -n $analysis_assembly_file ]]; then
      echo "$analysis_assembly_file"
      return 0
    fi
  done

  # If no files match the patterns, default to the first .fasta file
  analysis_assembly_file=$(find "$genome_assembly_main_abs" -type f -name "*.fasta" | head -n 1)
  if [[ -n $analysis_assembly_file ]]; then
    echo "$analysis_assembly_file"
  else
    echo "Error: No suitable assembly file found." >&2
    return 1
  fi
}

find_apptainer_sif_file() {
    local dir="$1"
    local pattern="$2"
    readarray -t matching_files < <(find "$dir" -name "$pattern" -print)
    if [ ${#matching_files[@]} -eq 0 ]; then
        echo "No matching .sif files found in $dir. Continuing with the next module in the pipeline." >&2
        exit 0
    elif [ ${#matching_files[@]} -gt 1 ]; then
        echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}" >&2
    fi
    echo "${matching_files[0]}"
}