#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# utils.sh

# Source theme/branding (colors, logo, display functions)
_utils_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${_utils_dir}/straincascade_theme.sh" ]]; then
    source "${_utils_dir}/straincascade_theme.sh"
fi

## straincascade info functions
usage() {
    print_header
    echo -e "${C_BOLD}Usage:${C_RESET} straincascade ${C_BLUE}[OPTIONS]${C_RESET}"
    
    print_section "Input Options"
    print_option "-i, --input FILE" "Input reads file or directory"
    print_option "-a, --assembly_input FILE" "Input assembly file or directory"
    print_option "-sr1, --short-reads-r1 FILE" "Illumina forward reads (R1)"
    print_option "-sr2, --short-reads-r2 FILE" "Illumina reverse reads (R2)"
    print_option "-ea, --external_assembly_dir" "External assembly directory"
    
    print_section "Output Options"
    print_option "-o, --output_dir DIR" "Output directory (default: .)"
    print_option "-r, --result_type TYPE" "all | main | R (default: main)"
    
    print_section "Pipeline Configuration"
    print_option "-s, --sequencing_type TYPE" "Sequencing platform type"
    print_sub "pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw | nano-corr | nano-hq"
    print_option "-e, --execution_mode MODE" "Execution mode (default: standard)"
    print_sub "minimal | efficient | standard | comprehensive | custom:SC1,SC2,..."
    print_option "-b, --bundle NAME" "Module bundle (mutually exclusive with -e)"
    print_sub "assembly | annotation | functional | phage"
    print_option "-t, --threads INT" "Number of threads (default: 32)"
    print_option "-sa, --selection_algorithm" "contig | continuity (default: contig)"
    print_option "-l, --locus-tag STRING" "Locus tag for annotation"
    print_option "--deterministic" "Reproducible mode (slower)"
    print_option "--heuristic" "Performance mode (default)"
    
    print_section "Maintenance"
    print_option "-us, --update-software" "Update StrainCascade"
    print_option "-uai, --update-images" "Update Apptainer images"
    print_option "-udb, --update-databases" "Update all databases"
    
    print_section "Help"
    print_option "-h, --help" "Show detailed help"
    print_option "-v, --version" "Show version"
    
    echo
    exit 1
}

help() {
    print_header
    echo -e "${C_BOLD}Usage:${C_RESET} straincascade ${C_BLUE}[OPTIONS]${C_RESET}"
    
    print_section "Input Options"
    print_option "-i, --input FILE" "Input reads file or directory"
    print_sub "Valid: .fasta, .fna, .fa, .fastq, .fastq.gz, .bam"
    print_sub "Mutually exclusive with -a"
    echo
    print_option "-a, --assembly_input FILE" "Input assembly file or directory"
    print_sub "Valid: .fasta, .fna, .fa"
    print_sub "Mutually exclusive with -i"
    echo
    print_option "-sr1, --short-reads-r1 FILE" "Illumina forward reads (R1)"
    print_option "-sr2, --short-reads-r2 FILE" "Illumina reverse reads (R2)"
    print_sub "Valid: .fastq, .fq, .fastq.gz, .fq.gz"
    print_sub "Use together for hybrid assembly (Unicycler)"
    echo
    print_option "-ea, --external_assembly_dir" "External assembly directory"
    
    print_section "Output Options"
    print_option "-o, --output_dir DIR" "Output directory (default: current directory)"
    print_option "-r, --result_type TYPE" "Result type to generate"
    print_sub "all | main | R"
    
    print_section "Sequencing Configuration"
    print_option "-s, --sequencing_type TYPE" "Sequencing platform type"
    print_sub "PacBio:   pacbio-raw | pacbio-corr | pacbio-hifi"
    print_sub "Nanopore: nano-raw | nano-corr | nano-hq"
    
    print_section "Execution Modes"
    print_option "-e, --execution_mode MODE" "Pipeline execution mode"
    echo -e "\n  ${C_WHITE}minimal${C_RESET}       Essential modules for quick analysis"
    print_sub "SC3, SC7, SC14, SC15, SC17, SC30"
    echo -e "\n  ${C_WHITE}efficient${C_RESET}     Balanced analysis for thorough results"
    print_sub "SC1-3, SC7-12, SC14-15, SC17, SC21, SC30"
    echo -e "\n  ${C_WHITE}standard${C_RESET}      Comprehensive analysis (default)"
    print_sub "SC1-6, SC7-15, SC17-30"
    echo -e "\n  ${C_WHITE}comprehensive${C_RESET} All available modules"
    echo -e "\n  ${C_WHITE}custom:...${C_RESET}    Custom module list (e.g., custom:SC1,SC2,SC3)"
    
    print_section "Module Bundles"
    print_option "-b, --bundle NAME" "Predefined module bundle (mutually exclusive with -e)"
    echo -e "\n  ${C_WHITE}assembly${C_RESET}      Genome assembly pipeline"
    print_sub "SC1-6, SC7-14"
    echo -e "\n  ${C_WHITE}annotation${C_RESET}    Genome annotation"
    print_sub "SC17-20"
    echo -e "\n  ${C_WHITE}functional${C_RESET}    Functional analysis"
    print_sub "SC15, SC22-25"
    echo -e "\n  ${C_WHITE}phage${C_RESET}         Phage identification"
    print_sub "SC26-29"
    
    print_section "Processing Options"
    print_option "-t, --threads INT" "Number of threads (default: 32)"
    print_option "-sa, --selection_algorithm" "Assembly selection: contig | continuity"
    print_option "-l, --locus-tag STRING" "Locus tag for annotation"
    print_sub "3-12 alphanumeric uppercase, starts with letter"
    
    print_section "Reproducibility"
    print_option "--deterministic" "Reproducible mode"
    print_sub "Controls entropy, fixed seeds, single-threaded"
    print_sub "${C_YELLOW}Significantly increases computation time${C_RESET}"
    print_option "--heuristic" "Performance mode (default)"
    print_sub "Optimized for speed over reproducibility"
    
    print_section "Maintenance"
    print_option "-us, --update-software" "Update StrainCascade to latest version"
    print_option "-uai, --update-images" "Update all Apptainer images"
    print_option "-udb, --update-databases" "Update all databases"
    
    print_section "General"
    print_option "-h, --help" "Show this help message"
    print_option "-v, --version" "Show version information"
    
    print_section "Available Modules"
    echo -e "  ${C_DIM}─────────────────────────────────────────────────────────────────${C_RESET}"
    echo -e "  ${C_BOLD}${C_WHITE}Assembly${C_RESET}"
    print_module "SC1" "Canu Correction/Trimming" "Long read correction"
    print_module "SC2" "LJA Assembly" "La Jolla Assembler"
    print_module "SC3" "SPAdes Assembly" "Multi-algorithm assembler"
    print_module "SC4" "Canu Assembly" "Oxford/PacBio assembler"
    print_module "SC5" "Flye Assembly" "Long read assembler"
    print_module "SC6" "Unicycler Assembly" "Hybrid/long-read assembly"
    print_module "SC7" "Assembly Evaluation 1" "Quality metrics (round 1)"
    print_module "SC8" "MAC2 Merging" "Assembly merging"
    print_module "SC9" "Assembly Evaluation 2" "Quality metrics (round 2)"
    print_module "SC10" "Circlator" "Contig circularization"
    print_module "SC11" "Assembly Evaluation 3" "Quality metrics (round 3)"
    print_module "SC12" "Assembly Polishing" "Error correction"
    print_module "SC13" "minimap2/BBMap Coverage" "Read mapping & coverage"
    print_module "SC14" "CheckM2 QC" "Quality control"
    echo -e "\n  ${C_BOLD}${C_WHITE}Taxonomy & Annotation${C_RESET}"
    print_module "SC15" "GTDB-Tk Taxonomy" "Taxonomic classification"
    print_module "SC16" "GTDB-Tk Tree" "Phylogenetic tree"
    print_module "SC17" "Bakta Annotation" "Genome annotation"
    print_module "SC18" "Prokka Annotation" "Prokaryotic annotation"
    print_module "SC19" "DeepFRI" "Deep learning function prediction"
    print_module "SC20" "MicrobeAnnotator" "Metabolic annotation"
    echo -e "\n  ${C_BOLD}${C_WHITE}Functional Analysis${C_RESET}"
    print_module "SC21" "PlasmidFinder" "Plasmid identification"
    print_module "SC22" "AMRFinderPlus" "AMR gene detection"
    print_module "SC23" "ResFinder" "Resistance genes"
    print_module "SC24" "dbCAN3" "CAZyme identification"
    print_module "SC25" "IslandPath" "Genomic islands"
    echo -e "\n  ${C_BOLD}${C_WHITE}Mobile Elements & Phages${C_RESET}"
    print_module "SC26" "VirSorter2" "Phage detection"
    print_module "SC27" "geNomad" "Phage & plasmid ID"
    print_module "SC28" "CRISPRCasFinder" "CRISPR-Cas systems"
    print_module "SC29" "ISEScan" "IS element detection"
    echo -e "\n  ${C_BOLD}${C_WHITE}Integration${C_RESET}"
    print_module "SC30" "Data Integration" "Result consolidation"
    echo -e "  ${C_DIM}─────────────────────────────────────────────────────────────────${C_RESET}"
    
    echo -e "\n${C_YELLOW}Note:${C_RESET} Some modules may be skipped based on sequencing type or input format."
    
    print_section "Examples"
    print_example "Long-read assembly" "straincascade -i longreads.fastq.gz -o output/ -e standard"
    print_example "Hybrid assembly with Illumina reads" "straincascade -i longreads.fastq.gz -sr1 R1.fastq.gz -sr2 R2.fastq.gz -o output/"
    print_example "Custom module selection" "straincascade -i reads.fastq.gz -o output/ -e 'custom:SC6,SC7,SC14'"
    
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
# Returns empty string if input is empty
get_absolute_path() {
    local path="$1"
    
    # Return empty string for empty input
    if [[ -z "$path" ]]; then
        echo ""
        return 0
    fi
    
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
  # Priority order: final assembly, then progressively earlier stages
  # Note: _final pattern matches both old (*_best_ev3_final) and new (*_assembly_*_final) naming
  local patterns=("_final" "_circularised" "_mac2.0_merged_assembly")

  # Check for files matching each pattern in order
  for pattern in "${patterns[@]}"; do
    analysis_assembly_file=$(find "$genome_assembly_main_abs" -maxdepth 1 -type f -name "*${pattern}*.fasta" | sort | head -n 1)
    if [[ -n $analysis_assembly_file ]]; then
      echo "$analysis_assembly_file"
      return 0
    fi
  done

  # If no files match the patterns, default to the first .fasta file in the main dir
  analysis_assembly_file=$(find "$genome_assembly_main_abs" -maxdepth 1 -type f -name "*.fasta" | sort | head -n 1)
  if [[ -n $analysis_assembly_file ]]; then
    echo "$analysis_assembly_file"
  else
    echo "Error: No suitable assembly file found." >&2
    return 1
  fi
}

# Function to find Apptainer SIF file
# Prefers _latest.sif files, falls back to pattern matching for backwards compatibility
# Handles both new naming (image_latest.sif) and old naming (image:v4.sif, image_v4.sif)
find_apptainer_sif_file() {
    local dir="$1"
    local pattern="$2"
    
    # Extract base name from pattern (e.g., 'straincascade_genome_assembly*.sif' -> 'straincascade_genome_assembly')
    local base_name="${pattern%%\**}"
    base_name="${base_name%%:*}"  # Remove any :version suffix
    base_name="${base_name%.sif}" # Remove .sif if present
    
    # First, try to find the _latest.sif file directly (new naming convention)
    local latest_file="$dir/${base_name}_latest.sif"
    if [ -f "$latest_file" ]; then
        echo "$latest_file"
        return 0
    fi
    
    # Fall back to pattern matching for backwards compatibility
    # This handles old naming like :v4.sif or _v4.sif
    readarray -t matching_files < <(find "$dir" -maxdepth 1 -name "$pattern" -type f -print 2>/dev/null)
    if [ ${#matching_files[@]} -eq 0 ]; then
        echo "No matching .sif files found for pattern '$pattern' in $dir. Skipping module." >&2
        exit 0
    elif [ ${#matching_files[@]} -gt 1 ]; then
        # If multiple files found, prioritize any with "_latest" in the name
        for file in "${matching_files[@]}"; do
            if [[ "$file" == *"_latest"* ]]; then
                echo "Warning: Multiple matching .sif files found. Prioritizing: $file" >&2
                echo "$file"
                return 0
            fi
        done
        echo "Warning: Multiple matching .sif files found. Using: ${matching_files[0]}" >&2
    fi
    echo "${matching_files[0]}"
}