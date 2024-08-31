#!/bin/bash

# Utils script:
# utils.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# straincascade info functions
usage() {
    echo "Usage: $0 (-i input | -a assembly) [-o output_directory] [-s seq_type] [-t threads] [-r result_type] [-e execution_mode | -b bundle] [-f force_overwrite] [-d databases_directory] [-l locus_tag] [-h] [-v]"
    echo "  -i  input               Input file or directory containing sequencing files (mutually exclusive with -a); valid file extensions include: .fasta, .fna, .fa, .fastq, .fastq.gz, and .bam"
    echo "  -a  assembly            Input assembly file or directory containing assembly files (mutually exclusive with -i); valid file extensions include: .fasta, .fna, .fa"
    echo "  -o  output_directory    Output directory (default: current working directory)"
    echo "  -s  seq_type            Sequence type (default: pacbio-hifi). Options are: pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw | nano-corr | nano-hq"
    echo "  -t  threads             Number of threads (default: 32)"
    echo "  -e  execution_mode      Execution mode (mutually exclusive with -b). Options are: minimal | efficient | standard | comprehensive | custom"
    echo "                          For custom mode, indicate the modules to run; e.g., \"custom:SC1,SC2,SC3\""
    echo "  -b  bundle              Bundle of modules to run (mutually exclusive with -e). Options are: assembly | annotation | functional | phage"
    echo "  -r  result_type         Result type (default: main). Options are: all, main, R"
    echo "  -f  force_overwrite     Force overwrite existing results (default: yes). Options are: yes | no"
    echo "  -d  databases_directory Path to databases directory (default: StrainCascade/databases)"
    echo "  -l  locus_tag           Specify your own locus tag for genome annotation (default: automatically generated); must contain between 3 and 12 alphanumeric uppercase characters and start with a letter"
    echo "  -h                      Show detailed help message"
    echo "  -v                      Show version information"
    exit 1
}

help() {
    echo "Usage: $0 (-i input | -a assembly) [-o output_directory] [-s seq_type] [-t threads] [-r result_type] [-e execution_mode | -b bundle] [-f force_overwrite] [-d databases_directory] [-l locus_tag] [-h] [-v]"
    echo ""
    echo "Options:"
    echo "  -i  input               Specify the input file or directory containing sequencing files (mutually exclusive with -a);"
    echo "                          Valid file extensions: .fasta, .fna, .fa, .fastq, .fastq.gz, .bam"
    echo ""
    echo "  -a  assembly            Specify the input assembly file or directory containing assembly files (mutually exclusive with -i);"
    echo "                          Valid file extensions: .fasta, .fna, .fa"
    echo ""
    echo "  -o  output_directory    Define the output directory (default: current working directory)"
    echo ""
    echo "  -s  seq_type            Set the sequence type (default: pacbio-hifi)."
    echo "                          Options: pacbio-raw | pacbio-corr | pacbio-hifi | nano-raw | nano-corr | nano-hq"
    echo ""
    echo "  -t  threads             Number of threads to use (default: 32)"
    echo ""
    echo "  -e  execution_mode      Select the execution mode (mutually exclusive with -b)."
    echo "                          Options: minimal | efficient | standard | comprehensive | custom"
    echo "                          For custom mode, specify the modules to run; e.g., \"custom:SC1,SC2,SC3\""
    echo ""
    echo "                          Execution Modes:"
    echo "                          - minimal: Runs essential modules for quick analysis"
    echo "                            Includes: SC3, SC6, SC13, SC14, SC16, SC26"
    echo "                          - efficient: A balanced set of modules for thorough analysis"
    echo "                            Includes: SC1, SC3, SC2, SC6, SC7, SC8, SC11, SC13, SC14, SC16, SC19, SC26"
    echo "                          - standard: Comprehensive analysis with most modules"
    echo "                            Includes: SC1, SC2, SC3, SC4, SC5, SC6, SC7, SC8, SC9, SC10, SC11, SC12,"
    echo "                                     SC13, SC14, SC16, SC17, SC18, SC19, SC20, SC21, SC22, SC23, SC24,"
    echo "                                     SC25, SC26"
    echo "                          - comprehensive: Runs all available modules"
    echo ""
    echo "  -b  bundle              Specify a bundle of modules to run (mutually exclusive with -e)."
    echo "                          Options: assembly | annotation | functional | phage"
    echo ""
    echo "                          Bundles:"
    echo "                          - assembly: Genome assembly-related modules"
    echo "                            Includes: SC1, SC2, SC3, SC4, SC5, SC6, SC7, SC8, SC9, SC10, SC11, SC12, SC13"
    echo "                          - annotation: Genome annotation-related modules"
    echo "                            Includes: SC16, SC17, SC18"
    echo "                          - functional: Functional analysis-related modules"
    echo "                            Includes: SC14, SC20, SC21, SC22, SC23"
    echo "                          - phage: Phage identification-related modules"
    echo "                            Includes: SC24, SC25"
    echo ""
    echo "  -r  result_type         Specify the type of results to generate (default: main)."
    echo "                          Options: all | main | R"
    echo ""
    echo "  -f  force_overwrite     Force overwrite existing results (default: yes)."
    echo "                          Options: yes | no"
    echo ""
    echo "  -d  databases_directory Specify the path to the databases directory (default: StrainCascade/databases)"
    echo ""
    echo "  -l  locus_tag           Specify a custom locus tag for genome annotation (default: automatically generated);"
    echo "                          Must contain 3-12 alphanumeric uppercase characters and start with a letter"
    echo ""
    echo "  -h                      Show this help message and exit"
    echo ""
    echo "  -v                      Show version information and exit"
    echo ""
    echo "Available Modules:"
    echo "  SC1  Canu Correction and Trimming            - Corrects and trims long reads using Canu"
    echo "  SC2  LJA Assembly                            - Assembles genomes using LJA"
    echo "  SC3  SPAdes Assembly                         - Assembles genomes using SPAdes"
    echo "  SC4  Canu Assembly                           - Assembles genomes using Canu"
    echo "  SC5  Flye Assembly                           - Assembles genomes using Flye"
    echo "  SC6  Assembly Evaluation 1                   - Evaluates assembly quality and selects the best assembly (round 1)"
    echo "  SC7  MAC2 Assembly Merging                   - Merges multiple assemblies using MAC2"
    echo "  SC8  Assembly Evaluation 2                   - Evaluates assembly quality and selects the best assembly (round 2)"
    echo "  SC9  Circlator Circularisation               - Circularizes contigs using Circlator"
    echo "  SC10 Assembly Evaluation 3                   - Evaluates assembly quality and selects the final assembly (round 3)"
    echo "  SC11 Arrow Medaka Polishing                  - Polishes assemblies using Arrow or Medaka when applicable"
    echo "  SC12 NGMLR BBMap Coverage                    - Maps reads to assemblies and assesses coverage"
    echo "  SC13 CheckM2 QC                              - Performs quality control using CheckM2"
    echo "  SC14 GTDB-Tk Taxonomy                        - Assigns taxonomy using GTDB-Tk"
    echo "  SC15 GTDB-Tk De Novo Tree                    - Constructs a de novo phylogenetic tree using GTDB-Tk"
    echo "  SC16 Bakta Annotation                        - Annotates genomes using Bakta"
    echo "  SC17 Prokka Annotation                       - Annotates genomes using Prokka"
    echo "  SC18 MicrobeAnnotator Annotation             - Annotates genomes using MicrobeAnnotator"
    echo "  SC19 PlasmidFinder Identification            - Identifies plasmids using PlasmidFinder"
    echo "  SC20 RGI Antimicrobial Resistance ID         - Identifies antimicrobial resistance genes using RGI"
    echo "  SC21 ResFinder Antimicrobial Resistance ID   - Identifies antimicrobial resistance genes using ResFinder"
    echo "  SC22 dbCAN3 CAZymes Identification           - Identifies CAZymes using dbCAN3"
    echo "  SC23 IslandPath Genomic Islands ID           - Identifies genomic islands using IslandPath"
    echo "  SC24 VirSorter2 Phage Identification         - Identifies phages using VirSorter2"
    echo "  SC25 DeepVirFinder Phage Identification      - Identifies phages using DeepVirFinder"
    echo "  SC26 Data Integration                        - Integrates results from various modules"
    echo ""
    echo "  CAVEAT:                        Some modules may be skipped automatically depending on specific parameters, such as the chosen seq_type or input file format."
    echo "                                 StrainCascade will detect and skip any modules that are incompatible with your configuration."
    exit 1
}


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
    log "$logs_dir" "pipeline.log" "Running $module_script module..."  # Log the start of the module

    bash "${script_dir}/${module_script}" "$script_dir" "$logs_dir" "$@"  # Always pass logs_dir and script_dir to the module script
    if [ $? -ne 0 ]; then
        log "$logs_dir" "pipeline.log" "Error: $module_script module failed."  # Log an error if the module fails
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
        echo "No matching .sif files found in $dir. Continuing with the next script in the pipeline."
        exit 0
    elif [ ${#matching_files[@]} -gt 1 ]; then
        echo "Warning: Multiple matching .sif files found. Using the first match: ${matching_files[0]}"
    fi
    echo "${matching_files[0]}"
}