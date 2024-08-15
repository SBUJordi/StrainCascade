#!/bin/bash

# Utils script:
# utils.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

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
    local module_name="$2"
    shift 2  # Shift the positional parameters to the left, so $@ contains the remaining arguments
    log "$base_dir/pipeline_logs" "pipeline.log" "Running $module_name module..."
    local script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    bash "${script_dir}/${module_script}" "$@"  # Run the module script from the dir of the utils.sh (script_dir) and wait for it to finish
    if [ $? -ne 0 ]; then
        log "$base_dir/pipeline_logs" "pipeline.log" "Error: $module_name module failed."
        exit 1
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
