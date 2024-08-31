#!/bin/bash

# StrainCascade_data_integration.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Check for the correct number of command line arguments
if [ "$#" -ne 9 ]; then
  echo "Usage: $0 <script_dir> <logs_dir> <utils_file> <apptainer_images_dir> <genome_assembly_main_abs> <taxonomic_classification_main_abs> <genome_annotation_main_abs> <functional_analysis_main_abs> <results_integration_abs>"
  exit 1
fi

script_dir=$1
logs_dir=$2
utils_file=$3
apptainer_images_dir=$4
genome_assembly_main_abs=$5
taxonomic_classification_main_abs=$6
genome_annotation_main_abs=$7
functional_analysis_main_abs=$8
results_integration_abs=$9

source "$utils_file"

