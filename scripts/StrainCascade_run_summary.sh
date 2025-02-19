#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_run_summary.sh
# Description: Generates detailed run summary and documentation for StrainCascade pipeline in both markdown and secure PDF formats

# set -euo pipefail

# Input parameters with descriptive names
readonly SCRIPT_DIR="$1"
readonly LOGS_DIR="$2"
readonly LOG_NAME="$3"
readonly UTILS_FILE="$4"
readonly APPTAINER_IMAGES_DIR="$5"
readonly INPUT_FILE="$6"
readonly INPUT_TYPE="$7"
readonly SAMPLE_NAME="$8"
readonly SEQUENCING_TYPE="$9"
readonly THREADS="${10}"
readonly DATABASES_DIR="${11}"
readonly MAIN_RESULTS_DIR="${12}"
readonly GENOME_ASSEMBLY_DIR="${13}"
readonly RESULTS_INTEGRATION_DIR="${14}"
readonly SELECTED_MODULES="${15}"
readonly VERSION="${16}"
readonly SELECTION_ALGORITHM="${17}"
readonly REPRODUCIBILITY_MODE="${18}"
readonly START_TIME="${19}"
readonly STOP_TIME="${20}"

# Source utility functions early to ensure availability
source "$UTILS_FILE"

# System information collection
readonly OS_INFO="$(uname -s)"
readonly OS_VERSION="$(cat /etc/os-release | grep "^PRETTY_NAME" | cut -d '=' -f 2 | tr -d '"')"
readonly KERNEL_VERSION="$(uname -r)"
readonly ARCHITECTURE="$(uname -m)"
readonly CPU_MODEL="$(grep "model name" /proc/cpuinfo | head -n1 | cut -d ':' -f2 | sed 's/^[ \t]*//')"
readonly TOTAL_MEMORY="$(free -h | awk '/^Mem:/{print $2}')"
readonly AVAILABLE_MEMORY="$(free -h | awk '/^Mem:/{print $7}')"

# Runtime calculation
start_seconds=$(date -d "$START_TIME" +%s 2>/dev/null || echo 0)
stop_seconds=$(date -d "$STOP_TIME" +%s 2>/dev/null || echo 0)

# Prevent negative durations
if [ "$stop_seconds" -lt "$start_seconds" ]; then
    duration=0
else
    duration=$((stop_seconds - start_seconds))
fi

hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

readonly RUN_DURATION_STRING=$(printf "%dh %02dm %02ds" $hours $minutes $seconds)

# Find assembly file
analysis_assembly_file=$(find_analysis_assembly_file "$GENOME_ASSEMBLY_DIR")
Rdata_file=$(find "${RESULTS_INTEGRATION_DIR}" -type f -name "StrainCascade_Results*.RData" | head -n 1)
HTML_file=$(find "${RESULTS_INTEGRATION_DIR}" -type f -name "StrainCascade_Analysis*.html" | head -n 1)


# Find all Apptainer images
readonly python_3_12_4_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'python_3.12.4*.sif')
echo "Found Python SIF: $python_3_12_4_sif"

readonly r_4_4_1_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'r_4.4.1*.sif')
echo "Found R SIF: $r_4_4_1_sif"

readonly straincascade_document_processing_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_document_processing*.sif')
echo "Found Document Processing SIF: $straincascade_document_processing_sif"

readonly straincascade_genome_assembly_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_genome_assembly*.sif')
echo "Found Genome Assembly SIF: $straincascade_genome_assembly_sif"

readonly straincascade_lja_genome_assembly_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_lja_genome_assembly*.sif')
echo "Found LJA Genome Assembly SIF: $straincascade_lja_genome_assembly_sif"

readonly straincascade_assembly_qc_refinement_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_assembly_qc_refinement*.sif')
echo "Found Assembly QC Refinement SIF: $straincascade_assembly_qc_refinement_sif"

readonly straincascade_genome_annotation_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_genome_annotation*.sif')
echo "Found Genome Annotation SIF: $straincascade_genome_annotation_sif"

readonly straincascade_taxonomic_functional_analysis_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_taxonomic_functional_analysis*.sif')
echo "Found Taxonomic Functional Analysis SIF: $straincascade_taxonomic_functional_analysis_sif"

readonly straincascade_crisprcas_phage_is_elements_sif=$(find_apptainer_sif_file "$APPTAINER_IMAGES_DIR" 'straincascade_crisprcas_phage_is_elements*.sif')
echo "Found CRISPR-Cas Phage IS Elements SIF: $straincascade_crisprcas_phage_is_elements_sif"

# Version collection
echo "Collecting tool version information..."

# Core tool versions
PYTHON_VERSION=$(apptainer exec "$python_3_12_4_sif" bash -c "python --version" 2>&1 | grep -E "Python" | awk '{print $2}')
echo "Python version: $PYTHON_VERSION"

R_VERSION=$(apptainer exec "$r_4_4_1_sif" bash -c "R --version" 2>&1 | grep -E "R version" | awk '{print $3}')
echo "R version: $R_VERSION"

# Assembly tools
LJA_VERSION="0.2"
echo "LJA version: $LJA_VERSION"

SPADES_VERSION=$(apptainer exec "$straincascade_genome_assembly_sif" bash -c "source activate spades_env && spades.py --version" 2>&1 | grep -E "SPAdes" | awk '{print $4}' | sed 's/^[a-zA-Z]//')
echo "SPAdes version: $SPADES_VERSION"

CANU_VERSION=$(apptainer exec "$straincascade_genome_assembly_sif" bash -c "source activate genome_assembly_env && canu --version" 2>&1 | awk '{print $2}' | sed 's/^[a-zA-Z]//')
echo "Canu version: $CANU_VERSION"

FLYE_VERSION=$(apptainer exec "$straincascade_genome_assembly_sif" bash -c "source activate genome_assembly_env && flye --version" 2>&1)
echo "Flye version: $FLYE_VERSION"


# QC tools
MAC2_VERSION="2.1 (with MUMmer 4.0.0rc1)"
echo "MAC2 version: $MAC2_VERSION"

CIRCLATOR_VERSION=$(apptainer exec "$straincascade_assembly_qc_refinement_sif" bash -c "source activate circlator_env && circlator version" 2>&1 | sed 's/^[a-zA-Z]//')
echo "Circlator version: $CIRCLATOR_VERSION"

ARROW_VERSION=$(apptainer exec "$straincascade_assembly_qc_refinement_sif" bash -c "source activate genomicconsensus_env && arrow --version" 2>&1 | sed 's/^[a-zA-Z]//')
echo "Arrow version: $ARROW_VERSION"

MEDAKA_VERSION="1.11.3"
echo "Medaka version: $MEDAKA_VERSION"

NGMLR_VERSION=$(apptainer exec "$straincascade_assembly_qc_refinement_sif" bash -c "source activate bbmap_ngmlr_env && ngmlr" 2>&1 | head -n 1 | awk '{print $2}' | sed 's/^[a-zA-Z]//')
echo "NGMLR version: $NGMLR_VERSION"

BBMAP_VERSION=$(apptainer exec "$straincascade_assembly_qc_refinement_sif" bash -c "source activate bbmap_ngmlr_env && bbmap.sh -v" 2>&1 | sed -n '3p' | awk '{print $2}' | sed 's/^[a-zA-Z]//')
echo "BBMap version: $BBMAP_VERSION"

QUAST_VERSION=$(apptainer exec "$straincascade_assembly_qc_refinement_sif" bash -c "source activate quast_env && quast.py --version" 2>&1 | grep -E "QUAST" | awk '{print $2}' | sed 's/^[a-zA-Z]//')
echo "QUAST version: $QUAST_VERSION"

CHECKM2_VERSION=$(apptainer exec "$straincascade_assembly_qc_refinement_sif" bash -c "source activate checkm2_env && checkm2 --version" 2>&1 | sed 's/^[a-zA-Z]//')
echo "CheckM2 version: $CHECKM2_VERSION"


# Genome annotation tools
BAKTA_VERSION=$(apptainer exec "$straincascade_genome_annotation_sif" bash -c "source activate bakta_env && bakta --version" 2>&1 | sed -n 's/^bakta \([0-9.]*\).*/\1/p')
echo "Bakta version: $BAKTA_VERSION"

PROKKA_VERSION=$(apptainer exec "$straincascade_genome_annotation_sif" bash -c "source activate prokka_env && prokka --version" 2>&1 | awk '{print $2}' | sed 's/^[a-zA-Z]//')
echo "Prokka version: $PROKKA_VERSION"

MICROBEANNOTATOR_VERSION=$(apptainer exec "$straincascade_genome_annotation_sif" bash -c "source activate microbeannotator_env && microbeannotator --version" 2>&1 | awk '{print $2}' | sed 's/^[a-zA-Z]//')
echo "MicrobeAnnotator version: $MICROBEANNOTATOR_VERSION"


# Taxonomy and functional analysis tools
GTDBTK_VERSION=$(apptainer exec "$straincascade_taxonomic_functional_analysis_sif" bash -c "source activate gtdbtk_env && gtdbtk -v" 2>&1 | awk '{print $3}' | sed 's/^[a-zA-Z]//')
echo "GTDB-Tk version: $GTDBTK_VERSION"

PLASMIDFINDER_VERSION="2.1.6"
echo "PlasmidFinder version: $PLASMIDFINDER_VERSION"

AMRFINDERPLUS_VERSION=$(apptainer exec "$straincascade_taxonomic_functional_analysis_sif" bash -c "source activate amrfinderplus_env && amrfinder --version" 2>&1 | sed 's/^[a-zA-Z]//')
echo "AMRFinderPlus version: $AMRFINDERPLUS_VERSION"

RESFINDER_VERSION=$(apptainer exec "$straincascade_taxonomic_functional_analysis_sif" bash -c "source activate resfinder_env && run_resfinder.py --version" 2>&1 | sed 's/^[a-zA-Z]//')
echo "ResFinder version: $RESFINDER_VERSION"

DBCAN_VERSION="4.1.4"
echo "dbCAN version: $DBCAN_VERSION"

ISLANDPATH_VERSION="1.0.6"
echo "IslandPath version: $ISLANDPATH_VERSION"

VIRSORTER2_VERSION=$(apptainer exec "$straincascade_crisprcas_phage_is_elements_sif" bash -c "source activate virsorter2_env && virsorter run -h" 2>&1 | head -n 1 | awk '{print $5}' | sed 's/^[a-zA-Z]//')
echo "VirSorter2 version: $VIRSORTER2_VERSION"

DEEPVIRFINDER_VERSION="1.0 (2020.11.21)"
echo "DeepVirFinder version: $DEEPVIRFINDER_VERSION"

CRISPRCASFINDER_VERSION=$(apptainer exec "$straincascade_crisprcas_phage_is_elements_sif" bash -c "source activate crisprcasfinder && /opt/CRISPRCasFinder/CRISPRCasFinder.pl -v" 2>&1 | sed -n '3p' | awk '{print $5}' | sed 's/^[a-zA-Z]//')
CRISPRCASFINDER_VERSION="${CRISPRCASFINDER_VERSION%,}"
echo "CRISPRCasFinder version: $CRISPRCASFINDER_VERSION"

ISESCAN_VERSION=$(apptainer exec "$straincascade_crisprcas_phage_is_elements_sif" bash -c "source activate isescan_env && isescan.py --version" 2>&1 | awk '{print $2}')
echo "ISEScan version: $ISESCAN_VERSION"

# Module definitions with versioning
declare -A MODULE_INFO
MODULE_INFO=(
    ["SC1"]="Canu Correction and Trimming (StrainCascade_Canu_correct_trim.sh). Using Canu v$CANU_VERSION"
    ["SC2"]="LJA Assembly (StrainCascade_LJA_assembly.sh). Using LJA v$LJA_VERSION"
    ["SC3"]="SPAdes Assembly (StrainCascade_SPAdes_assembly.sh). Using SPAdes v$SPADES_VERSION"
    ["SC4"]="Canu Assembly (StrainCascade_Canu_assembly.sh). Using Canu v$CANU_VERSION"
    ["SC5"]="Flye Assembly (StrainCascade_Flye_assembly.sh). Using Flye v$FLYE_VERSION"
    ["SC6"]="Assembly Evaluation 1 (StrainCascade_assembly_evaluation1.sh). Using QUAST v$QUAST_VERSION for assembly assessment and StrainCascade's assembly selection algorithm for assembly selection with python v$PYTHON_VERSION"
    ["SC7"]="MAC2.0 Assembly Merging (StrainCascade_MAC2_assembly_merging.sh). Using MAC.2 v$MAC2_VERSION"
    ["SC8"]="Assembly Evaluation 2 (StrainCascade_assembly_evaluation2.sh). Using QUAST v$QUAST_VERSION for assembly assessment and StrainCascade's assembly selection algorithm for assembly selection with python v$PYTHON_VERSION"
    ["SC9"]="Circlator Circularisation (StrainCascade_Circlator_circularisation.sh). Using Circlator v$CIRCLATOR_VERSION"
    ["SC10"]="Assembly Evaluation 3 (StrainCascade_assembly_evaluation3.sh). Using QUAST v$QUAST_VERSION for assembly assessment and StrainCascade's assembly selection algorithm for assembly selection with python v$PYTHON_VERSION"
    ["SC11"]="Arrow / Medaka Polishing (StrainCascade_arrow_medaka_polishing.sh). Using either Arrow v$ARROW_VERSION or Medaka v$MEDAKA_VERSION"
    ["SC12"]="NGMLR BBMap Coverage (StrainCascade_NGMLR_BBMap_coverage.sh). Using either NGMLR v$NGMLR_VERSION or BBMap v$BBMAP_VERSION"
    ["SC13"]="CheckM2 QC (StrainCascade_CheckM2_QC.sh). Using CheckM2 v$CHECKM2_VERSION"
    ["SC14"]="GTDB-Tk Taxonomy (StrainCascade_GTDB-Tk_taxonomy.sh). Using GTDB-Tk v$GTDBTK_VERSION"
    ["SC15"]="GTDB-Tk De Novo Tree (StrainCascade_GTDB-Tk_de_novo_tree.sh). Using GTDB-Tk v$GTDBTK_VERSION"
    ["SC16"]="Bakta Annotation (StrainCascade_Bakta_annotation.sh). Using Bakta v$BAKTA_VERSION"
    ["SC17"]="Prokka Annotation (StrainCascade_Prokka_annotation.sh). Using Prokka v$PROKKA_VERSION"
    ["SC18"]="MicrobeAnnotator Annotation (StrainCascade_MicrobeAnnotator_annotation.sh). Using MicrobeAnnotator v$MICROBEANNOTATOR_VERSION"
    ["SC19"]="PlasmidFinder Identification (StrainCascade_PlasmidFinder_identification.sh). Using PlasmidFinder v$PLASMIDFINDER_VERSION"
    ["SC20"]="AMRFinderPlus Antimicrobial Resistance Identification (StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh). Using AMRFinderPlus v$AMRFINDERPLUS_VERSION"
    ["SC21"]="ResFinder Antimicrobial Resistance Identification (StrainCascade_ResFinder_antimicrobial_resistance_identification.sh). Using ResFinder v$RESFINDER_VERSION"
    ["SC22"]="dbCAN3 CAZymes Identification (StrainCascade_dbCAN3_CAZymes_identification.sh). Using dbCAN3 v$DBCAN_VERSION"
    ["SC23"]="IslandPath Genomic Islands Identification (StrainCascade_IslandPath_genomic_islands_identification.sh). Using IslandPath-DIMOB v$ISLANDPATH_VERSION"
    ["SC24"]="VirSorter2 Phage Identification (StrainCascade_VirSorter2_phage_identification.sh). Using VirSorter2 v$VIRSORTER2_VERSION"
    ["SC25"]="DeepVirFinder Phage Identification (StrainCascade_DeepVirFinder_phage_identification.sh). Using DeepVirFinder v$DEEPVIRFINDER_VERSION"
    ["SC26"]="CRISPRCasFinder CRISPRCas Identification (StrainCascade_CRISPRCasFinder_identification.sh). Using CRISPRCasFinder v$CRISPRCASFINDER_VERSION"
    ["SC27"]="ISEScan Insertion Sequence Elements (StrainCascade_ISEScan_IS_elements_identification.sh). Using ISEScan v$ISESCAN_VERSION"
    ["SC28"]="Data Integration (StrainCascade_data_integration.sh). Using R v$R_VERSION"
)

# Check source exists and copy logo with error handling
if [ -f "${SCRIPT_DIR}/../assets/straincascade_logo_colour.png" ]; then
    mkdir -p "${RESULTS_INTEGRATION_DIR}"
    cp "${SCRIPT_DIR}/../assets/straincascade_logo_colour.png" "${RESULTS_INTEGRATION_DIR}/" || {
        echo "Error: Failed to copy logo file"
    }
else
    echo "Error: Logo file not found"
fi

readonly LOGO="${RESULTS_INTEGRATION_DIR}/straincascade_logo_colour.png"
readonly MD_FILE="${RESULTS_INTEGRATION_DIR}/StrainCascade_run_documentation.md"
{
    # Title and Header
    cat << EOF
---
title: "StrainCascade Analysis Report"
author: "StrainCascade v${VERSION} | Developed by Sebastian B.U. Jordi et al."
date: "$(date "+%Y-%m-%d %H:%M:%S")"
logo: "straincascade_logo_colour.png"
---

# StrainCascade Run Summary

## Executive Summary

This report documents the execution of the StrainCascade v${VERSION} for analyzing bacterial isolate whole genome sequencing data. The analysis was performed on sample "${SAMPLE_NAME}" using the following configuration: input type = ${INPUT_TYPE}, sequencing type = ${SEQUENCING_TYPE}, assembly selection algorithm = ${SELECTION_ALGORITHM}, reproducibility mode = ${REPRODUCIBILITY_MODE}. For further details, please refer to https://github.com/SBUJordi/StrainCascade.

EOF

    # System Information
    printf "## Local System Configuration and Resources\n\n"
    printf "| Component | Specification |\n"
    printf "|:----------|---------------:|\n"
    printf "| **Operating System** | \`%s\` |\n" "$OS_INFO"
    printf "| **OS Version** | \`%s\` |\n" "$OS_VERSION"
    printf "| **Kernel Version** | \`%s\` |\n" "$KERNEL_VERSION"
    printf "| **Architecture** | \`%s\` |\n" "$ARCHITECTURE"
    printf "| **CPU Model** | \`%s\` |\n" "$CPU_MODEL"
    printf "| **Available Memory** | \`%s\` |\n" "$AVAILABLE_MEMORY"
    printf "| **CPU Threads Used** | \`%s\` |\n" "$THREADS"
    printf "| **Total Runtime** | \`%s\` |\n\n" "$RUN_DURATION_STRING"

    # Function to print versioned tables
    print_version_table() {
        local title="$1"
        shift
        printf "### %s\n\n" "$title"
        printf "| Software | Version |\n"
        printf "|:---------|---------:|\n"
        while [ $# -gt 0 ]; do
            local software="$1"
            local version="$2"
            printf "| %s | \`%s\` |\n" "$software" "$version"
            shift 2
        done
        printf "\n"
    }

    # Software Versions
    printf "## Software Versions\n\n"
    
    print_version_table "Core Tools" \
        "R" "$R_VERSION" \
        "Python" "$PYTHON_VERSION"

    print_version_table "Assembly Tools" \
        "La Jolla Assembler (LJA)" "$LJA_VERSION" \
        "SPAdes" "$SPADES_VERSION" \
        "Canu" "$CANU_VERSION" \
        "Flye" "$FLYE_VERSION"

    print_version_table "QC & Refinement Tools" \
        "MAC.2" "$MAC2_VERSION" \
        "Circlator" "$CIRCLATOR_VERSION" \
        "Arrow (GenomicConsensus)" "$ARROW_VERSION" \
        "Medaka" "$MEDAKA_VERSION" \
        "NGMLR" "$NGMLR_VERSION" \
        "BBMAP" "$BBMAP_VERSION" \
        "QUAST" "$QUAST_VERSION" \
        "CheckM2" "$CHECKM2_VERSION"

    print_version_table "Taxonomy Tools" \
        "GTDB-Tk" "v$GTDBTK_VERSION"

    print_version_table "Genome Annotation Tools" \
        "Bakta" "$BAKTA_VERSION" \
        "Prokka" "$PROKKA_VERSION" \
        "MicrobeAnnotator" "$MICROBEANNOTATOR_VERSION"

    print_version_table "Functional Analysis Tools" \
        "GTDB-Tk" "$GTDBTK_VERSION" \
        "PlasmidFinder" "$PLASMIDFINDER_VERSION" \
        "AMRFinderPlus" "$AMRFINDERPLUS_VERSION" \
        "ResFinder" "$RESFINDER_VERSION" \
        "dbCAN3" "$DBCAN_VERSION" \
        "IslandPath" "$ISLANDPATH_VERSION" \
        "VirSorter2" "$VIRSORTER2_VERSION" \
        "DeepVirFinder" "$DEEPVIRFINDER_VERSION" \
        "CRISPRCasFinder" "$CRISPRCASFINDER_VERSION" \
        "ISEScan" "$ISESCAN_VERSION"

    # Executed Modules
    printf "## Executed SC Modules\n\n"
    printf "| Module ID | Description |\n"
    printf "|:----------|-------------:|\n"
    for module in $SELECTED_MODULES; do
        for key in "${!MODULE_INFO[@]}"; do
            if [[ "${MODULE_INFO[$key],,}" == *"${module,,}"* ]]; then
                printf "| \`%s\` | %s |\n" "$key" "${MODULE_INFO[$key]}"
            fi
        done
    done
    printf "\n"

    # Data Integrity Section
    printf "## Data Integrity\n\n"

    # Combined Files Hash Table
    printf "### Data Files\n\n"
    printf "| Description                  |\n"
    printf "|:-----------------------------|\n"

    # Input file
    if [[ -f "$INPUT_FILE" ]]; then
        filename="$(basename "$INPUT_FILE")"
        hash="$(sha256sum "$INPUT_FILE" | cut -d' ' -f1)"
    else
        filename="N/A"
        filepath="N/A"
        hash="N/A (file not found)"
    fi
    printf "| **Input File**               |\n"
    printf "| File Name: \`%s\`              |\n" "$filename"
    printf "| Hash: \`%s\`                   |\n" "$hash"
    printf "| ---                           |\n"

    # Assembly file
    if [[ -f "$analysis_assembly_file" ]]; then
        filename="$(basename "$analysis_assembly_file")"
        filepath="$analysis_assembly_file"
        hash="$(sha256sum "$analysis_assembly_file" | cut -d' ' -f1)"
    else
        filename="N/A"
        filepath="N/A"
        hash="N/A (file not found)"
    fi
    printf "| **Analysis Assembly File**   |\n"
    printf "| File Name: \`%s\`              |\n" "$filename"
    printf "| Hash: \`%s\`                   |\n" "$hash"
    printf "| ---                           |\n"

    # RData file
    if [[ -f "$Rdata_file" ]]; then
        filename="$(basename "$Rdata_file")"
        filepath="$Rdata_file"
        hash="$(sha256sum "$Rdata_file" | cut -d' ' -f1)"
    else
        filename="N/A"
        filepath="N/A"
        hash="N/A (file not found)"
    fi
    printf "| **R Result File**            |\n"
    printf "| File Name: \`%s\`              |\n" "$filename"
    printf "| Hash: \`%s\`                   |\n" "$hash"
    printf "| ---                           |\n"

    # HTML file
    if [[ -f "$HTML_file" ]]; then
        filename="$(basename "$HTML_file")"
        filepath="$HTML_file"
        hash="$(sha256sum "$HTML_file" | cut -d' ' -f1)"
    else
        filename="N/A"
        filepath="N/A"
        hash="N/A (file not found)"
    fi
    printf "| **Result Summary File**      |\n"
    printf "| File Name: \`%s\`              |\n" "$filename"
    printf "| Hash: \`%s\`                   |\n" "$hash"
    printf "| ---                           |\n"

    # Apptainer Images Hash Table
    printf "### Apptainer Images\n\n"
    printf "| File Hashes (SHA256) |\n"
    printf "|:------------|\n"
    for file in "$python_3_12_4_sif" \
            "$r_4_4_1_sif" \
            "$straincascade_genome_assembly_sif" \
            "$straincascade_lja_genome_assembly_sif" \
            "$straincascade_assembly_qc_refinement_sif" \
            "$straincascade_genome_annotation_sif" \
            "$straincascade_taxonomic_functional_analysis_sif" \
            "$straincascade_crisprcas_phage_is_elements_sif"; do
        filename="$(basename "$file")"
        hash="$(sha256sum "$file" | cut -d' ' -f1)"
        printf "| **\`%s\`** |\n" "$filename"
        printf "| \`%s\` |\n" "$hash"
        printf "| --- |\n"
    done
    printf "\n"

    # Scripts Directory Hash
    printf "### StrainCascade Software\n\n"
    printf "| File Hashes (SHA256) |\n"
    printf "|:------------|\n"
    printf "| **\`%s\`** |\n" "$SCRIPT_DIR"
    printf "| \`%s\` |\n" "$(tar --mode=644 --mtime='2023-01-01 00:00:00 UTC' -cf - "$SCRIPT_DIR" | sha256sum | awk '{print $1}')"
    printf "| --- |\n\n"

    # Hash Verification Instructions
    printf "### Hash Verification Instructions\n\n"
    printf "To verify hashes, use these commands:\n\n"
    printf -- "- For files and Apptainer images:\n"
    printf "  \`\`\`bash\n"
    printf "  sha256sum path_to_file | awk '{print \$1}'\n"
    printf "  \`\`\`\n\n"
    printf -- "- For scripts directory:\n"
    printf "  \`\`\`bash\n"
    printf -- "  tar --mode=644 --mtime='2023-01-01 00:00:00 UTC'\n"
    printf -- "-cf - path_to_scripts_dir | sha256sum | awk '{print \$1}'\n"
    printf "  \`\`\`\n\n"

    # Add page break
    #printf "<div style=\"page-break-after: always;\"></div>\n\n"
    
    # Add log section
    #printf "### StrainCascade run log\n\n"
    #printf "\`\`\`\n"
    #cat "$LOGS_DIR/$LOG_NAME" 2>/dev/null || printf "Log file not found or not readable\n"
    #printf "\`\`\`\n\n"
} > "$MD_FILE"

# Write LaTeX template
# Create template file
readonly TEMPLATE_FILE="${RESULTS_INTEGRATION_DIR}/straincascade_template.tex"

cat << 'EOF' > "$TEMPLATE_FILE"
%%
% Copyright (c) 2017 - 2024, Pascal Wagler;
% Copyright (c) 2014 - 2024, John MacFarlane
%
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%
% - Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
%
% - Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% - Neither the name of John MacFarlane nor the names of other
% contributors may be used to endorse or promote products derived
% from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%%

%%
% This is the Eisvogel pandoc LaTeX template.
%
% For usage information and examples visit the official GitHub page:
% https://github.com/Wandmalfarbe/pandoc-latex-template
%%

% Options for packages loaded elsewhere
\PassOptionsToPackage{unicode$for(hyperrefoptions)$,$hyperrefoptions$$endfor$}{hyperref}
\PassOptionsToPackage{hyphens}{url}
\PassOptionsToPackage{dvipsnames,svgnames,x11names,table}{xcolor}
$if(CJKmainfont)$
\PassOptionsToPackage{space}{xeCJK}
$endif$
%
\documentclass[
$if(fontsize)$
  $fontsize$,
$endif$
$if(papersize)$
  $papersize$paper,
$else$
  paper=a4,
$endif$
$if(beamer)$
  ignorenonframetext,
$if(handout)$
  handout,
$endif$
$if(aspectratio)$
  aspectratio=$aspectratio$,
$endif$
$if(babel-lang)$
  $babel-lang$,
$endif$
$endif$
$for(classoption)$
  $classoption$$sep$,
$endfor$
  ,captions=tableheading
]{$if(beamer)$$documentclass$$else$$if(book)$scrbook$else$scrartcl$endif$$endif$}
$if(beamer)$
$if(background-image)$
\usebackgroundtemplate{%
  \includegraphics[width=\paperwidth]{$background-image$}%
}
% In beamer background-image does not work well when other images are used, so this is the workaround
\pgfdeclareimage[width=\paperwidth,height=\paperheight]{background}{$background-image$}
\usebackgroundtemplate{\pgfuseimage{background}}
$endif$
\usepackage{pgfpages}
\setbeamertemplate{caption}[numbered]
\setbeamertemplate{caption label separator}{: }
\setbeamercolor{caption name}{fg=normal text.fg}
\beamertemplatenavigationsymbols$if(navigation)$$navigation$$else$empty$endif$
$for(beameroption)$
\setbeameroption{$beameroption$}
$endfor$
% Prevent slide breaks in the middle of a paragraph
\widowpenalties 1 10000
\raggedbottom
$if(section-titles)$
\setbeamertemplate{part page}{
  \centering
  \begin{beamercolorbox}[sep=16pt,center]{part title}
    \usebeamerfont{part title}\insertpart\par
  \end{beamercolorbox}
}
\setbeamertemplate{section page}{
  \centering
  \begin{beamercolorbox}[sep=12pt,center]{section title}
    \usebeamerfont{section title}\insertsection\par
  \end{beamercolorbox}
}
\setbeamertemplate{subsection page}{
  \centering
  \begin{beamercolorbox}[sep=8pt,center]{subsection title}
    \usebeamerfont{subsection title}\insertsubsection\par
  \end{beamercolorbox}
}
\AtBeginPart{
  \frame{\partpage}
}
\AtBeginSection{
  \ifbibliography
  \else
    \frame{\sectionpage}
  \fi
}
\AtBeginSubsection{
  \frame{\subsectionpage}
}
$endif$
$endif$
$if(beamerarticle)$
\usepackage{beamerarticle} % needs to be loaded first
$endif$
\usepackage{amsmath,amssymb}
$if(linestretch)$
\usepackage{setspace}
$else$
% Use setspace anyway because we change the default line spacing.
% The spacing is changed early to affect the titlepage and the TOC.
\usepackage{setspace}
\setstretch{1.2}
$endif$
\usepackage{iftex}
\ifPDFTeX
  \usepackage[$if(fontenc)$$fontenc$$else$T1$endif$]{fontenc}
  \usepackage[utf8]{inputenc}
  \usepackage{textcomp} % provide euro and other symbols
\else % if luatex or xetex
$if(mathspec)$
  \ifXeTeX
    \usepackage{mathspec} % this also loads fontspec
  \else
    \usepackage{unicode-math} % this also loads fontspec
  \fi
$else$
  \usepackage{unicode-math} % this also loads fontspec
$endif$
  \defaultfontfeatures{Scale=MatchLowercase}$-- must come before Beamer theme
  \defaultfontfeatures[\rmfamily]{Ligatures=TeX,Scale=1}
\fi
$if(fontfamily)$
$else$
$-- Set default font before Beamer theme so the theme can override it
\usepackage{lmodern}
$endif$
$-- Set Beamer theme before user font settings so they can override theme
$if(beamer)$
$if(theme)$
\usetheme[$for(themeoptions)$$themeoptions$$sep$,$endfor$]{$theme$}
$endif$
$if(colortheme)$
\usecolortheme{$colortheme$}
$endif$
$if(fonttheme)$
\usefonttheme{$fonttheme$}
$endif$
$if(mainfont)$
\usefonttheme{serif} % use mainfont rather than sansfont for slide text
$endif$
$if(innertheme)$
\useinnertheme{$innertheme$}
$endif$
$if(outertheme)$
\useoutertheme{$outertheme$}
$endif$
$endif$
$-- User font settings (must come after default font and Beamer theme)
$if(fontfamily)$
\usepackage[$for(fontfamilyoptions)$$fontfamilyoptions$$sep$,$endfor$]{$fontfamily$}
$endif$
\ifPDFTeX\else
  % xetex/luatex font selection
$if(mainfont)$
  $if(mainfontfallback)$
    \ifLuaTeX
      \usepackage{luaotfload}
      \directlua{luaotfload.add_fallback("mainfontfallback",{
        $for(mainfontfallback)$"$mainfontfallback$"$sep$,$endfor$
      })}
    \fi
  $endif$
  \setmainfont[$for(mainfontoptions)$$mainfontoptions$$sep$,$endfor$$if(mainfontfallback)$,RawFeature={fallback=mainfontfallback}$endif$]{$mainfont$}
$endif$
$if(sansfont)$
  $if(sansfontfallback)$
    \ifLuaTeX
      \usepackage{luaotfload}
      \directlua{luaotfload.add_fallback("sansfontfallback",{
        $for(sansfontfallback)$"$sansfontfallback$"$sep$,$endfor$
      })}
    \fi
  $endif$
  \setsansfont[$for(sansfontoptions)$$sansfontoptions$$sep$,$endfor$$if(sansfontfallback)$,RawFeature={fallback=sansfontfallback}$endif$]{$sansfont$}
$endif$
$if(monofont)$
  $if(monofontfallback)$
    \ifLuaTeX
      \usepackage{luaotfload}
      \directlua{luaotfload.add_fallback("monofontfallback",{
        $for(monofontfallback)$"$monofontfallback$"$sep$,$endfor$
      })}
    \fi
  $endif$
  \setmonofont[$for(monofontoptions)$$monofontoptions$$sep$,$endfor$$if(monofontfallback)$,RawFeature={fallback=monofontfallback}$endif$]{$monofont$}
$endif$
$for(fontfamilies)$
  \newfontfamily{$fontfamilies.name$}[$for(fontfamilies.options)$$fontfamilies.options$$sep$,$endfor$]{$fontfamilies.font$}
$endfor$
$if(mathfont)$
$if(mathspec)$
  \ifXeTeX
    \setmathfont(Digits,Latin,Greek)[$for(mathfontoptions)$$mathfontoptions$$sep$,$endfor$]{$mathfont$}
  \else
    \setmathfont[$for(mathfontoptions)$$mathfontoptions$$sep$,$endfor$]{$mathfont$}
  \fi
$else$
  \setmathfont[$for(mathfontoptions)$$mathfontoptions$$sep$,$endfor$]{$mathfont$}
$endif$
$endif$
$if(CJKmainfont)$
  \ifXeTeX
    \usepackage{xeCJK}
    \setCJKmainfont[$for(CJKoptions)$$CJKoptions$$sep$,$endfor$]{$CJKmainfont$}
    $if(CJKsansfont)$
      \setCJKsansfont[$for(CJKoptions)$$CJKoptions$$sep$,$endfor$]{$CJKsansfont$}
    $endif$
    $if(CJKmonofont)$
      \setCJKmonofont[$for(CJKoptions)$$CJKoptions$$sep$,$endfor$]{$CJKmonofont$}
    $endif$
  \fi
$endif$
$if(luatexjapresetoptions)$
  \ifLuaTeX
    \usepackage[$for(luatexjapresetoptions)$$luatexjapresetoptions$$sep$,$endfor$]{luatexja-preset}
  \fi
$endif$
$if(CJKmainfont)$
  \ifLuaTeX
    \usepackage[$for(luatexjafontspecoptions)$$luatexjafontspecoptions$$sep$,$endfor$]{luatexja-fontspec}
    \setmainjfont[$for(CJKoptions)$$CJKoptions$$sep$,$endfor$]{$CJKmainfont$}
  \fi
$endif$
\fi
$if(zero-width-non-joiner)$
%% Support for zero-width non-joiner characters.
\makeatletter
\def\zerowidthnonjoiner{%
  % Prevent ligatures and adjust kerning, but still support hyphenating.
  \texorpdfstring{%
    \TextOrMath{\nobreak\discretionary{-}{}{\kern.03em}%
      \ifvmode\else\nobreak\hskip\z@skip\fi}{}%
  }{}%
}
\makeatother
\ifPDFTeX
  \DeclareUnicodeCharacter{200C}{\zerowidthnonjoiner}
\else
  \catcode`^^^^200c=\active
  \protected\def ^^^^200c{\zerowidthnonjoiner}
\fi
%% End of ZWNJ support
$endif$
% Use upquote if available, for straight quotes in verbatim environments
\IfFileExists{upquote.sty}{\usepackage{upquote}}{}
\IfFileExists{microtype.sty}{% use microtype if available
  \usepackage[$for(microtypeoptions)$$microtypeoptions$$sep$,$endfor$]{microtype}
  \UseMicrotypeSet[protrusion]{basicmath} % disable protrusion for tt fonts
}{}
$if(indent)$
$else$
\makeatletter
\@ifundefined{KOMAClassName}{% if non-KOMA class
  \IfFileExists{parskip.sty}{%
    \usepackage{parskip}
  }{% else
    \setlength{\parindent}{0pt}
    \setlength{\parskip}{6pt plus 2pt minus 1pt}}
}{% if KOMA class
  \KOMAoptions{parskip=half}}
\makeatother
$endif$
$if(verbatim-in-note)$
\usepackage{fancyvrb}
$endif$
\usepackage{xcolor}
\definecolor{default-linkcolor}{HTML}{A50000}
\definecolor{default-filecolor}{HTML}{A50000}
\definecolor{default-citecolor}{HTML}{4077C0}
\definecolor{default-urlcolor}{HTML}{4077C0}
$if(footnotes-pretty)$
% load footmisc in order to customize footnotes (footmisc has to be loaded before hyperref, cf. https://tex.stackexchange.com/a/169124/144087)
\usepackage[hang,flushmargin,bottom,multiple]{footmisc}
\setlength{\footnotemargin}{0.8em} % set space between footnote nr and text
\setlength{\footnotesep}{\baselineskip} % set space between multiple footnotes
\setlength{\skip\footins}{0.3cm} % set space between page content and footnote
\setlength{\footskip}{0.9cm} % set space between footnote and page bottom
$endif$
$if(geometry)$
$if(beamer)$
\geometry{$for(geometry)$$geometry$$sep$,$endfor$}
$else$
\usepackage[$for(geometry)$$geometry$$sep$,$endfor$]{geometry}
$endif$
$else$
$if(beamer)$
$else$
\usepackage[margin=2.5cm,includehead=true,includefoot=true,centering,$for(geometry)$$geometry$$sep$,$endfor$]{geometry}
$endif$
$endif$
$if(titlepage-logo)$
\usepackage[export]{adjustbox}
\usepackage{graphicx}
$endif$
$if(beamer)$
\newif\ifbibliography
$endif$
$if(listings)$
\usepackage{listings}
\newcommand{\passthrough}[1]{#1}
\lstset{defaultdialect=[5.3]Lua}
\lstset{defaultdialect=[x86masm]Assembler}
$endif$
$if(listings-no-page-break)$
\usepackage{etoolbox}
\BeforeBeginEnvironment{lstlisting}{\par\noindent\begin{minipage}{\linewidth}}
\AfterEndEnvironment{lstlisting}{\end{minipage}\par\addvspace{\topskip}}
$endif$
$if(lhs)$
\lstnewenvironment{code}{\lstset{language=Haskell,basicstyle=\small\ttfamily}}{}
$endif$
$if(highlighting-macros)$
$highlighting-macros$

% Workaround/bugfix from jannick0.
% See https://github.com/jgm/pandoc/issues/4302#issuecomment-360669013)
% or https://github.com/Wandmalfarbe/pandoc-latex-template/issues/2
%
% Redefine the verbatim environment 'Highlighting' to break long lines (with
% the help of fvextra). Redefinition is necessary because it is unlikely that
% pandoc includes fvextra in the default template.
\usepackage{fvextra}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,fontsize=$if(code-block-font-size)$$code-block-font-size$$else$\small$endif$,commandchars=\\\{\}}

$endif$
$if(tables)$
\usepackage{longtable,booktabs,array}
$if(multirow)$
\usepackage{multirow}
$endif$
\usepackage{calc} % for calculating minipage widths
$if(beamer)$
\usepackage{caption}
% Make caption package work with longtable
\makeatletter
\def\fnum@table{\tablename~\thetable}
\makeatother
$else$
% Correct order of tables after \paragraph or \subparagraph
\usepackage{etoolbox}
\makeatletter
\patchcmd\longtable{\par}{\if@noskipsec\mbox{}\fi\par}{}{}
\makeatother
% Allow footnotes in longtable head/foot
\IfFileExists{footnotehyper.sty}{\usepackage{footnotehyper}}{\usepackage{footnote}}
\makesavenoteenv{longtable}
$endif$
$endif$
% add backlinks to footnote references, cf. https://tex.stackexchange.com/questions/302266/make-footnote-clickable-both-ways
$if(footnotes-disable-backlinks)$
$else$
\usepackage{footnotebackref}
$endif$
$if(graphics)$
\usepackage{graphicx}
\makeatletter
\newsavebox\pandoc@box
\newcommand*\pandocbounded[1]{% scales image to fit in text height/width
  \sbox\pandoc@box{#1}%
  \Gscale@div\@tempa{\textheight}{\dimexpr\ht\pandoc@box+\dp\pandoc@box\relax}%
  \Gscale@div\@tempb{\linewidth}{\wd\pandoc@box}%
  \ifdim\@tempb\p@<\@tempa\p@\let\@tempa\@tempb\fi% select the smaller of both
  \ifdim\@tempa\p@<\p@\scalebox{\@tempa}{\usebox\pandoc@box}%
  \else\usebox{\pandoc@box}%
  \fi%
}
% Set default figure placement to htbp
% Make use of float-package and set default placement for figures to H.
% The option H means 'PUT IT HERE' (as  opposed to the standard h option which means 'You may put it here if you like').
\usepackage{float}
\floatplacement{figure}{$if(float-placement-figure)$$float-placement-figure$$else$H$endif$}
\makeatother
$endif$
$if(svg)$
\usepackage{svg}
$endif$
$if(strikeout)$
$-- also used for underline
\ifLuaTeX
  \usepackage{luacolor}
  \usepackage[soul]{lua-ul}
\else
\usepackage{soul}
$if(beamer)$
  \makeatletter
  \let\HL\hl
  \renewcommand\hl{% fix for beamer highlighting
    \let\set@color\beamerorig@set@color
    \let\reset@color\beamerorig@reset@color
    \HL}
  \makeatother
$endif$
$if(CJKmainfont)$
  \ifXeTeX
    % soul's \st doesn't work for CJK:
    \usepackage{xeCJKfntef}
    \renewcommand{\st}[1]{\sout{#1}}
  \fi
$endif$
\fi
$endif$
\setlength{\emergencystretch}{3em} % prevent overfull lines
\providecommand{\tightlist}{%
  \setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
$if(numbersections)$
\setcounter{secnumdepth}{$if(secnumdepth)$$secnumdepth$$else$5$endif$}
$else$
\setcounter{secnumdepth}{-\maxdimen} % remove section numbering
$endif$
$if(subfigure)$
\usepackage{subcaption}
$endif$
$if(beamer)$
$else$
$if(block-headings)$
% Make \paragraph and \subparagraph free-standing
\makeatletter
\ifx\paragraph\undefined\else
  \let\oldparagraph\paragraph
  \renewcommand{\paragraph}{
    \@ifstar
      \xxxParagraphStar
      \xxxParagraphNoStar
  }
  \newcommand{\xxxParagraphStar}[1]{\oldparagraph*{#1}\mbox{}}
  \newcommand{\xxxParagraphNoStar}[1]{\oldparagraph{#1}\mbox{}}
\fi
\ifx\subparagraph\undefined\else
  \let\oldsubparagraph\subparagraph
  \renewcommand{\subparagraph}{
    \@ifstar
      \xxxSubParagraphStar
      \xxxSubParagraphNoStar
  }
  \newcommand{\xxxSubParagraphStar}[1]{\oldsubparagraph*{#1}\mbox{}}
  \newcommand{\xxxSubParagraphNoStar}[1]{\oldsubparagraph{#1}\mbox{}}
\fi
\makeatother
$endif$
$endif$
$if(pagestyle)$
\pagestyle{$pagestyle$}
$endif$
$if(csl-refs)$
% definitions for citeproc citations
\NewDocumentCommand\citeproctext{}{}
\NewDocumentCommand\citeproc{mm}{%
  \begingroup\def\citeproctext{#2}\cite{#1}\endgroup}
\makeatletter
 % allow citations to break across lines
 \let\@cite@ofmt\@firstofone
 % avoid brackets around text for \cite:
 \def\@biblabel#1{}
 \def\@cite#1#2{{#1\if@tempswa , #2\fi}}
\makeatother
\newlength{\cslhangindent}
\setlength{\cslhangindent}{1.5em}
\newlength{\csllabelwidth}
\setlength{\csllabelwidth}{3em}
\newenvironment{CSLReferences}[2] % #1 hanging-indent, #2 entry-spacing
  {\begin{list}{}{%
   \setlength{\itemindent}{0pt}
   \setlength{\leftmargin}{0pt}
   \setlength{\parsep}{0pt}
   % turn on hanging indent if param 1 is 1
   \ifodd #1
    \setlength{\leftmargin}{\cslhangindent}
    \setlength{\itemindent}{-1\cslhangindent}
   \fi
   % set entry spacing
   \setlength{\itemsep}{#2\baselineskip}}}
  {\end{list}}
\usepackage{calc}
\newcommand{\CSLBlock}[1]{\hfill\break\parbox[t]{\linewidth}{\strut\ignorespaces#1\strut}}
\newcommand{\CSLLeftMargin}[1]{\parbox[t]{\csllabelwidth}{\strut#1\strut}}
\newcommand{\CSLRightInline}[1]{\parbox[t]{\linewidth - \csllabelwidth}{\strut#1\strut}}
\newcommand{\CSLIndent}[1]{\hspace{\cslhangindent}#1}
$endif$
$if(lang)$
\ifLuaTeX
\usepackage[bidi=basic]{babel}
\else
\usepackage[bidi=default]{babel}
\fi
$if(babel-lang)$
\babelprovide[main,import]{$babel-lang$}
$if(mainfont)$
\ifPDFTeX
\else
\babelfont{rm}[$for(mainfontoptions)$$mainfontoptions$$sep$,$endfor$$if(mainfontfallback)$,RawFeature={fallback=mainfontfallback}$endif$]{$mainfont$}
\fi
$endif$
$endif$
$for(babel-otherlangs)$
\babelprovide[import]{$babel-otherlangs$}
$endfor$
$for(babelfonts/pairs)$
\babelfont[$babelfonts.key$]{rm}{$babelfonts.value$}
$endfor$
% get rid of language-specific shorthands (see #6817):
\let\LanguageShortHands\languageshorthands
\def\languageshorthands#1{}
$if(selnolig-langs)$
\ifLuaTeX
  \usepackage[$for(selnolig-langs)$$it$$sep$,$endfor$]{selnolig} % disable illegal ligatures
\fi
$endif$
$endif$
$for(header-includes)$
$header-includes$
$endfor$
$if(dir)$
\ifPDFTeX
  \TeXXeTstate=1
  \newcommand{\RL}[1]{\beginR #1\endR}
  \newcommand{\LR}[1]{\beginL #1\endL}
  \newenvironment{RTL}{\beginR}{\endR}
  \newenvironment{LTR}{\beginL}{\endL}
\fi
$endif$
$if(natbib)$
\usepackage[$natbiboptions$]{natbib}
\bibliographystyle{$if(biblio-style)$$biblio-style$$else$plainnat$endif$}
$endif$
$if(biblatex)$
\usepackage[$if(biblio-style)$style=$biblio-style$,$endif$$for(biblatexoptions)$$biblatexoptions$$sep$,$endfor$]{biblatex}
$for(bibliography)$
\addbibresource{$bibliography$}
$endfor$
$endif$
$if(nocite-ids)$
\nocite{$for(nocite-ids)$$it$$sep$, $endfor$}
$endif$
$if(csquotes)$
\usepackage{csquotes}
$endif$
\usepackage{bookmark}
\IfFileExists{xurl.sty}{\usepackage{xurl}}{} % add URL line breaks if available
\urlstyle{$if(urlstyle)$$urlstyle$$else$same$endif$}
$if(links-as-notes)$
% Make links footnotes instead of hotlinks:
\DeclareRobustCommand{\href}[2]{#2\footnote{\url{#1}}}
$endif$
$if(verbatim-in-note)$
\VerbatimFootnotes % allow verbatim text in footnotes
$endif$
\hypersetup{
$if(title-meta)$
  pdftitle={$title-meta$},
$endif$
$if(author-meta)$
  pdfauthor={$author-meta$},
$endif$
$if(lang)$
  pdflang={$lang$},
$endif$
$if(subject)$
  pdfsubject={$subject$},
$endif$
$if(keywords)$
  pdfkeywords={$for(keywords)$$keywords$$sep$, $endfor$},
$endif$
$if(colorlinks)$
  colorlinks=true,
  linkcolor={$if(linkcolor)$$linkcolor$$else$default-linkcolor$endif$},
  filecolor={$if(filecolor)$$filecolor$$else$default-filecolor$endif$},
  citecolor={$if(citecolor)$$citecolor$$else$default-citecolor$endif$},
  urlcolor={$if(urlcolor)$$urlcolor$$else$default-urlcolor$endif$},
$else$
$if(boxlinks)$
$else$
  hidelinks,
$endif$
$endif$
  breaklinks=true,
  pdfcreator={LaTeX via pandoc with the Eisvogel template}}
$if(title)$
\title{$title$$if(thanks)$\thanks{$thanks$}$endif$}
$endif$
$if(subtitle)$
$if(beamer)$
$else$
\usepackage{etoolbox}
\makeatletter
\providecommand{\subtitle}[1]{% add subtitle to \maketitle
  \apptocmd{\@title}{\par {\large #1 \par}}{}{}
}
\makeatother
$endif$
\subtitle{$subtitle$}
$endif$
\author{$for(author)$$author$$sep$ \and $endfor$}
\date{$date$}
$if(beamer)$
$if(institute)$
\institute{$for(institute)$$institute$$sep$ \and $endfor$}
$endif$
$if(titlegraphic)$
\titlegraphic{\includegraphics$if(titlegraphicoptions)$[$for(titlegraphicoptions)$$titlegraphicoptions$$sep$, $endfor$]$endif${$titlegraphic$}}
$endif$
$if(logo)$
\logo{\includegraphics{$logo$}}
$endif$
$endif$



%%
%% added
%%

$if(page-background)$
\usepackage[pages=all]{background}
$endif$

%
% for the background color of the title page
%
$if(titlepage)$
\usepackage{pagecolor}
\usepackage{afterpage}
$if(titlepage-background)$
\usepackage{tikz}
$endif$
$if(geometry)$
$else$
\usepackage[margin=2.5cm,includehead=true,includefoot=true,centering]{geometry}
$endif$
$endif$

%
% break urls
%
\PassOptionsToPackage{hyphens}{url}

%
% When using babel or polyglossia with biblatex, loading csquotes is recommended
% to ensure that quoted texts are typeset according to the rules of your main language.
%
\usepackage{csquotes}

%
% captions
%
\definecolor{caption-color}{HTML}{777777}
$if(beamer)$
$else$
\usepackage[font={stretch=1.2}, textfont={color=caption-color}, position=top, skip=4mm, labelfont=bf, singlelinecheck=false, justification=$if(caption-justification)$$caption-justification$$else$raggedright$endif$]{caption}
\setcapindent{0em}
$endif$

%
% blockquote
%
\definecolor{blockquote-border}{RGB}{221,221,221}
\definecolor{blockquote-text}{RGB}{119,119,119}
\usepackage{mdframed}
\newmdenv[rightline=false,bottomline=false,topline=false,linewidth=3pt,linecolor=blockquote-border,skipabove=\parskip]{customblockquote}
\renewenvironment{quote}{\begin{customblockquote}\list{}{\rightmargin=0em\leftmargin=0em}%
\item\relax\color{blockquote-text}\ignorespaces}{\unskip\unskip\endlist\end{customblockquote}}

%
% Source Sans Pro as the default font family
% Source Code Pro for monospace text
%
% 'default' option sets the default
% font family to Source Sans Pro, not \sfdefault.
%
\ifnum 0\ifxetex 1\fi\ifluatex 1\fi=0 % if pdftex
  $if(fontfamily)$
  $else$
  \usepackage[default]{sourcesanspro}
  \usepackage{sourcecodepro}
  $endif$
\else % if not pdftex
  $if(mainfont)$
  $else$
  \usepackage[default]{sourcesanspro}
  \usepackage{sourcecodepro}

  % XeLaTeX specific adjustments for straight quotes: https://tex.stackexchange.com/a/354887
  % This issue is already fixed (see https://github.com/silkeh/latex-sourcecodepro/pull/5) but the
  % fix is still unreleased.
  % TODO: Remove this workaround when the new version of sourcecodepro is released on CTAN.
  \ifxetex
    \makeatletter
    \defaultfontfeatures[\ttfamily]
      { Numbers   = \sourcecodepro@figurestyle,
        Scale     = \SourceCodePro@scale,
        Extension = .otf }
    \setmonofont
      [ UprightFont    = *-\sourcecodepro@regstyle,
        ItalicFont     = *-\sourcecodepro@regstyle It,
        BoldFont       = *-\sourcecodepro@boldstyle,
        BoldItalicFont = *-\sourcecodepro@boldstyle It ]
      {SourceCodePro}
    \makeatother
  \fi
  $endif$
\fi

%
% heading color
%
\definecolor{heading-color}{RGB}{40,40,40}
$if(beamer)$
$else$
\addtokomafont{section}{\color{heading-color}}
$endif$
% When using the classes report, scrreprt, book,
% scrbook or memoir, uncomment the following line.
%\addtokomafont{chapter}{\color{heading-color}}

%
% variables for title, author and date
%
$if(beamer)$
$else$
\usepackage{titling}
\title{$title$}
\author{$for(author)$$author$$sep$, $endfor$}
\date{$date$}
$endif$

%
% Table configuration
%
$if(tables)$

% Colors
\definecolor{table-row-color}{HTML}{F5F5F5}
\definecolor{table-rule-color}{HTML}{999999}

% Basic table settings
\arrayrulecolor{table-rule-color}
\setlength\heavyrulewidth{0.3ex}
\renewcommand{\arraystretch}{1.3}

% Required packages
\usepackage{tabularx}
\usepackage{ltablex}
\usepackage{booktabs}
\usepackage{array}
\usepackage{seqsplit}  % For splitting long strings
\usepackage{xstring}   % For string manipulation

% Table width and positioning
\newlength{\tablewidth}
\setlength{\tablewidth}{0.8\textwidth}  % Set to 80% of text width

% Enhanced column types with improved wrapping
\newcolumntype{L}{>{\raggedright\arraybackslash\hspace{0pt}\hsize=\hsize plus 1fill\linewidth=\hsize}p{\dimexpr\tablewidth/\real{3}}}
\newcolumntype{C}{>{\centering\arraybackslash\hspace{0pt}\hsize=\hsize plus 1fill\linewidth=\hsize}p{\dimexpr\tablewidth/\real{3}}}
\newcolumntype{R}{>{\raggedleft\arraybackslash\hspace{0pt}\hsize=\hsize plus 1fill\linewidth=\hsize}p{\dimexpr\tablewidth/\real{3}}}

% Configure longtable alignment - left-aligned
\AtBeginEnvironment{longtable}{%
    \setlength\LTleft{0pt}           % Set left margin to 0
    \setlength\LTright{\fill}        % Fill space on right side
}

% Optimize tabular environment with word wrap handling
\let\oldtabular\tabular
\let\endoldtabular\endtabular
\renewenvironment{tabular}[1]{%
    \tabularx{\tablewidth}{#1}%
}{%
    \endtabularx
}

% Define default column type for tables with improved spacing
\AtBeginEnvironment{longtable}{%
    \renewcommand{\arraystretch}{1.3}%
    \renewcommand{\tabcolsep}{0.5em}%
    \setlength{\emergencystretch}{3em}  % Allows emergency stretching
}

% Command for handling very long strings
\newcommand{\wrapcell}[1]{%
    \begin{minipage}[t]{\linewidth}%
        \raggedright\seqsplit{#1}%
    \end{minipage}%
}

% Optional row colors
$if(table-use-row-colors)$
\usepackage{etoolbox}
\AtBeginEnvironment{longtable}{%
    \rowcolors{2}{}{table-row-color!100}%
}
\preto{\toprule}{\hiderowcolors}
\appto{\endhead}{\showrowcolors}
\appto{\endfirsthead}{\showrowcolors}
$endif$

$endif$

% Spacing settings
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt plus 2pt minus 1pt}
\setlength{\emergencystretch}{3em}

%best pandoc templates
% remove paragraph indentation
%
\setlength{\parindent}{0pt}
\setlength{\parskip}{6pt plus 2pt minus 1pt}
\setlength{\emergencystretch}{3em}  % prevent overfull lines

%
%
% Listings
%
%

$if(listings)$

%
% general listing colors
%
\definecolor{listing-background}{HTML}{F7F7F7}
\definecolor{listing-rule}{HTML}{B3B2B3}
\definecolor{listing-numbers}{HTML}{B3B2B3}
\definecolor{listing-text-color}{HTML}{000000}
\definecolor{listing-keyword}{HTML}{435489}
\definecolor{listing-keyword-2}{HTML}{1284CA} % additional keywords
\definecolor{listing-keyword-3}{HTML}{9137CB} % additional keywords
\definecolor{listing-identifier}{HTML}{435489}
\definecolor{listing-string}{HTML}{00999A}
\definecolor{listing-comment}{HTML}{8E8E8E}

\lstdefinestyle{eisvogel_listing_style}{
  language         = java,
$if(listings-disable-line-numbers)$
  xleftmargin      = 0.6em,
  framexleftmargin = 0.4em,
$else$
  numbers          = left,
  xleftmargin      = 2.7em,
  framexleftmargin = 2.5em,
$endif$
  backgroundcolor  = \color{listing-background},
  basicstyle       = \color{listing-text-color}\linespread{1.0}%
                      \lst@ifdisplaystyle%
                      $if(code-block-font-size)$$code-block-font-size$$else$\small$endif$%
                      \fi\ttfamily{},
  breaklines       = true,
  frame            = single,
  framesep         = 0.19em,
  rulecolor        = \color{listing-rule},
  frameround       = ffff,
  tabsize          = 4,
  numberstyle      = \color{listing-numbers},
  aboveskip        = 1.0em,
  belowskip        = 0.1em,
  abovecaptionskip = 0em,
  belowcaptionskip = 1.0em,
  keywordstyle     = {\color{listing-keyword}\bfseries},
  keywordstyle     = {[2]\color{listing-keyword-2}\bfseries},
  keywordstyle     = {[3]\color{listing-keyword-3}\bfseries\itshape},
  sensitive        = true,
  identifierstyle  = \color{listing-identifier},
  commentstyle     = \color{listing-comment},
  stringstyle      = \color{listing-string},
  showstringspaces = false,
  escapeinside     = {/*@}{@*/}, % Allow LaTeX inside these special comments
  literate         =
  {á}{{\'a}}1 {é}{{\'e}}1 {í}{{\'i}}1 {ó}{{\'o}}1 {ú}{{\'u}}1
  {Á}{{\'A}}1 {É}{{\'E}}1 {Í}{{\'I}}1 {Ó}{{\'O}}1 {Ú}{{\'U}}1
  {à}{{\`a}}1 {è}{{\`e}}1 {ì}{{\`i}}1 {ò}{{\`o}}1 {ù}{{\`u}}1
  {À}{{\`A}}1 {È}{{\`E}}1 {Ì}{{\`I}}1 {Ò}{{\`O}}1 {Ù}{{\`U}}1
  {ä}{{\"a}}1 {ë}{{\"e}}1 {ï}{{\"i}}1 {ö}{{\"o}}1 {ü}{{\"u}}1
  {Ä}{{\"A}}1 {Ë}{{\"E}}1 {Ï}{{\"I}}1 {Ö}{{\"O}}1 {Ü}{{\"U}}1
  {â}{{\^a}}1 {ê}{{\^e}}1 {î}{{\^i}}1 {ô}{{\^o}}1 {û}{{\^u}}1
  {Â}{{\^A}}1 {Ê}{{\^E}}1 {Î}{{\^I}}1 {Ô}{{\^O}}1 {Û}{{\^U}}1
  {œ}{{\oe}}1 {Œ}{{\OE}}1 {æ}{{\ae}}1 {Æ}{{\AE}}1 {ß}{{\ss}}1
  {ç}{{\c c}}1 {Ç}{{\c C}}1 {ø}{{\o}}1 {å}{{\r a}}1 {Å}{{\r A}}1
  {€}{{\EUR}}1 {£}{{\pounds}}1 {«}{{\guillemotleft}}1
  {»}{{\guillemotright}}1 {ñ}{{\~n}}1 {Ñ}{{\~N}}1 {¿}{{?`}}1
  {…}{{\ldots}}1 {≥}{{>=}}1 {≤}{{<=}}1 {„}{{\glqq}}1 {“}{{\grqq}}1
  {”}{{''}}1
}
\lstset{style=eisvogel_listing_style}

%
% Java (Java SE 12, 2019-06-22)
%
\lstdefinelanguage{Java}{
  morekeywords={
    % normal keywords (without data types)
    abstract,assert,break,case,catch,class,continue,default,
    do,else,enum,exports,extends,final,finally,for,if,implements,
    import,instanceof,interface,module,native,new,package,private,
    protected,public,requires,return,static,strictfp,super,switch,
    synchronized,this,throw,throws,transient,try,volatile,while,
    % var is an identifier
    var
  },
  morekeywords={[2] % data types
    % primitive data types
    boolean,byte,char,double,float,int,long,short,
    % String
    String,
    % primitive wrapper types
    Boolean,Byte,Character,Double,Float,Integer,Long,Short
    % number types
    Number,AtomicInteger,AtomicLong,BigDecimal,BigInteger,DoubleAccumulator,DoubleAdder,LongAccumulator,LongAdder,Short,
    % other
    Object,Void,void
  },
  morekeywords={[3] % literals
    % reserved words for literal values
    null,true,false,
  },
  sensitive,
  morecomment  = [l]//,
  morecomment  = [s]{/*}{*/},
  morecomment  = [s]{/**}{*/},
  morestring   = [b]",
  morestring   = [b]',
}

\lstdefinelanguage{XML}{
  morestring      = [b]",
  moredelim       = [s][\bfseries\color{listing-keyword}]{<}{\ },
  moredelim       = [s][\bfseries\color{listing-keyword}]{</}{>},
  moredelim       = [l][\bfseries\color{listing-keyword}]{/>},
  moredelim       = [l][\bfseries\color{listing-keyword}]{>},
  morecomment     = [s]{<?}{?>},
  morecomment     = [s]{<!--}{-->},
  commentstyle    = \color{listing-comment},
  stringstyle     = \color{listing-string},
  identifierstyle = \color{listing-identifier}
}
$endif$

%
% header and footer
%
$if(beamer)$
$else$
$if(disable-header-and-footer)$
$else$
\usepackage[headsepline,footsepline]{scrlayer-scrpage}

\newpairofpagestyles{eisvogel-header-footer}{
  \clearpairofpagestyles
  \ihead*{$if(header-left)$$header-left$$else$$title$$endif$}
  \chead*{$if(header-center)$$header-center$$else$$endif$}
  \ohead*{$if(header-right)$$header-right$$else$$date$$endif$}
  \ifoot*{$if(footer-left)$$footer-left$$else$$for(author)$$author$$sep$, $endfor$$endif$}
  \cfoot*{$if(footer-center)$$footer-center$$else$$endif$}
  \ofoot*{$if(footer-right)$$footer-right$$else$\thepage$endif$}
  \addtokomafont{pageheadfoot}{\upshape}
}
\pagestyle{eisvogel-header-footer}

$if(book)$
\deftripstyle{ChapterStyle}{}{}{}{}{\pagemark}{}
\renewcommand*{\chapterpagestyle}{ChapterStyle}
$endif$

$if(page-background)$
\backgroundsetup{
scale=1,
color=black,
opacity=$if(page-background-opacity)$$page-background-opacity$$else$0.2$endif$,
angle=0,
contents={%
  \includegraphics[width=\paperwidth,height=\paperheight]{$page-background$}
  }%
}
$endif$
$endif$
$endif$

%%
%% end added
%%

\begin{document}

%%
%% begin titlepage
%%
$if(beamer)$
$else$
$if(titlepage)$
\begin{titlepage}
$if(titlepage-background)$
\newgeometry{top=2cm, right=4cm, bottom=3cm, left=4cm}
$else$
\newgeometry{left=6cm}
$endif$
$if(titlepage-color)$
\definecolor{titlepage-color}{HTML}{$titlepage-color$}
\newpagecolor{titlepage-color}\afterpage{\restorepagecolor}
$endif$
$if(titlepage-background)$
\tikz[remember picture,overlay] \node[inner sep=0pt] at (current page.center){\includegraphics[width=\paperwidth,height=\paperheight]{$titlepage-background$}};
$endif$
\newcommand{\colorRule}[3][black]{\textcolor[HTML]{#1}{\rule{#2}{#3}}}
\begin{flushleft}
\noindent
\\[-1em]
\color[HTML]{$if(titlepage-text-color)$$titlepage-text-color$$else$5F5F5F$endif$}
\makebox[0pt][l]{\colorRule[$if(titlepage-rule-color)$$titlepage-rule-color$$else$435488$endif$]{1.3\textwidth}{$if(titlepage-rule-height)$$titlepage-rule-height$$else$4$endif$pt}}
\par
\noindent

$if(titlepage-background)$
% The titlepage with a background image has other text spacing and text size
{
  \setstretch{2}
  \vfill
  \vskip -8em
  \noindent {\huge \textbf{\textsf{$title$}}}
  $if(subtitle)$
  \vskip 1em
  {\Large \textsf{$subtitle$}}
  $endif$
  \vskip 2em
  \noindent {\Large \textsf{$for(author)$$author$$sep$, $endfor$} \vskip 0.6em \textsf{$date$}}
  \vfill
}
$else$
{
  \setstretch{1.4}
  \vfill
  \noindent {\huge \textbf{\textsf{$title$}}}
  $if(subtitle)$
  \vskip 1em
  {\Large \textsf{$subtitle$}}
  $endif$
  \vskip 2em
  \noindent {\Large \textsf{$for(author)$$author$$sep$, $endfor$}}
  \vfill
}
$endif$

$if(titlepage-logo)$
\noindent
\includegraphics[width=$if(logo-width)$$logo-width$$else$35mm$endif$, left]{$titlepage-logo$}
$endif$

$if(titlepage-background)$
$else$
\textsf{$date$}
$endif$
\end{flushleft}
\end{titlepage}
\restoregeometry
\pagenumbering{arabic}
$endif$
$endif$

%%
%% end titlepage
%%

$if(has-frontmatter)$
\frontmatter
$endif$
$if(title)$
$if(beamer)$
\frame{\titlepage}
% don't generate the default title
% $else$
% \maketitle
$endif$
$if(abstract)$
\begin{abstract}
$abstract$
\end{abstract}
$endif$
$endif$

$if(first-chapter)$
\setcounter{chapter}{$first-chapter$}
\addtocounter{chapter}{-1}
$endif$

$for(include-before)$
$include-before$

$endfor$
$if(toc)$
$if(toc-title)$
\renewcommand*\contentsname{$toc-title$}
$endif$
$if(beamer)$
\begin{frame}[allowframebreaks]
$if(toc-title)$
  \frametitle{$toc-title$}
$endif$
  \setcounter{tocdepth}{$toc-depth$}
  \tableofcontents
\end{frame}
$if(toc-own-page)$
\newpage
$endif$
$else$
{
$if(colorlinks)$
\hypersetup{linkcolor=$if(toccolor)$$toccolor$$else$$endif$}
$endif$
\setcounter{tocdepth}{$toc-depth$}
\tableofcontents
$if(toc-own-page)$
\newpage
$endif$
}
$endif$
$endif$
$if(lof)$
\listoffigures
$endif$
$if(lot)$
\listoftables
$endif$
$if(linestretch)$
\setstretch{$linestretch$}
$endif$
$if(has-frontmatter)$
\mainmatter
$endif$
$body$

$if(has-frontmatter)$
\backmatter
$endif$
$if(natbib)$
$if(bibliography)$
$if(biblio-title)$
$if(has-chapters)$
\renewcommand\bibname{$biblio-title$}
$else$
\renewcommand\refname{$biblio-title$}
$endif$
$endif$
$if(beamer)$
\begin{frame}[allowframebreaks]{$biblio-title$}
  \bibliographytrue
$endif$
  \bibliography{$for(bibliography)$$bibliography$$sep$,$endfor$}
$if(beamer)$
\end{frame}
$endif$

$endif$
$endif$
$if(biblatex)$
$if(beamer)$
\begin{frame}[allowframebreaks]{$biblio-title$}
  \bibliographytrue
  \printbibliography[heading=none]
\end{frame}
$else$
\printbibliography$if(biblio-title)$[title=$biblio-title$]$endif$
$endif$

$endif$
$for(include-after)$
$include-after$

$endfor$
\end{document}
EOF

# Set permissions
chmod 644 "$TEMPLATE_FILE"

# Generate PDF with pandoc using container
readonly PDF_FILE="${RESULTS_INTEGRATION_DIR}/StrainCascade_run_documentation${SAMPLE_NAME}.pdf"

apptainer exec --bind "${RESULTS_INTEGRATION_DIR}:/input_dir" "$straincascade_document_processing_sif" \
    bash -c "
        # Create temp dir in /tmp (which should always exist in container)
        TEMP_DIR=/tmp/pandoc_\$\$ && \
        mkdir -p \"\${TEMP_DIR}\" && \
        if [ ! -d \"\${TEMP_DIR}\" ]; then
            echo \"Failed to create temporary directory\" >&2
            exit 1
        fi && \
        chmod 777 \"\${TEMP_DIR}\" && \
        
        # Set TMPDIR environment variable for pandoc
        export TMPDIR=\"\${TEMP_DIR}\" && \
        
        # Run pandoc using the temporary directory
        pandoc \"/input_dir/$(basename "${MD_FILE}")\" \
            -f markdown \
            -t pdf \
            --pdf-engine=xelatex \
            --template=\"/input_dir/$(basename "${TEMPLATE_FILE}")\" \
            -V mainfont=\"Source Sans 3\" \
            -V mainfontoptions=\"Path=/usr/share/fonts/opentype/source-sans/OTF/,Extension=.otf,UprightFont=SourceSans3-Regular,BoldFont=SourceSans3-Bold,ItalicFont=SourceSans3-It,BoldItalicFont=SourceSans3-BoldIt\" \
            -V sansfont=\"Source Sans 3\" \
            -V sansfontoptions=\"Path=/usr/share/fonts/opentype/source-sans/OTF/,Extension=.otf,UprightFont=SourceSans3-Regular,BoldFont=SourceSans3-Bold,ItalicFont=SourceSans3-It,BoldItalicFont=SourceSans3-BoldIt\" \
            --output=\"\${TEMP_DIR}/output.pdf\" && \
            
        # Copy result to mounted directory
        cp \"\${TEMP_DIR}/output.pdf\" \"/input_dir/$(basename "${PDF_FILE}")\" && \
        
        # Cleanup
        rm -rf \"\${TEMP_DIR}\"
    "
    
# Encrypt PDF using container
apptainer exec --bind "${RESULTS_INTEGRATION_DIR}:/input_dir" "$straincascade_document_processing_sif" \
    bash -c '
        # Check if source PDF exists
        if [ ! -f "/input_dir/$(basename '"$PDF_FILE"')" ]; then
            echo "Error: Source PDF not found" >&2
            exit 1
        fi

        # Run encryption with error checking
        qpdf --encrypt "" "$(openssl rand -base64 32)" 256 \
            --print=full \
            --modify=none \
            --extract=y \
            --accessibility=y \
            -- "/input_dir/$(basename '"$PDF_FILE"')" "/input_dir/$(basename '"$PDF_FILE"').encrypted" || {
                echo "Error: PDF encryption failed" >&2
                exit 1
            }

        # Verify encrypted file exists before moving
        if [ -f "/input_dir/$(basename '"$PDF_FILE"').encrypted" ]; then
            mv "/input_dir/$(basename '"$PDF_FILE"').encrypted" "/input_dir/$(basename '"$PDF_FILE"')"
        else
            echo "Error: Encrypted PDF not created" >&2
            exit 1
        fi
    '

# Clean up the temporary file
rm "$MD_FILE"
rm "$TEMPLATE_FILE"

# Generate verification QR code using container
readonly QR_FILE="${RESULTS_INTEGRATION_DIR}/verification_qr.png"
apptainer exec --bind "${RESULTS_INTEGRATION_DIR}:/input_dir" "$straincascade_document_processing_sif" \
    qrencode -o "/input_dir/$(basename "$QR_FILE")" -l H "$(sha256sum "$PDF_FILE" | cut -d' ' -f1)"

# Final output
log "$LOGS_DIR" "$LOG_NAME" "Documentation generation completed successfully:"
log "$LOGS_DIR" "$LOG_NAME" "- PDF: $PDF_FILE"
log "$LOGS_DIR" "$LOG_NAME" "- QR Code: $QR_FILE"

