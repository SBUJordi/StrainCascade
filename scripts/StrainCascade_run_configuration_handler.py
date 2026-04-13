#!/usr/bin/env python3

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_run_configuration_handler.py

import sys
import difflib
import re

# Dictionary mapping script names to human-readable descriptions with SC prefix
script_info = {
    "StrainCascade_Canu_correct_trim.sh": "SC1 Canu Correction and Trimming",
    "StrainCascade_LJA_assembly.sh": "SC2 LJA Assembly",
    "StrainCascade_SPAdes_assembly.sh": "SC3 SPAdes Assembly",
    "StrainCascade_Canu_assembly.sh": "SC4 Canu Assembly",
    "StrainCascade_Flye_assembly.sh": "SC5 Flye Assembly",
    "StrainCascade_Unicycler_assembly.sh": "SC6 Unicycler Assembly",
    "StrainCascade_assembly_evaluation1.sh": "SC7 Assembly Evaluation 1",
    "StrainCascade_MAC2_assembly_merging.sh": "SC8 MAC2 Assembly Merging",
    "StrainCascade_assembly_evaluation2.sh": "SC9 Assembly Evaluation 2",
    "StrainCascade_Circlator_circularisation.sh": "SC10 Circlator Circularisation",
    "StrainCascade_assembly_evaluation3.sh": "SC11 Assembly Evaluation 3",
    "StrainCascade_assembly_polishing.sh": "SC12 Assembly Polishing",
    "StrainCascade_minimap2_BBMap_coverage.sh": "SC13 minimap2 BBMap Coverage",
    "StrainCascade_CheckM2_QC.sh": "SC14 CheckM2 QC",
    "StrainCascade_GTDB-Tk_taxonomy.sh": "SC15 GTDB-Tk Taxonomy",
    "StrainCascade_GTDB-Tk_de_novo_tree.sh": "SC16 GTDB-Tk De Novo Tree",
    "StrainCascade_Bakta_annotation.sh": "SC17 Bakta Annotation",
    "StrainCascade_Prokka_annotation.sh": "SC18 Prokka Annotation",
    "StrainCascade_DeepFRI_annotation.sh": "SC19 DeepFRI Annotation",
    "StrainCascade_MicrobeAnnotator_annotation.sh": "SC20 MicrobeAnnotator Annotation",
    "StrainCascade_PlasmidFinder_identification.sh": "SC21 PlasmidFinder Identification",
    "StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh": "SC22 AMRFinderPlus Antimicrobial Resistance Identification",
    "StrainCascade_ResFinder_antimicrobial_resistance_identification.sh": "SC23 ResFinder Antimicrobial Resistance Identification",
    "StrainCascade_dbCAN3_CAZymes_identification.sh": "SC24 dbCAN3 CAZymes Identification",
    "StrainCascade_IslandPath_genomic_islands_identification.sh": "SC25 IslandPath Genomic Islands Identification",
    "StrainCascade_VirSorter2_phage_identification.sh": "SC26 VirSorter2 Phage Identification",
    "StrainCascade_geNomad_phage_identification.sh": "SC27 geNomad Phage and Plasmid Identification",
    "StrainCascade_CRISPRCasFinder_identification.sh": "SC28 CRISPRCasFinder CRISPRCas Identification",
    "StrainCascade_ISEScan_IS_elements_identification.sh": "SC29 ISEScan IS Elements Identification",
    "StrainCascade_data_integration.sh": "SC30 Data Integration"
}

# Reverse dictionary for human-readable names and SC numbers
script_lookup = {v.split(" ", 1)[1].lower(): k for k, v in script_info.items()}
script_number_lookup = {v.split(" ")[0]: k for k, v in script_info.items()}

# Predefined bundles and execution modes
bundles = {
    "assembly": [
        "StrainCascade_Canu_correct_trim.sh",
        "StrainCascade_LJA_assembly.sh",
        "StrainCascade_SPAdes_assembly.sh",
        "StrainCascade_Canu_assembly.sh",
        "StrainCascade_Flye_assembly.sh",
        "StrainCascade_Unicycler_assembly.sh",
        "StrainCascade_assembly_evaluation1.sh",
        "StrainCascade_MAC2_assembly_merging.sh",
        "StrainCascade_assembly_evaluation2.sh",
        "StrainCascade_Circlator_circularisation.sh",
        "StrainCascade_assembly_evaluation3.sh",
        "StrainCascade_assembly_polishing.sh",
        "StrainCascade_minimap2_BBMap_coverage.sh",
        "StrainCascade_CheckM2_QC.sh"
    ],
    "annotation": [
        "StrainCascade_Bakta_annotation.sh",
        "StrainCascade_Prokka_annotation.sh",
        "StrainCascade_DeepFRI_annotation.sh",
        "StrainCascade_MicrobeAnnotator_annotation.sh"
    ],
    "functional": [
        "StrainCascade_GTDB-Tk_taxonomy.sh",
        "StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh",
        "StrainCascade_ResFinder_antimicrobial_resistance_identification.sh",
        "StrainCascade_dbCAN3_CAZymes_identification.sh",
        "StrainCascade_IslandPath_genomic_islands_identification.sh"
    ],
    "phage": [
        "StrainCascade_VirSorter2_phage_identification.sh",
        "StrainCascade_geNomad_phage_identification.sh",
        "StrainCascade_CRISPRCasFinder_identification.sh",
        "StrainCascade_ISEScan_IS_elements_identification.sh"
    ]
}

execution_modes = {
    "minimal": [
        "StrainCascade_SPAdes_assembly.sh",
        "StrainCascade_assembly_evaluation1.sh",
        "StrainCascade_CheckM2_QC.sh",
        "StrainCascade_GTDB-Tk_taxonomy.sh",
        "StrainCascade_Bakta_annotation.sh",
        "StrainCascade_data_integration.sh"
    ],
    "efficient": [
        "StrainCascade_Canu_correct_trim.sh",
        "StrainCascade_LJA_assembly.sh",
        "StrainCascade_SPAdes_assembly.sh",
        "StrainCascade_assembly_evaluation1.sh",
        "StrainCascade_MAC2_assembly_merging.sh",
        "StrainCascade_assembly_evaluation2.sh",
        "StrainCascade_Circlator_circularisation.sh",
        "StrainCascade_assembly_evaluation3.sh",
        "StrainCascade_assembly_polishing.sh",
        "StrainCascade_CheckM2_QC.sh",
        "StrainCascade_GTDB-Tk_taxonomy.sh",
        "StrainCascade_Bakta_annotation.sh",
        "StrainCascade_PlasmidFinder_identification.sh",
        "StrainCascade_data_integration.sh"
    ],
    "standard": [
        "StrainCascade_Canu_correct_trim.sh",
        "StrainCascade_LJA_assembly.sh",
        "StrainCascade_SPAdes_assembly.sh",
        "StrainCascade_Canu_assembly.sh",
        "StrainCascade_Flye_assembly.sh",
        "StrainCascade_Unicycler_assembly.sh",
        "StrainCascade_assembly_evaluation1.sh",
        "StrainCascade_MAC2_assembly_merging.sh",
        "StrainCascade_assembly_evaluation2.sh",
        "StrainCascade_Circlator_circularisation.sh",
        "StrainCascade_assembly_evaluation3.sh",
        "StrainCascade_assembly_polishing.sh",
        "StrainCascade_minimap2_BBMap_coverage.sh",
        "StrainCascade_CheckM2_QC.sh",
        "StrainCascade_GTDB-Tk_taxonomy.sh",
        "StrainCascade_Bakta_annotation.sh",
        "StrainCascade_Prokka_annotation.sh",
        "StrainCascade_DeepFRI_annotation.sh",
        "StrainCascade_MicrobeAnnotator_annotation.sh",
        "StrainCascade_PlasmidFinder_identification.sh",
        "StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh",
        "StrainCascade_ResFinder_antimicrobial_resistance_identification.sh",
        "StrainCascade_dbCAN3_CAZymes_identification.sh",
        "StrainCascade_IslandPath_genomic_islands_identification.sh",
        "StrainCascade_VirSorter2_phage_identification.sh",
        "StrainCascade_geNomad_phage_identification.sh",
        "StrainCascade_CRISPRCasFinder_identification.sh",
        "StrainCascade_ISEScan_IS_elements_identification.sh",
        "StrainCascade_data_integration.sh"
    ],
    "comprehensive": list(script_info.keys())
}

def find_best_match(input_string, script_info):
    input_string = input_string.lower().strip()

    # Check if input is a direct match with human-readable names or SC numbers
    if input_string in script_lookup:
        return [script_lookup[input_string]]
    elif input_string.isdigit() and f"SC{input_string}" in script_number_lookup:
        return [script_number_lookup[f"SC{input_string}"]]

    # Fuzzy matching
    all_scripts = list(script_info.keys())
    best_match = difflib.get_close_matches(input_string, all_scripts, n=1, cutoff=0.6)

    return best_match if best_match else None

def preprocess_input(input_string):
    # Remove spaces and convert to uppercase
    input_string = input_string.strip().upper()
    
    # Check for SCx or SCxx pattern
    if re.search(r'SC\d{1,2}$', input_string):
        return input_string
    
    # Check for pure number
    if input_string.isdigit():
        return f"SC{input_string}"
    
    # Return as is for text inputs
    return input_string

def process_custom_input(inputs, script_info):
    matched_scripts = []
    seen = {}
    info_messages = []

    for input_string in inputs:
        processed_input = preprocess_input(input_string)
        
        if processed_input.startswith('SC') and re.fullmatch(r'SC\d+', processed_input):
            # Direct lookup for SC numbers (e.g., SC7, SC12)
            script_number = processed_input
            if script_number in script_number_lookup:
                script = script_number_lookup[script_number]
                matched_scripts.append(script)
                seen[script] = seen.get(script, 0) + 1
            else:
                info_messages.append(f"Warning: No match found for {input_string}")
        else:
            # Use existing fuzzy matching for text inputs
            matches = find_best_match(processed_input, script_info)
            if matches:
                script = matches[0]
                matched_scripts.append(script)
                seen[script] = seen.get(script, 0) + 1
            else:
                info_messages.append(f"Warning: No match found for {input_string}")

    # Check for duplicates
    for script, count in seen.items():
        if count > 1:
            info_messages.append(f"Warning: Module '{script_info[script]}' appears {count} times.")

    # Sort scripts by SC number
    def sc_sort_key(x):
        sc_label = script_info[x].split(" ")[0][2:]  # e.g., '6', '12'
        return int(sc_label) if sc_label.isdigit() else 0
    sorted_scripts = sorted(matched_scripts, key=sc_sort_key)

    if matched_scripts != sorted_scripts:
        info_messages.append("Note: Script numbers were not in increasing order. They have been rearranged.")

    return info_messages, sorted_scripts

def main():
    if len(sys.argv) != 4:
        print("Error: Incorrect number of arguments.")
        sys.exit(1)

    execution_mode = sys.argv[1]
    bundle = sys.argv[2]
    custom_args = [arg.strip() for arg in sys.argv[3].split(',')] if sys.argv[3] else []

    selected_modules = []
    info_messages = []

    if execution_mode and bundle:
        print("Error: Execution mode and bundle are mutually exclusive. Please specify only one.")
        sys.exit(1)

    if execution_mode:
        if execution_mode in execution_modes:
            selected_modules = execution_modes[execution_mode]
            info_messages.append(f"Selected execution mode: {execution_mode}")
        elif execution_mode == "custom":
            if not custom_args:
                print("Error: Custom execution mode requires arguments.")
                sys.exit(1)
            info_messages, selected_modules = process_custom_input(custom_args, script_info)
            info_messages.insert(0, "Selected custom execution mode")
        else:
            print(f"Error: Unknown execution mode '{execution_mode}'.")
            sys.exit(1)
    elif bundle:
        if bundle in bundles:
            selected_modules = bundles[bundle]
            info_messages.append(f"Selected bundle: {bundle}")
        else:
            print(f"Error: Unknown bundle '{bundle}'.")
            sys.exit(1)

    else:
        print("Error: Neither execution mode nor bundle specified.")
        sys.exit(1)

    # Prepare output strings
    info_output = ",".join(info_messages)
    scripts_output = " ".join(selected_modules)
    
    # Create a new list with human-readable names
    human_readable_output = [script_info[module] for module in selected_modules]
    human_readable_string = ", ".join(human_readable_output)
    
    # Print outputs (which can be captured separately by the shell script)
    print(info_output)
    print(scripts_output)
    print(human_readable_string)

if __name__ == "__main__":
    main()