#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_installation.sh
# Description: Installation script for StrainCascade

set -euo pipefail

# Constants
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly INSTALLATION_UTILS_FILE="$SCRIPT_DIR/StrainCascade_installation_utils.sh"

# Main script execution
source "$INSTALLATION_UTILS_FILE"

# Check environment setup
check_main_directory() {
    if [[ ! -f "$INSTALLATION_UTILS_FILE" ]]; then
        echo "Error: Please launch the installation from within the scripts directory"
        exit 1
    fi
    cd "$SCRIPT_DIR" || { echo "Error: Failed to change to the script's directory"; exit 1; }
}

# Execute checks
check_main_directory
check_disk_space
check_apptainer_installation
check_existing_installations

# Full or individual component installation
read -p "Do you want to do a full installation of StrainCascade? This will overwrite any existing installations. (y/n) " full_install_choice
case "$full_install_choice" in
    y|Y)
        echo "Proceeding with full installation..."
        echo "Apptainer images will be installed in: $APPTAINER_IMAGES_DIR"
        echo "Databases will be installed in: $DEFAULT_DB_LOCATION"
        
        pull_apptainer_images
        add_scripts_to_path

        mkdir -p "$DEFAULT_DB_LOCATION"
        if ! get_apptainer_images "$APPTAINER_IMAGES_DIR"; then
            echo "Some Apptainer images were not found. Continuing with database installation."
        fi
        echo "Installing all databases..."
        install_all_databases "$DEFAULT_DB_LOCATION"
        ;;
    n|N)
        read -p "Do you want to update the StrainCascade scripts? (y/n) " update_choice
        case "$update_choice" in
            y|Y)
                add_scripts_to_path
                update_scripts
                ;;
            n|N) echo "Proceeding with the current version of scripts." ;;
            *) echo "Invalid choice. Proceeding with the current version of scripts." ;;
        esac

        read -p "Do you want to install Apptainer images? (y/n) " install_apptainer_choice
        case "$install_apptainer_choice" in
            y|Y)
                pull_apptainer_images
                add_scripts_to_path
                ;;
            n|N) echo "Skipping Apptainer image installation." ;;
            *) echo "Invalid choice. Skipping Apptainer image installation." ;;
        esac

        read -p "Do you want to install databases? (y/n) " install_db_choice
        case "$install_db_choice" in
            y|Y)
                if ! get_apptainer_images "$APPTAINER_IMAGES_DIR"; then
                    echo "Some Apptainer images were not found. Continuing without database installation."
                else
                    mkdir -p "$DEFAULT_DB_LOCATION"
                    install_databases "$DEFAULT_DB_LOCATION"
                    add_scripts_to_path
                fi
                ;;
            n|N) echo "Skipping database installation." ;;
            *) echo "Invalid choice. Skipping database installation." ;;
        esac
        ;;
    *) echo "Invalid choice. Exiting."; exit 1 ;;
esac

echo "Installation process completed"
