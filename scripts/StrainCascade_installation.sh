#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_installation.sh
# Description: Installation script for StrainCascade.
#
# Supports both interactive and non-interactive operation. See --help.

set -euo pipefail

# Constants
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly INSTALLATION_UTILS_FILE="$SCRIPT_DIR/StrainCascade_installation_utils.sh"

# Validate environment before sourcing utilities
if [[ ! -f "$INSTALLATION_UTILS_FILE" ]]; then
    echo "Error: Please launch the installation from within the scripts directory"
    exit 1
fi
cd "$SCRIPT_DIR" || { echo "Error: Failed to change to the script's directory"; exit 1; }

# Source utility functions (also defines vlog, retry_with_backoff, ask_yes_no, etc.)
# shellcheck disable=SC1090
source "$INSTALLATION_UTILS_FILE"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
print_install_help() {
    cat <<'EOF'
StrainCascade installation script

Usage:
  ./scripts/StrainCascade_installation.sh [OPTIONS]

Phase selection (combine freely; default = interactive prompts):
  --full              Run all phases (scripts PATH + images + databases)
  --scripts           Add scripts directory to PATH only
  --images            Pull Apptainer images only
  --databases         Install databases only

Behaviour:
  --skip-existing     Skip already-installed images/databases (default for --full)
  --force             Force reinstall (overrides --skip-existing)
  -y, --yes           Non-interactive: auto-accept all prompts
  -h, --help          Show this help

Verbosity is controlled by the STRAINCASCADE_VERBOSITY environment variable:
  0 = quiet     (only errors and section headers; output goes to log file)
  1 = normal    (default)
  2 = verbose   (full subcommand output on terminal)

Examples:
  # Interactive installation (recommended for first-time users)
  ./scripts/StrainCascade_installation.sh

  # Non-interactive full install, idempotent (safe to re-run after a crash)
  ./scripts/StrainCascade_installation.sh --full --yes

  # Force reinstall everything
  ./scripts/StrainCascade_installation.sh --full --yes --force

  # Quiet full install with logged output
  STRAINCASCADE_VERBOSITY=0 ./scripts/StrainCascade_installation.sh --full --yes
EOF
}

run_full=0
run_scripts_phase=0
run_images_phase=0
run_databases_phase=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --full)            run_full=1; shift ;;
        --scripts)         run_scripts_phase=1; shift ;;
        --images)          run_images_phase=1; shift ;;
        --databases)       run_databases_phase=1; shift ;;
        --skip-existing)   SC_SKIP_EXISTING=1; shift ;;
        --force)           SC_FORCE=1; SC_SKIP_EXISTING=0; shift ;;
        -y|--yes)          SC_NON_INTERACTIVE=1; shift ;;
        -h|--help)         print_install_help; exit 0 ;;
        *)
            echo "Error: unknown option '$1'" >&2
            echo "Try './scripts/StrainCascade_installation.sh --help' for usage." >&2
            exit 2
            ;;
    esac
done

# --full implies --skip-existing unless the user explicitly asked for --force.
# This makes a re-run after a transient failure idempotent.
if [ "$run_full" = "1" ] && [ "$SC_FORCE" != "1" ]; then
    SC_SKIP_EXISTING=1
fi

export STRAINCASCADE_VERBOSITY SC_NON_INTERACTIVE SC_SKIP_EXISTING SC_FORCE

# ---------------------------------------------------------------------------
# Pre-flight checks (always run)
# ---------------------------------------------------------------------------
check_main_directory
check_disk_space
check_required_tools
check_apptainer_installation
check_existing_installations

setup_apptainer_env
trap 'cleanup_apptainer_env' EXIT

# ---------------------------------------------------------------------------
# Phase dispatch
# ---------------------------------------------------------------------------
any_phase_flag=$(( run_full + run_scripts_phase + run_images_phase + run_databases_phase ))

if [ "$any_phase_flag" -eq 0 ]; then
    # Legacy interactive flow (no flags supplied)
    read -r -p "Do you want to do a full installation of StrainCascade? This will overwrite any existing installations unless --skip-existing was supplied. (y/n) " full_install_choice
    case "$full_install_choice" in
        y|Y)
            run_full=1
            ;;
        n|N)
            read -r -p "Do you want to update the StrainCascade scripts? (y/n) " update_choice
            case "$update_choice" in
                y|Y)
                    add_scripts_to_path
                    update_scripts
                    ;;
                *) echo "Proceeding with the current version of scripts." ;;
            esac

            read -r -p "Do you want to install Apptainer images? (y/n) " install_apptainer_choice
            case "$install_apptainer_choice" in
                y|Y) run_images_phase=1 ;;
                *) echo "Skipping Apptainer image installation." ;;
            esac

            read -r -p "Do you want to install databases? (y/n) " install_db_choice
            case "$install_db_choice" in
                y|Y) run_databases_phase=1 ;;
                *) echo "Skipping database installation." ;;
            esac
            ;;
        *) echo "Invalid choice. Exiting."; exit 1 ;;
    esac
fi

if [ "$run_full" = "1" ]; then
    vlog 0 "Proceeding with full installation..."
    vlog 1 "Apptainer images will be installed in: $APPTAINER_IMAGES_DIR"
    vlog 1 "Databases will be installed in: $DEFAULT_DB_LOCATION"
    if [ "$SC_SKIP_EXISTING" = "1" ]; then
        vlog 1 "Mode: --skip-existing (already-installed components will be skipped)"
    elif [ "$SC_FORCE" = "1" ]; then
        vlog 1 "Mode: --force (existing components will be overwritten)"
    fi

    add_scripts_to_path
    pull_apptainer_images

    mkdir -p "$DEFAULT_DB_LOCATION"
    if ! get_apptainer_images "$APPTAINER_IMAGES_DIR"; then
        vlog 0 "Some Apptainer images were not found. Continuing with database installation."
    fi
    vlog 0 "Installing all databases..."
    install_all_databases "$DEFAULT_DB_LOCATION"
else
    if [ "$run_scripts_phase" = "1" ]; then
        add_scripts_to_path
    fi
    if [ "$run_images_phase" = "1" ]; then
        pull_apptainer_images
        # Re-add to PATH in case scripts dir wasn't yet picked up
        add_scripts_to_path
    fi
    if [ "$run_databases_phase" = "1" ]; then
        if ! get_apptainer_images "$APPTAINER_IMAGES_DIR"; then
            vlog 0 "Some Apptainer images were not found. Continuing without database installation."
        else
            mkdir -p "$DEFAULT_DB_LOCATION"
            if [ "$SC_NON_INTERACTIVE" = "1" ]; then
                install_all_databases "$DEFAULT_DB_LOCATION"
            else
                install_databases "$DEFAULT_DB_LOCATION"
            fi
            add_scripts_to_path
        fi
    fi
fi

vlog 0 "Installation process completed"
