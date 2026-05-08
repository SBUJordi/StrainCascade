#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_installation_utils.sh
# Description: Utility functions for StrainCascade installation

# Define required directories (only scripts must pre-exist; others are created during installation)
required_dirs=("scripts")

# Define minimum disk space in GB (last measurement: 534G)
min_disk_space=580

# Default database location
readonly DEFAULT_DB_LOCATION="../databases"

# Default Apptainer images directory
readonly APPTAINER_IMAGES_DIR="../apptainer_images"

# ----------------------------------------------------------------------------
# Installer behaviour flags (set by StrainCascade_installation.sh argument
# parser; safe defaults so functions also work when called from other scripts).
#   STRAINCASCADE_VERBOSITY  0=quiet  1=normal (default)  2=verbose
#   SC_NON_INTERACTIVE       1 to suppress all prompts (auto-yes)
#   SC_SKIP_EXISTING         1 to skip already-installed images/databases
#   SC_FORCE                 1 to force reinstall (overrides SC_SKIP_EXISTING)
# ----------------------------------------------------------------------------
: "${STRAINCASCADE_VERBOSITY:=1}"
: "${SC_NON_INTERACTIVE:=0}"
: "${SC_SKIP_EXISTING:=0}"
: "${SC_FORCE:=0}"
export STRAINCASCADE_VERBOSITY SC_NON_INTERACTIVE SC_SKIP_EXISTING SC_FORCE

# Verbosity-aware logging.
# Levels: 0=always (errors/section headers), 1=normal info, 2=debug.
vlog() {
    local level="$1"; shift
    if [ "$STRAINCASCADE_VERBOSITY" -ge "$level" ]; then
        echo "$@"
    fi
}

# Run a command, routing its stdout/stderr according to verbosity.
#   Verbosity 0: stdout+stderr -> log file only
#   Verbosity 1: stdout -> log file, stderr -> terminal + log file
#   Verbosity 2: stdout+stderr -> terminal + log file
# Returns the command's exit code.
sc_install_log_dir() {
    local img_dir
    img_dir="$(cd "$APPTAINER_IMAGES_DIR" 2>/dev/null && pwd)" || img_dir="$APPTAINER_IMAGES_DIR"
    local log_dir="$img_dir/.install_logs"
    mkdir -p "$log_dir" 2>/dev/null || true
    echo "$log_dir"
}

vrun() {
    local log_dir
    log_dir="$(sc_install_log_dir)"
    local log_file="$log_dir/install_$(date +%Y%m%d).log"
    {
        echo "----- $(date '+%Y-%m-%d %H:%M:%S') CMD: $* -----"
    } >> "$log_file" 2>/dev/null || true

    if [ "$STRAINCASCADE_VERBOSITY" -ge 2 ]; then
        # Verbose: tee everything to terminal AND log
        "$@" 2>&1 | tee -a "$log_file"
        return "${PIPESTATUS[0]}"
    elif [ "$STRAINCASCADE_VERBOSITY" -le 0 ]; then
        # Quiet: log only
        "$@" >> "$log_file" 2>&1
        return $?
    else
        # Normal: stderr to terminal + log, stdout to log only
        "$@" 2> >(tee -a "$log_file" >&2) >> "$log_file"
        return $?
    fi
}

# Generic retry-with-exponential-backoff wrapper.
# Usage: retry_with_backoff <max_attempts> <initial_sleep_s> -- <cmd> [args...]
# Or:    retry_with_backoff <max_attempts> <initial_sleep_s> <cmd_string>
retry_with_backoff() {
    local max_attempts="$1"; shift
    local sleep_seconds="$1"; shift
    # Allow '--' delimiter for clarity
    if [ "${1:-}" = "--" ]; then shift; fi

    local attempt=1
    local rc=0
    while [ "$attempt" -le "$max_attempts" ]; do
        if [ "$attempt" -gt 1 ]; then
            vlog 1 "  Retry attempt $attempt/$max_attempts (waiting ${sleep_seconds}s)..."
            sleep "$sleep_seconds"
            sleep_seconds=$(( sleep_seconds * 2 ))
        fi
        if [ $# -eq 1 ]; then
            # shellcheck disable=SC2086
            bash -c "$1"
            rc=$?
        else
            "$@"
            rc=$?
        fi
        if [ "$rc" -eq 0 ]; then
            return 0
        fi
        # Fail-fast on signal kills (128 + signal number):
        #   137 = SIGKILL (typically OOM-killer)
        #   143 = SIGTERM (typically SLURM walltime / cancellation)
        # These are deterministic resource/control failures; retrying just
        # wastes hours of compute. Surface them immediately so the caller
        # can raise memory / walltime instead.
        if [ "$rc" -eq 137 ] || [ "$rc" -eq 143 ]; then
            vlog 0 "  Command was killed by a signal (exit $rc); not retrying (likely OOM or walltime)."
            return "$rc"
        fi
        attempt=$(( attempt + 1 ))
    done
    return "$rc"
}

# Prompt the user (y/n) honouring SC_NON_INTERACTIVE.
# Usage: ask_yes_no "Question text" <default y|n>
# Returns 0 for yes, 1 for no.
ask_yes_no() {
    local prompt="$1"
    local default="${2:-n}"
    local hint="[y/N]"
    [ "$default" = "y" ] && hint="[Y/n]"

    if [ "$SC_NON_INTERACTIVE" = "1" ]; then
        # Non-interactive: assume default
        [ "$default" = "y" ] && return 0 || return 1
    fi

    local reply
    read -r -p "$prompt $hint " reply
    reply="${reply:-$default}"
    case "$reply" in
        y|Y|yes|YES) return 0 ;;
        *) return 1 ;;
    esac
}

# Function to check if script is in the main directory
check_main_directory() {
    for dir in "${required_dirs[@]}"; do
        if [ ! -d "../$dir" ]; then
            echo "Error: Directory '../$dir' not found. Please ensure you are running this script from the scripts directory."
            exit 1
        fi
    done
    # Create optional directories if they don't exist yet (fresh installation)
    mkdir -p "../databases" "../apptainer_images"
}

# Function to check disk space (POSIX-compatible, works across Linux distributions)
# Measures the parent directory (../) where images and databases are installed,
# not the scripts directory, to correctly reflect the target filesystem.
check_disk_space() {
    local install_root="$(cd ".." && pwd)"
    local available_kb
    available_kb=$(df -P "$install_root" 2>/dev/null | awk 'NR==2 {print $4}')
    if [[ -z "$available_kb" ]] || ! [[ "$available_kb" =~ ^[0-9]+$ ]]; then
        echo "Warning: Could not determine available disk space. Continuing installation."
        return 0
    fi
    local available_gb=$((available_kb / 1048576))
    if [ "$available_gb" -lt "$min_disk_space" ]; then
        echo "Warning: Insufficient disk space for full StrainCascade installation."
        echo "Available disk space: ${available_gb}GB"
        echo "Recommended minimal disk space: ${min_disk_space}GB"
        read -p "Do you want to continue anyway? (y/n) " choice
        case "$choice" in
            y|Y) echo "Continuing with installation despite low disk space." ;;
            *) echo "Installation aborted due to insufficient disk space."; exit 0 ;;
        esac
    else
        echo "Sufficient disk space available: ${available_gb}GB"
        echo "Recommended minimal disk space: ${min_disk_space}GB"
    fi

    # Warn if the system's default temp directory is on a small/separate filesystem.
    # Apptainer needs ~20 GB in its temp dir for SIF conversions. The pull function
    # redirects APPTAINER_TMPDIR, but this early warning helps users understand
    # potential issues with other tools that still use /tmp.
    local min_tmp_gb=20
    local tmp_dir="${TMPDIR:-/tmp}"
    local tmp_kb
    tmp_kb=$(df -P "$tmp_dir" 2>/dev/null | awk 'NR==2 {print $4}')
    if [[ -n "$tmp_kb" ]] && [[ "$tmp_kb" =~ ^[0-9]+$ ]]; then
        local tmp_gb=$((tmp_kb / 1048576))
        if [ "$tmp_gb" -lt "$min_tmp_gb" ]; then
            echo "Warning: Temp directory ($tmp_dir) has only ${tmp_gb}GB free (recommend >= ${min_tmp_gb}GB)."
            echo "  Apptainer image pulls will use a local temp dir to avoid this constraint."
        fi
    fi
}

# Function to check if Apptainer is installed
check_apptainer_installation() {
    if ! command -v apptainer &> /dev/null; then
        echo "Error: Apptainer is not installed."
        echo "Please install Apptainer and try again."
        exit 1
    else
        echo "Apptainer is installed."
    fi
}

# Function to check if required command-line tools are available
check_required_tools() {
    local missing_tools=""
    local required_tools=("wget" "curl" "git" "tar")

    for tool in "${required_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+="  - $tool\n"
        fi
    done

    if [[ -n "$missing_tools" ]]; then
        echo "Error: The following required tools are not installed:"
        echo -e "$missing_tools"
        echo "Please install them and try again."
        exit 1
    fi
    echo "All required tools are available."
}

# Function to update StrainCascade scripts.
#
# Behaviour:
#   * If the parent directory is a git checkout (has a .git/ folder), perform
#     a fast-forward-only `git pull` against the configured remote. This is
#     the safe, standard way to update a cloned public repo and keeps the
#     git history intact for users who want to track upstream.
#   * Otherwise (e.g. the user installed from a tarball with no .git), fall
#     back to the previous clone-into-temp-and-copy approach.
#
# IMPORTANT: this function NO LONGER deletes the user's .git directory.
update_scripts() {
    local current_dir
    current_dir=$(pwd)
    local parent_dir
    parent_dir="$(dirname "$current_dir")"

    echo "Updating StrainCascade scripts..."

    if [ -d "${parent_dir}/.git" ] && command -v git >/dev/null 2>&1; then
        echo "Detected git checkout at ${parent_dir}; updating via 'git pull --ff-only'..."
        if ( cd "$parent_dir" && git fetch --quiet && git pull --ff-only ); then
            echo "StrainCascade has been updated successfully via git."
            echo "Please restart any running shells / sourced configs for changes to take effect."
            return 0
        else
            echo "Warning: 'git pull --ff-only' failed (likely diverging history or local changes)." >&2
            echo "Refusing to overwrite local modifications. Resolve manually with:" >&2
            echo "    cd ${parent_dir} && git status && git pull --ff-only" >&2
            return 1
        fi
    fi

    # Non-git installation: fall back to clone-and-copy.
    local temp_dir="$parent_dir/script_update_temp"
    mkdir -p "$temp_dir"
    cd "$temp_dir" || return 1

    if ! git clone --depth 1 https://github.com/SBUJordi/StrainCascade.git; then
        echo "Failed to clone the repository. Update aborted." >&2
        cd "$current_dir" || true
        rm -rf "$temp_dir"
        return 1
    fi

    cp -R StrainCascade/scripts/* "$current_dir"
    mkdir -p "${current_dir}/../assets"
    cp -R "StrainCascade/assets/." "${current_dir}/../assets/"

    cd "$current_dir" || true
    rm -rf "$temp_dir"

    echo "StrainCascade scripts have been updated successfully."
    echo "Please restart any running shells / sourced configs for changes to take effect."
    return 0
}

# Redirect Apptainer/Singularity temp and cache directories to the same filesystem
# as the installation target. By default Apptainer uses /tmp for SIF conversion and
# layer extraction, which on many HPC systems is a small tmpfs (RAM-backed).
# This function must be called once before any apptainer pull/exec operations.
setup_apptainer_env() {
    mkdir -p "$APPTAINER_IMAGES_DIR"
    local img_dir
    img_dir="$(cd "$APPTAINER_IMAGES_DIR" && pwd)"
    export APPTAINER_TMPDIR="$img_dir/.apptainer_tmp"
    export APPTAINER_CACHEDIR="$img_dir/.apptainer_cache"
    # Also set legacy Singularity variables for older runtime versions
    export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
    export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
    mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

    # Force TMPDIR=/tmp inside the container. On SLURM clusters,
    # the host TMPDIR is set to a node-local path like /scratch/local/<jobid>
    # which is NOT bind-mounted into the container by default. Apptainer would
    # otherwise propagate that variable into the container, causing tools that
    # honour $TMPDIR (amrfinder_update, mamba/conda, snakemake) to crash with:
    #   "Error creating a temporary directory in /scratch/local/<jobid>"
    # The APPTAINERENV_ / SINGULARITYENV_ prefix translates to setting the
    # variable INSIDE the container only, regardless of the host value.
    export APPTAINERENV_TMPDIR="/tmp"
    export SINGULARITYENV_TMPDIR="/tmp"
}

# Clean up Apptainer temporary build artefacts. Call at the end of all
# apptainer operations (the cache directory is kept for potential re-pulls).
cleanup_apptainer_env() {
    if [[ -n "${APPTAINER_TMPDIR:-}" && -d "${APPTAINER_TMPDIR}" ]]; then
        rm -rf "$APPTAINER_TMPDIR"
    fi
}

# Function to pull Apptainer images
# Honours SC_SKIP_EXISTING / SC_FORCE / SC_NON_INTERACTIVE.
# Failures are collected and retried once at the end; the function does not
# abort the whole installation on a single image failure.
pull_apptainer_images() {
    mkdir -p "$APPTAINER_IMAGES_DIR"
    # Add the list of docker images to be pulled (always use :latest tag)
    docker_images=(
        "sbujordi/straincascade_genome_assembly:latest"
        "sbujordi/straincascade_lja_genome_assembly:latest"
        "sbujordi/straincascade_assembly_qc_refinement:latest"
        "sbujordi/straincascade_genome_annotation:latest"
        "sbujordi/straincascade_taxonomic_functional_analysis:latest"
        "sbujordi/straincascade_crisprcas_phage_is_elements:latest"
        "sbujordi/python_3.12.4:latest"
        "sbujordi/r_4.4.1:latest"
        "sbujordi/straincascade_document_processing:latest"
    )

    local total=${#docker_images[@]}
    local current=0
    local failed_images=()
    local skipped=0
    local pulled=0

    _pull_one_image() {
        local image="$1"
        local image_file="$2"
        local force_flag="$3"
        # apptainer pull writes mostly progress to stderr; route through vrun.
        if [ -n "$force_flag" ]; then
            vrun apptainer pull "$force_flag" "$image_file" "docker://$image"
        else
            vrun apptainer pull "$image_file" "docker://$image"
        fi
    }

    for image in "${docker_images[@]}"; do
        current=$((current + 1))
        local image_basename
        image_basename=$(basename "$image" | tr ':' '_')
        local image_file="$APPTAINER_IMAGES_DIR/${image_basename}.sif"
        local force_flag=""

        if [ -f "$image_file" ]; then
            if [ "$SC_FORCE" = "1" ]; then
                force_flag="--force"
            elif [ "$SC_SKIP_EXISTING" = "1" ]; then
                vlog 1 "[$current/$total] Skipping existing image: $image_basename.sif"
                skipped=$((skipped + 1))
                continue
            elif [ "$SC_NON_INTERACTIVE" = "1" ]; then
                # Default for interactive non-flagged runs: skip if present
                vlog 1 "[$current/$total] Skipping existing image (non-interactive default): $image_basename.sif"
                skipped=$((skipped + 1))
                continue
            else
                if ask_yes_no "Image file $image_file already exists. Overwrite?" "y"; then
                    force_flag="--force"
                else
                    vlog 1 "Skipping $image_file..."
                    skipped=$((skipped + 1))
                    continue
                fi
            fi
        fi

        vlog 0 "[$current/$total] Pulling $image..."
        if retry_with_backoff 3 30 -- _pull_one_image "$image" "$image_file" "$force_flag"; then
            pulled=$((pulled + 1))
        else
            vlog 0 "Warning: Failed to pull $image after 3 attempts; will retry at end of phase." >&2
            failed_images+=("$image")
        fi
    done

    # Final retry pass for any failures (one more attempt with longer backoff).
    if [ "${#failed_images[@]}" -gt 0 ]; then
        vlog 0 "Retrying ${#failed_images[@]} failed image pull(s) once more..."
        local still_failed=()
        for image in "${failed_images[@]}"; do
            local image_basename
            image_basename=$(basename "$image" | tr ':' '_')
            local image_file="$APPTAINER_IMAGES_DIR/${image_basename}.sif"
            local force_flag=""
            [ -f "$image_file" ] && force_flag="--force"
            if retry_with_backoff 2 120 -- _pull_one_image "$image" "$image_file" "$force_flag"; then
                pulled=$((pulled + 1))
            else
                still_failed+=("$image")
            fi
        done
        failed_images=("${still_failed[@]}")
    fi

    echo ""
    echo "=========================================="
    echo "  Apptainer Image Pull Summary"
    echo "=========================================="
    echo "  Total: $total | Pulled: $pulled | Skipped: $skipped | Failed: ${#failed_images[@]}"
    if [ "${#failed_images[@]}" -gt 0 ]; then
        echo "  FAILED:"
        for image in "${failed_images[@]}"; do echo "    - $image"; done
        echo "  See $(sc_install_log_dir) for details. You can re-run with --skip-existing"
        echo "  to retry only the failed images without re-pulling successful ones."
    fi
    echo "=========================================="

    # Return non-zero if everything failed; otherwise success so the
    # installer continues to the database phase.
    if [ "${#failed_images[@]}" -eq "$total" ] && [ "$total" -gt 0 ]; then
        return 1
    fi
    return 0
}

# Function to add scripts directory to PATH in a portable way.
#
# Gold-standard behaviour (cf. rustup, mise, fnm):
#   1. ALWAYS print the exact export line the user needs.
#   2. If a non-base Conda env is active, scope the change to that env via
#      activate.d/deactivate.d hooks (this is opt-in by virtue of activating
#      the env before running the installer; matches existing behaviour).
#   3. Otherwise, ASK before modifying the user's shell init file.
#      With SC_NON_INTERACTIVE=1, default to auto-add to the rc file matching
#      $SHELL (creating it if missing) so CI / Dockerfile installs still work.
add_scripts_to_path() {
    local sc_dir
    sc_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    # Helper: prepend without duplication (used for current session only)
    add_path_once() {
        if [[ ":$PATH:" != *":$1:"* ]]; then
            PATH="$1:$PATH"
            export PATH
        fi
    }

    echo ""
    echo "=========================================="
    echo "  PATH Configuration"
    echo "=========================================="
    echo "  StrainCascade scripts directory:"
    echo "    $sc_dir"
    echo ""
    echo "  To use the 'straincascade' command, this directory must be on your PATH."
    echo "  Add the following line to your shell init file (~/.bashrc, ~/.zshrc, ..):"
    echo ""
    echo "      export PATH=\"$sc_dir:\$PATH\""
    echo ""

    # Conda-active path: scope to the env (existing, well-tested behaviour).
    if [[ -n "${CONDA_DEFAULT_ENV:-}" && "${CONDA_DEFAULT_ENV:-}" != "base" ]]; then
        echo "  Detected active Conda environment: '$CONDA_DEFAULT_ENV'"
        echo "  PATH will be scoped to this environment via activate.d/deactivate.d hooks."

        local conda_env_path="${CONDA_PREFIX:?CONDA_PREFIX is not set - cannot configure PATH for conda environment}"
        local activate_script="$conda_env_path/etc/conda/activate.d/env_vars.sh"
        local deactivate_script="$conda_env_path/etc/conda/deactivate.d/env_vars.sh"

        mkdir -p "$(dirname "$activate_script")" "$(dirname "$deactivate_script")"
        {
            echo 'add_path_once() { if [[ ":$PATH:" != *":$1:"* ]]; then PATH="$1:$PATH"; fi }'
            echo "add_path_once \"$sc_dir\""
        } > "$activate_script"
        {
            echo 'remove_path() { PATH=$(echo "$PATH" | sed -e "s|:$1||g" -e "s|$1:||g" -e "s|$1||g"); }'
            echo "remove_path \"$sc_dir\""
        } > "$deactivate_script"

        add_path_once "$sc_dir"
        echo "  Done. Take effect on next 'conda activate $CONDA_DEFAULT_ENV'."
        echo "=========================================="
        return 0
    fi

    # Always update the current session immediately.
    add_path_once "$sc_dir"

    # Determine which rc file to target based on $SHELL.
    local target_rc=""
    case "${SHELL:-}" in
        */zsh)  target_rc="$HOME/.zshrc" ;;
        */bash) target_rc="$HOME/.bashrc" ;;
        */fish) target_rc="$HOME/.config/fish/config.fish" ;;
        *)
            # Fall back to the first existing rc, else ~/.bashrc
            for rc in "$HOME/.bashrc" "$HOME/.zshrc" "$HOME/.profile"; do
                if [ -f "$rc" ]; then target_rc="$rc"; break; fi
            done
            [ -z "$target_rc" ] && target_rc="$HOME/.bashrc"
            ;;
    esac

    # Check for already-present marker
    local marker="# Added by StrainCascade_installation.sh"
    if [ -f "$target_rc" ] && grep -qF "$sc_dir" "$target_rc"; then
        echo "  StrainCascade scripts directory is already in $target_rc."
        echo "=========================================="
        return 0
    fi

    local do_modify=0
    if [ "$SC_NON_INTERACTIVE" = "1" ]; then
        # In non-interactive mode, auto-add (so CI / Dockerfiles work).
        do_modify=1
        echo "  Non-interactive mode: auto-appending PATH export to $target_rc"
    else
        if ask_yes_no "  Append the export line above to $target_rc now?" "n"; then
            do_modify=1
        else
            echo "  Skipped automatic modification. Please add the line above manually."
        fi
    fi

    if [ "$do_modify" = "1" ]; then
        # Create file if it doesn't exist (e.g. fresh container / HPC account).
        mkdir -p "$(dirname "$target_rc")"
        if [ ! -f "$target_rc" ]; then
            touch "$target_rc"
            echo "  Created $target_rc"
        fi
        if [[ "$target_rc" == *config.fish ]]; then
            {
                echo ""
                echo "$marker"
                echo "if not contains $sc_dir \$PATH"
                echo "    set -gx PATH $sc_dir \$PATH"
                echo "end"
            } >> "$target_rc"
        else
            {
                echo ""
                echo "$marker"
                echo "if [[ \":\$PATH:\" != *\":$sc_dir:\"* ]]; then export PATH=\"$sc_dir:\$PATH\"; fi"
            } >> "$target_rc"
        fi
        echo "  Added PATH entry to $target_rc"
        echo "  Restart your terminal or run: source $target_rc"
    fi
    echo "=========================================="
}

# Function to check if a container image exists
check_apptainer_image_exists() {
    local image_name=$1
    if ! apptainer inspect "$image_name" &>/dev/null; then
        echo "Error: Container image $image_name not found." >&2
        return 1
    fi
    return 0
}

# Function to check if a database appears to be already installed
is_database_installed() {
    local db_name=$1
    local db_dir=$2

    case "$db_name" in
        "bakta_db")
            [[ -f "$db_dir/$db_name/db/version.json" ]]
            ;;
        "prokka_db")
            # Prokka DB is embedded in the container; listing always runs
            return 1
            ;;
        "resfinder_db")
            # resfinder also installs pointfinder_db and disinfinder_db
            [[ -d "$db_dir/resfinder_db" ]] && [[ "$(ls -A "$db_dir/resfinder_db" 2>/dev/null)" ]] && \
            [[ -d "$db_dir/pointfinder_db" ]] && [[ -d "$db_dir/disinfinder_db" ]]
            ;;
        "microbeannotator_db")
            # microbeannotator.db is created during Step 12 but is incomplete
            # if the process was killed before finishing (Parsing RefSeq KO
            # numbers). Use a separate sentinel file written only after a clean
            # exit to distinguish a complete install from a partial one.
            [[ -f "$db_dir/$db_name/.straincascade_installed" ]]
            ;;
        *)
            [[ -d "$db_dir/$db_name" ]] && [[ "$(ls -A "$db_dir/$db_name" 2>/dev/null)" ]]
            ;;
    esac
}

# Function to install individual databases
install_individual_database() {
    local db_name=$1
    local db_dir=$2
    
    echo "Installing $db_name database..."

    case "$db_name" in
        "checkm2_db")
            db_url="https://zenodo.org/api/records/5571251/files/checkm2_database.tar.gz/content"
            mkdir -p "$db_dir/$db_name"
            local archive="$db_dir/$db_name/checkm2_database.tar.gz"
            # Resumable download (-C -) with built-in curl retry; then extract.
            # Avoids the previous 'curl | tar' pipe which discarded any partial
            # download on a transient network drop.
            if retry_with_backoff 3 60 -- curl -L --fail \
                    --retry 10 --retry-all-errors --retry-delay 30 \
                    --connect-timeout 30 -C - \
                    -o "$archive" "$db_url" \
               && tar -xzf "$archive" -C "$db_dir/$db_name"; then
                rm -f "$archive"
                if [ -d "$db_dir/$db_name/CheckM2_database" ]; then
                    mv "$db_dir/$db_name/CheckM2_database/"* "$db_dir/$db_name/"
                    rmdir "$db_dir/$db_name/CheckM2_database"
                fi
                echo "$db_name database installed successfully in $db_dir/$db_name"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
                return 1
            fi
        ;;
        "gtdbtk_db") 
            mkdir -p "$db_dir/$db_name"
            primary_url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
            backup_url="https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz"
            download_success=false
            if ! command -v wget &> /dev/null; then
                echo "wget is required but not installed. Aborting." >&2
                return 1
            fi
            # Resumable (-c) wget with retries, wrapped in our backoff helper for
            # extra resilience on long-running multi-GB downloads.
            if retry_with_backoff 3 60 -- wget -c --tries=10 --waitretry=30 \
                    --read-timeout=60 --timeout=30 \
                    -P "$db_dir/$db_name" "$primary_url" \
               || retry_with_backoff 3 60 -- wget -c --tries=10 --waitretry=30 \
                    --read-timeout=60 --timeout=30 \
                    -P "$db_dir/$db_name" "$backup_url"; then
                download_success=true
            fi
            if [ "$download_success" = true ]; then
                archive_file="$db_dir/$db_name/gtdbtk_r220_data.tar.gz"
                if [ -f "$archive_file" ]; then
                    if tar -xzf "$archive_file" -C "$db_dir/$db_name"; then
                        rm "$archive_file"
                        echo "$db_name database installed successfully in $db_dir/$db_name"
                        return 0
                    else
                        echo "Error occurred during $db_name database extraction." >&2
                        return 1
                    fi
                else
                    echo "Error: Downloaded file not found at $archive_file." >&2
                    return 1
                fi
            else
                echo "Error occurred during $db_name database download from both primary and backup URLs." >&2
                return 1
            fi
        ;;
        "bakta_db")
            if ! check_apptainer_image_exists "$straincascade_genome_annotation"; then
                echo "Error: Required Apptainer image (straincascade_genome_annotation) for Bakta installation not found." >&2
                return 1
            fi

            echo "Installing Bakta database (with xz compression support)..."
            
            mkdir -p "$db_dir/$db_name"

            # Function to handle retries
            try_command() {
                local max_attempts=2
                local attempt=1
                local cmd="$1"
                
                while [ $attempt -le $max_attempts ]; do
                    echo "Attempt $attempt of $max_attempts..."
                    if eval "$cmd"; then
                        return 0
                    else
                        if [ $attempt -eq $max_attempts ]; then
                            return 1
                        fi
                        echo "Attempt $attempt failed, retrying..."
                        ((attempt++))
                        sleep 5
                    fi
                done
            }

            # Download and install Bakta database v6.0 using the container.
            # NOTE: bakta_db download (v1.11+) always attempts an internal
            # amrfinder_update at the end and exits non-zero when that sub-step
            # fails — even though the 32 GB database itself was downloaded and
            # extracted correctly. We therefore ignore the exit code and check
            # for version.json as the true success indicator. This prevents
            # re-downloading 32 GB simply because amrfinder_update was flaky.
            # We still retry up to 3 times for genuine download failures.
            local _bakta_attempt=1
            while [[ $_bakta_attempt -le 3 ]]; do
                if [[ -f "$db_dir/$db_name/db/version.json" ]]; then
                    echo "Bakta database already present (version.json found). Skipping download."
                    break
                fi
                echo "Attempt $_bakta_attempt/3: Downloading Bakta database (full type, xz compressed)..."
                apptainer exec \
                    --bind "$db_dir/$db_name":/mnt/db \
                    "${straincascade_genome_annotation}" \
                    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate bakta_env && \
                            bakta_db download --output /mnt/db --type full" || true
                if [[ -f "$db_dir/$db_name/db/version.json" ]]; then
                    break
                fi
                if [[ $_bakta_attempt -lt 3 ]]; then
                    echo "Download did not complete (version.json absent); retrying in 120s..."
                    sleep 120
                fi
                _bakta_attempt=$(( _bakta_attempt + 1 ))
            done

            if [ ! -f "$db_dir/$db_name/db/version.json" ]; then
                echo "Error: Bakta database download failed after 3 attempts (version.json not found)." >&2
                return 1
            fi
            echo "Bakta database downloaded successfully."

            if ! check_apptainer_image_exists "$straincascade_taxonomic_functional_analysis"; then
                echo "Error: Required Apptainer image (straincascade_taxonomic_functional_analysis) for populating AMRFinderPlus database not found." >&2
                return 1
            fi

            # Verify that amrfinderplus-db directory was created by bakta_db
            if [ ! -d "$db_dir/$db_name/db/amrfinderplus-db" ]; then
                echo "Warning: amrfinderplus-db directory not found. Creating it..." >&2
                mkdir -p "$db_dir/$db_name/db/amrfinderplus-db"
            fi

            # Populate AMRFinderPlus database using the SAME container that runs Bakta
            # This ensures the database version matches the AMRFinderPlus binary in bakta_env
            # Note: bakta_db creates the directory structure, but amrfinder_update populates it
            echo "Populating AMRFinderPlus database (using Bakta container for version compatibility)..."
            if ! try_command "apptainer exec \
                        --bind \"$db_dir/$db_name/db/amrfinderplus-db\":/data \
                        \"$straincascade_genome_annotation\" \
                        /usr/bin/bash -c \"
                            source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate bakta_env && \
                            amrfinder_update --force_update --database /data\""; then
                echo "Error: Failed to populate AMRFinderPlus database after 2 attempts" >&2
                return 1
            fi

            # Verify complete installation
            if [ ! -d "$db_dir/$db_name/db/amrfinderplus-db" ] || [ ! -f "$db_dir/$db_name/db/version.json" ]; then
                echo "Warning: Bakta database installation verification failed - check manually" >&2
            fi

            # Rebuild Diamond databases to ensure version compatibility
            # This prevents "Diamond failed! diamond-error-code=1" errors caused by version mismatches
            # Only attempt rebuild if source FASTA files exist (db v6.0 ships pre-built .dmnd files)
            if [ -f "$db_dir/$db_name/db/psc.faa" ] && [ -f "$db_dir/$db_name/db/pscc.faa" ]; then
                echo "Rebuilding Diamond databases for version compatibility..."
                if ! apptainer exec \
                    --bind "$db_dir/$db_name/db":/mnt/db \
                    "${straincascade_genome_annotation}" \
                    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate bakta_env && \
                            cd /mnt/db && \
                            echo 'Rebuilding psc.dmnd...' && \
                            diamond makedb --in psc.faa --db psc && \
                            echo 'Rebuilding pscc.dmnd...' && \
                            diamond makedb --in pscc.faa --db pscc && \
                            echo 'Diamond databases rebuilt successfully'"; then
                    echo "Warning: Failed to rebuild Diamond databases. Bakta may fail with Diamond errors." >&2
                    echo "You can manually rebuild by running:" >&2
                    echo "  apptainer exec --bind $db_dir/$db_name/db:/mnt/db <genome_annotation.sif> bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate bakta_env && cd /mnt/db && diamond makedb --in psc.faa --db psc && diamond makedb --in pscc.faa --db pscc'" >&2
                fi
            elif [ -f "$db_dir/$db_name/db/psc.dmnd" ] && [ -f "$db_dir/$db_name/db/pscc.dmnd" ]; then
                echo "Diamond databases already present (pre-built by bakta_db). Skipping rebuild."
            else
                echo "Warning: Neither FASTA source files nor pre-built Diamond databases found in $db_dir/$db_name/db/" >&2
                echo "Bakta may fail with Diamond errors. Check database integrity manually." >&2
            fi

            echo "Bakta database installation completed successfully"
            return 0
        ;;
        "prokka_db")
            if ! check_apptainer_image_exists "$straincascade_genome_annotation"; then
                return 1
            fi
            if output=$(apptainer exec --bind "$db_dir":/opt/conda/envs/prokka_env/db "$straincascade_genome_annotation" bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate prokka_env && prokka --listdb" 2>&1); then
                echo "Prokka databases listed successfully:"
                echo "$output"
                echo "Information: The Prokka databases are installed in the Apptainer image, NOT the local databases directory."
                return 0
            else
                echo "Error occurred during $db_name database listing." >&2
                echo "Error output: $output" >&2
                return 1
            fi
        ;;
        "microbeannotator_db")
            if ! check_apptainer_image_exists "$straincascade_genome_annotation"; then
                return 1
            fi
            # Limit threads to avoid FTP connection limits (NCBI servers typically limit to 3-5 connections)
            local max_threads=4
            local available_threads=$(( ($(getconf _NPROCESSORS_ONLN) * 4) / 5 ))
            available_threads=$(( available_threads < 1 ? 1 : available_threads ))
            local threads=$(( available_threads < max_threads ? available_threads : max_threads ))

            # MicrobeAnnotator's database builder pulls ~150 GB from NCBI FTP and
            # is the slowest step of the whole installation. Network drops over
            # 12-24h are common; wrap in retry-with-backoff so a transient
            # failure doesn't waste the entire run. The builder reuses files
            # already present on disk, so retries effectively resume.

            # PATCH: InterPro moved its FTP structure (interpro.xml.gz moved to
            # current_release/ and ftp:// is less robust than https://).
            # We extract the Python file, patch it, and bind-mount it over the read-only file.
            local patch_dir="$db_dir/$db_name/patch"
            mkdir -p "$patch_dir"
            local py_file="/opt/conda/envs/microbeannotator_env/lib/python3.7/site-packages/microbeannotator/database/conversion_database_creator.py"
            local patched_py="$patch_dir/conversion_database_creator.py"
            
            apptainer exec "$straincascade_genome_annotation" cat "$py_file" | \
                sed 's|ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz|https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro.xml.gz|g' > "$patched_py"

            # On success, write a sentinel file so --skip-existing can reliably
            # distinguish a complete install from a partially-built one
            # (microbeannotator.db is created mid-build and will exist even after
            # an OOM kill).
            if retry_with_backoff 4 600 -- apptainer exec \
               --bind "$db_dir":/data \
               --bind "$patched_py":"$py_file" \
               "$straincascade_genome_annotation" \
               microbeannotator_db_builder -d /data/microbeannotator_db -m blast -t $threads --no_aspera; then
                touch "$db_dir/$db_name/.straincascade_installed"
                echo "$db_name database installed successfully in $db_dir/microbeannotator_db"
                return 0
            else
                echo "Error occurred during $db_name database installation after multiple retries." >&2
                echo "  See log file under $(sc_install_log_dir) for details." >&2
                echo "  You can resume by re-running: ./scripts/StrainCascade_installation.sh --databases --yes --skip-existing" >&2
                return 1
            fi
        ;;
        "plasmidfinder_db")
            if ! check_apptainer_image_exists "$straincascade_assembly_qc_refinement"; then
                return 1
            fi
            if apptainer exec \
                        --bind "$db_dir":/mnt \
                        "$straincascade_assembly_qc_refinement" \
                        /usr/bin/bash -c "
                            source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate plasmidfinder_env && \
                            cd /mnt && \
                            git clone https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git"; then
                echo "$db_name database installed successfully in $db_dir/plasmidfinder_db"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
                return 1
            fi
        ;;
        "resfinder_db")
            mkdir -p "$db_dir/resfinder_db" "$db_dir/pointfinder_db" "$db_dir/disinfinder_db"
            if git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git "$db_dir/resfinder_db" && \
               git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git "$db_dir/pointfinder_db" && \
               git clone https://bitbucket.org/genomicepidemiology/disinfinder_db.git "$db_dir/disinfinder_db"; then
                echo "resfinder_db, pointfinder_db and disinfinder_db databases installed successfully in $db_dir"
                return 0
            else
                echo "Error occurred during resfinder_db, pointfinder_db, or disinfinder_db database download." >&2
                return 1
            fi
        ;;
        "amrfinderplus_db")
            if ! check_apptainer_image_exists "$straincascade_taxonomic_functional_analysis"; then
                return 1
            fi
            mkdir -p "$db_dir/$db_name"
            if apptainer exec \
                        --bind "$db_dir/$db_name":/data \
                        "$straincascade_taxonomic_functional_analysis" \
                        /usr/bin/bash -c "
                            source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate amrfinderplus_env && \
                            amrfinder_update --force_update --database /data"; then
                echo "$db_name database installed successfully in $db_dir/$db_name"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
                return 1
            fi
        ;;
        "dbcan3_db")
            if ! check_apptainer_image_exists "$straincascade_taxonomic_functional_analysis"; then
                return 1
            fi
            mkdir -p "$db_dir/$db_name"
            # dbCAN v5: use built-in database download from AWS S3
            # The original bcb.unl.edu server is permanently offline (cyberattack);
            # v5 downloads pre-built databases from AWS S3 via run_dbcan database
            if apptainer exec \
                        --bind "$db_dir/$db_name":/mnt/ \
                        "$straincascade_taxonomic_functional_analysis" \
                        /usr/bin/bash -c '
                            set -euo pipefail
                            source /opt/conda/etc/profile.d/conda.sh
                            conda activate dbcan_env
                            cd /mnt
                            test -d db && rm -rf db
                            mkdir db
                            echo "Downloading dbCAN v5 databases from AWS S3..."
                            run_dbcan database --db_dir /mnt/db --aws_s3 --cgc
                            echo "dbCAN v5 database download complete."'; then
                echo "$db_name database installed successfully in $db_dir/$db_name"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
                return 1
            fi
        ;;
        "virsorter2_db")
            if ! check_apptainer_image_exists "$straincascade_crisprcas_phage_is_elements"; then
                echo "Missing Apptainer image" >&2
                return 1
            fi

            # Run installation in container
            if apptainer exec \
                --bind "$db_dir":/mnt \
                --bind "$HOME":/home/user \
                "$straincascade_crisprcas_phage_is_elements" \
                /bin/bash -c '
                    config_dir=/mnt/virsorter2_db/config &&
                    mkdir -p "$config_dir" &&
                    unset CONDA_DEFAULT_ENV CONDA_PREFIX CONDA_PROMPT_MODIFIER CONDA_SHLVL &&
                    export PATH=/opt/conda/envs/virsorter2_env/bin:$PATH &&
                    export CONDA_PREFIX=/opt/conda/envs/virsorter2_env &&
                    export CONDA_DEFAULT_ENV=virsorter2_env &&
                    cd /mnt &&
                    # Let virsorter create config in home dir
                    virsorter setup -d virsorter2_db -j $(($(nproc) - 1)) &&
                    # Then move the config to the desired location
                    mv /home/user/.virsorter/template-config.yaml "$config_dir/" &&
                    rm -rf /home/user/.virsorter
                '; then
                echo "VirSorter2 database installed successfully"
                return 0
            else
                echo "Installation failed" >&2
                return 1
            fi
        ;;
        "genomad_db")
            if ! check_apptainer_image_exists "$straincascade_crisprcas_phage_is_elements"; then
                echo "Missing Apptainer image" >&2
                return 1
            fi

            echo "Installing geNomad database..."
            mkdir -p "$db_dir/$db_name"

            # Download geNomad database using the container
            if apptainer exec \
                --bind "$db_dir":/mnt \
                "$straincascade_crisprcas_phage_is_elements" \
                bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                         conda activate genomad_env && \
                         cd /mnt && \
                         genomad download-database ."; then
                echo "geNomad database installed successfully in $db_dir/genomad_db"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
                return 1
            fi
        ;;
        "deepfri_db")
            echo "Installing DeepFri database (CPU models)..."
            mkdir -p "$db_dir/$db_name"
            
            db_url="https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz"
            
            if ! command -v wget &> /dev/null; then
                echo "wget is required but not installed. Aborting." >&2
                return 1
            fi
            
            # Resumable (-c) download with extra retry layer.
            if retry_with_backoff 3 60 -- wget -c --tries=10 --waitretry=30 \
                    --read-timeout=60 --timeout=30 \
                    -P "$db_dir/$db_name" "$db_url"; then
                archive_file="$db_dir/$db_name/newest_trained_models.tar.gz"
                if [ -f "$archive_file" ]; then
                    echo "Extracting DeepFri models..."
                    if tar -xzf "$archive_file" -C "$db_dir/$db_name"; then
                        rm "$archive_file"
                        echo "$db_name database installed successfully in $db_dir/$db_name"
                        return 0
                    else
                        echo "Error: Failed to extract DeepFri models" >&2
                        return 1
                    fi
                else
                    echo "Error: Downloaded archive file not found" >&2
                    return 1
                fi
            else
                echo "Error occurred during $db_name database download." >&2
                return 1
            fi
        ;;
        *)
            echo "Unknown database: $db_name" >&2
            return 1
        ;;
    esac
}

# Function to install databases with user interaction
install_databases() {
    local db_location=$1

    # List of available databases
    databases=("checkm2_db" "gtdbtk_db" "bakta_db" "prokka_db" "microbeannotator_db" "plasmidfinder_db" "amrfinderplus_db" "resfinder_db" "dbcan3_db" "virsorter2_db" "genomad_db" "deepfri_db")

    local failed_dbs=""
    local succeeded_dbs=""
    local skipped_dbs=""
    local fail_count=0
    local success_count=0
    local skip_count=0

    for db in "${databases[@]}"; do
        local prompt_msg="Do you want to install $db? (y/n) "
        if is_database_installed "$db" "$db_location"; then
            prompt_msg="$db appears to be already installed. Reinstall? (y/n) "
        fi
        read -p "$prompt_msg" choice
        case "$choice" in
            y|Y )
                if install_individual_database "$db" "$db_location"; then
                    succeeded_dbs+="$db "
                    success_count=$((success_count + 1))
                else
                    failed_dbs+="$db "
                    fail_count=$((fail_count + 1))
                fi
                ;;
            * )
                echo "Skipping $db installation."
                skipped_dbs+="$db "
                skip_count=$((skip_count + 1))
                ;;
        esac
    done

    echo ""
    echo "=========================================="
    echo "  Database Installation Summary"
    echo "=========================================="
    echo "  Succeeded: $success_count | Failed: $fail_count | Skipped: $skip_count"
    if [[ -n "$succeeded_dbs" ]]; then echo "  Succeeded: $succeeded_dbs"; fi
    if [[ -n "$failed_dbs" ]]; then echo "  FAILED:    $failed_dbs"; fi
    if [[ -n "$skipped_dbs" ]]; then echo "  Skipped:   $skipped_dbs"; fi
    echo "=========================================="
}

# Function to install all databases without user interaction
install_all_databases() {
    local db_location=$1

    databases=("checkm2_db" "gtdbtk_db" "bakta_db" "prokka_db" "microbeannotator_db" "plasmidfinder_db" "amrfinderplus_db" "resfinder_db" "dbcan3_db" "virsorter2_db" "genomad_db" "deepfri_db")

    local failed_dbs=""
    local succeeded_dbs=""
    local skipped_dbs=""
    local fail_count=0
    local success_count=0
    local skip_count=0
    local total=${#databases[@]}
    local current=0

    for db in "${databases[@]}"; do
        current=$((current + 1))
        echo ""
        # Skip already-installed databases unless --force was supplied.
        if [ "$SC_FORCE" != "1" ] && [ "$SC_SKIP_EXISTING" = "1" ] \
           && is_database_installed "$db" "$db_location"; then
            echo "[$current/$total] Skipping $db (already installed; use --force to reinstall)."
            skipped_dbs+="$db "
            skip_count=$((skip_count + 1))
            continue
        fi
        echo "[$current/$total] Installing $db..."
        if install_individual_database "$db" "$db_location"; then
            succeeded_dbs+="$db "
            success_count=$((success_count + 1))
        else
            failed_dbs+="$db "
            fail_count=$((fail_count + 1))
        fi
    done

    echo ""
    echo "=========================================="
    echo "  Database Installation Summary"
    echo "=========================================="
    echo "  Total: $total | Succeeded: $success_count | Skipped: $skip_count | Failed: $fail_count"
    if [[ -n "$skipped_dbs" ]]; then echo "  Skipped: $skipped_dbs"; fi
    if [[ -n "$failed_dbs" ]]; then
        echo "  FAILED: $failed_dbs"
        echo "  Re-run with --databases --skip-existing to retry only the failed ones,"
        echo "  or inspect logs at $(sc_install_log_dir)."
    fi
    echo "=========================================="
}

# Function to get necessary Apptainer images for database installation
# Looks for images with _latest.sif suffix (from :latest docker tags)
get_apptainer_images() {
    local apptainer_img_dir="$1"
    local images=("straincascade_assembly_qc_refinement" "straincascade_genome_annotation" "straincascade_taxonomic_functional_analysis" "straincascade_crisprcas_phage_is_elements")

    for image in "${images[@]}"; do
        # First try to find _latest.sif, then fall back to any matching pattern
        local image_file="$apptainer_img_dir/${image}_latest.sif"
        if [ ! -f "$image_file" ]; then
            # Fall back to glob pattern for backwards compatibility
            local matching_files=("$apptainer_img_dir/${image}"*.sif)
            if [ ${#matching_files[@]} -eq 0 ] || [ ! -f "${matching_files[0]}" ]; then
                echo "No matching .sif file for $image database installation found in $apptainer_img_dir. Continuing with the next installation step."
                return 1
            elif [ ${#matching_files[@]} -gt 1 ]; then
                echo "Warning: Multiple matching .sif files for $image found. Using: ${matching_files[0]}"
            fi
            image_file="${matching_files[0]}"
        fi
        declare -g "${image}=${image_file}"
    done

    return 0
}

# Function to check for existing installations
check_existing_installations() {
    local existing_installations=false

    # Check for existing Apptainer images
    if [ -d "$APPTAINER_IMAGES_DIR" ] && [ "$(ls -A "$APPTAINER_IMAGES_DIR" 2>/dev/null)" ]; then
        existing_installations=true
    fi

    # Check for existing databases
    if [ -d "$DEFAULT_DB_LOCATION" ] && [ "$(ls -A "$DEFAULT_DB_LOCATION" 2>/dev/null)" ]; then
        existing_installations=true
    fi

    if $existing_installations; then
        echo "Warning: Existing Apptainer images or databases detected."
        echo "A full installation will overwrite these existing installations."
    fi
}