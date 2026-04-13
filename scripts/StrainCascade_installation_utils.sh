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
check_disk_space() {
    local available_kb
    available_kb=$(df -P . 2>/dev/null | awk 'NR==2 {print $4}')
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

# Function to update StrainCascade scripts
update_scripts() {
    local current_dir=$(pwd)
    local parent_dir="$(dirname "$current_dir")"
    local temp_dir="$parent_dir/script_update_temp"

    echo "Updating StrainCascade scripts..."

    # Create temporary directory
    mkdir -p "$temp_dir"
    cd "$temp_dir"

    # Clone the latest version of the public repository
    if ! git clone https://github.com/SBUJordi/StrainCascade.git; then
        echo "Failed to clone the repository. Update aborted."
        cd "$current_dir"
        rm -rf "$temp_dir"
        return 1
    fi

    # Copy new scripts to the original location
    cp -R StrainCascade/scripts/* "$current_dir"
    mkdir -p "${current_dir}/../assets"
    cp -R "StrainCascade/assets/." "${current_dir}/../assets/"
    
    # Clean up
    cd "$current_dir"
    rm -rf "$temp_dir"
    
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

    echo "StrainCascade scripts have been updated successfully."
    echo "Please restart the installation script for the changes to take effect."
    exit 0
}

# Function to pull Apptainer images
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

    for image in "${docker_images[@]}"; do
        # Convert docker image name to SIF filename (replace : with _)
        local image_basename=$(basename "$image" | tr ':' '_')
        local image_file="$APPTAINER_IMAGES_DIR/${image_basename}.sif"
        
        if [ -f "$image_file" ]; then
            read -p "Image file $image_file already exists. Overwrite? [Y/n] " choice
            case "$choice" in
                n|N ) echo "Skipping $image_file..."; continue;;
                * ) echo "Overwriting $image_file..."; force_flag="--force";;
            esac
        else
            force_flag=""
        fi
        echo "Pulling $image..."
        apptainer pull $force_flag "$image_file" "docker://$image"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to pull $image"
            exit 1
        fi
    done

    echo "All chosen images pulled successfully."
}

# Function to add scripts directory to PATH in a portable way
add_scripts_to_path() {
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    
    # Function to add path without duplication
    add_path_once() {
        if [[ ":$PATH:" != *":$1:"* ]]; then
            PATH="$1:$PATH"
        fi
    }

    # Function to remove path
    remove_path() {
        PATH=$(echo $PATH | sed -e "s|:$1||g" -e "s|$1:||g" -e "s|$1||g")
    }

    # Check if a non-base Conda environment is active
    if [[ -n "$CONDA_DEFAULT_ENV" && "$CONDA_DEFAULT_ENV" != "base" ]]; then
        echo "Conda environment '$CONDA_DEFAULT_ENV' is active."
        
        # Get the path to the active Conda environment
        conda_env_path=$(conda info --envs | grep '*' | awk '{print $NF}')
        
        # Add the script directory to the Conda environment's activation script
        activate_script="$conda_env_path/etc/conda/activate.d/env_vars.sh"
        deactivate_script="$conda_env_path/etc/conda/deactivate.d/env_vars.sh"
        
        mkdir -p "$(dirname "$activate_script")" "$(dirname "$deactivate_script")"
        
        echo "Adding scripts directory to PATH for Conda environment '$CONDA_DEFAULT_ENV'."
        echo "add_path_once() { if [[ \":\$PATH:\" != *\":\$1:\"* ]]; then PATH=\"\$1:\$PATH\"; fi }" > "$activate_script"
        echo "add_path_once \"$script_dir\"" >> "$activate_script"
        
        echo "remove_path() { PATH=\$(echo \$PATH | sed -e \"s|:\$1||g\" -e \"s|\$1:||g\" -e \"s|\$1||g\"); }" > "$deactivate_script"
        echo "remove_path \"$script_dir\"" >> "$deactivate_script"
        
        # Immediately update the current session's PATH
        add_path_once "$script_dir"
        
        echo "Scripts directory added to PATH for Conda environment '$CONDA_DEFAULT_ENV'."
        echo "The change will take effect in new terminal sessions or after reactivating the environment."
    else
        echo "No specific Conda environment active. Adding to shell config files."
        
        declare -a shell_configs=("$HOME/.bashrc" "$HOME/.zshrc" "$HOME/.config/fish/config.fish")

        for config in "${shell_configs[@]}"; do
            if [[ -f "$config" ]]; then
                if grep -q "$script_dir" "$config"; then
                    echo "Scripts directory already in PATH in $config."
                else
                    echo "Adding scripts directory to PATH in $config."
                    echo "# Added by StrainCascade_installation.sh" >> "$config"
                    echo "if [[ \":\$PATH:\" != *\":$script_dir:\"* ]]; then export PATH=\"$script_dir:\$PATH\"; fi" >> "$config"
                fi
            fi
        done

        # Immediately update the current session's PATH
        add_path_once "$script_dir"

        echo "Please restart your terminal or source your shell configuration file to update your PATH."
    fi
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
            if curl -L --connect-timeout 30 --max-time 600 "$db_url" | tar -xz -C "$db_dir/$db_name"; then
                mv "$db_dir/$db_name/CheckM2_database/"* "$db_dir/$db_name/"
                rmdir "$db_dir/$db_name/CheckM2_database"
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
            if wget --timeout=30 --tries=3 -P "$db_dir/$db_name" "$primary_url" || wget --timeout=30 --tries=3 -P "$db_dir/$db_name" "$backup_url"; then
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

            # Download and install Bakta database v6.0 using the container
            # This automatically handles xz decompression and creates the amrfinderplus-db directory structure
            if ! apptainer exec \
                --bind "$db_dir/$db_name":/mnt/db \
                "${straincascade_genome_annotation}" \
                bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                        conda activate bakta_env && \
                        echo \"Downloading Bakta database (full type, xz compressed)...\" && \
                        bakta_db download --output /mnt/db --type full"; then
                echo "Error: Failed to download and setup Bakta database" >&2
                return 1
            fi

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

            if apptainer exec --bind "$db_dir":/data "$straincascade_genome_annotation" \
               microbeannotator_db_builder -d /data/microbeannotator_db -m blast -t $threads --no_aspera; then
                echo "$db_name database installed successfully in $db_dir/microbeannotator_db"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
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
            
            # Download the trained models
            if wget --timeout=30 --tries=3 -P "$db_dir/$db_name" "$db_url"; then
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
    local fail_count=0
    local success_count=0
    local total=${#databases[@]}
    local current=0

    for db in "${databases[@]}"; do
        current=$((current + 1))
        echo ""
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
    echo "  Total: $total | Succeeded: $success_count | Failed: $fail_count"
    if [[ -n "$failed_dbs" ]]; then echo "  FAILED: $failed_dbs"; fi
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