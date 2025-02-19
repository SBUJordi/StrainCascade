#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_installation_utils.sh
# Description: Utility functions for StrainCascade installation

# Define required directories
required_dirs=("scripts" "databases" "apptainer_images")

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
}

# Function to check disk space
check_disk_space() {
    available_space=$(df -BG --output=avail . | tail -1 | sed 's/G//')
    if [ "$available_space" -lt "$min_disk_space" ]; then
        echo "Warning: Insufficient disk space for full StrainCascade installation."
        echo "Available disk space: ${available_space}GB"
        echo "Recommended minimal disk space: ${min_disk_space}GB"
    else
        echo "Sufficient disk space available: ${available_space}GB"
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

# Function to update StrainCascade scripts
update_scripts() {
    local current_dir=$(pwd)
    local parent_dir="$(dirname "$current_dir")"
    local temp_dir="$parent_dir/script_update_temp"

    echo "Updating StrainCascade scripts..."

    # Create temporary directory
    mkdir -p "$temp_dir"
    cd "$temp_dir"

    # Clone the latest version of the repository (pre-release with read-only access token)
    if ! git clone https://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@github.com/SBUJordi/StrainCascade.git; then
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
    # Add the list of docker images to be pulled
    docker_images=(
        "sbujordi/straincascade_genome_assembly:v3"
        "sbujordi/straincascade_lja_genome_assembly:v1"
        "sbujordi/straincascade_assembly_qc_refinement:v2"
        "sbujordi/straincascade_genome_annotation:v2"
        "sbujordi/straincascade_taxonomic_functional_analysis:v2"
        "sbujordi/straincascade_crisprcas_phage_is_elements:v1"
        "sbujordi/python_3.12.4:v1"
        "sbujordi/r_4.4.1:v2"
        "sbujordi/straincascade_document_processing:v1"
    )

    for image in "${docker_images[@]}"; do
        local image_file="$APPTAINER_IMAGES_DIR/$(basename $image).sif"
        if [ -f "$image_file" ]; then
            read -p "Image file $image_file already exists. Do you want to overwrite it? (y/n) " choice
            case "$choice" in
                y|Y ) echo "Overwriting $image_file..."; force_flag="--force";;
                n|N ) echo "Skipping $image_file..."; continue;;
                * ) echo "Invalid choice. Skipping $image_file..."; continue;;
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

            echo "Installing Bakta database..."
            
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

            # Download and install Bakta database using the container
            if ! apptainer exec \
                --bind "$db_dir/$db_name":/mnt/db \
                "${straincascade_genome_annotation}" \
                bash -c "source /opt/conda/etc/profile.d/conda.sh && \
                        conda activate bakta_env && \
                        echo \"Downloading full Bakta database...\" && \
                        bakta_db download --output /mnt/db --type full"; then
                echo "Error: Failed to download and setup Bakta databases" >&2
                return 1
            fi

            if ! check_apptainer_image_exists "$straincascade_taxonomic_functional_analysis"; then
                echo "Error: Required Apptainer image (straincascade_taxonomic_functional_analysis) for installing/updating AMRPlusFinder database for Bakta not found." >&2
                return 1
            fi

            # Create the directory for AMRPlusFinder database if not present yet
            mkdir -p "$db_dir/$db_name/db/amrfinderplus-db"

            # Install AMRPlusFinder database using the container with retry
            if ! try_command "apptainer exec \
                        --bind \"$db_dir/$db_name/db/amrfinderplus-db\":/data \
                        \"$straincascade_taxonomic_functional_analysis\" \
                        /usr/bin/bash -c \"
                            source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate amrfinderplus_env && \
                            amrfinder_update --force_update --database /data\""; then
                echo "Error: Failed to install AMRPlusFinder database after 2 attempts" >&2
                return 1
            fi

            # Verify installation
            if [ ! -d "$db_dir/$db_name/db/amrfinderplus-db" ] || [ ! -f "$db_dir/$db_name/db/version.json" ]; then
                echo "Warning: Bakta database installation verification failed - check manually" >&2
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
            local threads=$(( ($(getconf _NPROCESSORS_ONLN) * 4) / 5 ))  # Multiply by 4 and divide by 5 to approximate division by 1.25
            threads=$(( threads < 1 ? 1 : threads ))

            if apptainer exec --bind "$db_dir":/data "$straincascade_genome_annotation" \
               microbeannotator_db_builder -d /data/microbeannotator_db -m blast -t $threads; then
                echo "$db_name database installed successfully in $db_dir/microbeannotator_db"
                return 0
            else
                echo "Error occurred during $db_name database installation." >&2
                echo "If this error is related to Aspera (Connect) you can ignore it; StrainCascade does not use Aspera Connect"
                return 0
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
            max_cpus=$(($(nproc) - 1))
            mkdir -p "$db_dir/$db_name"
            if apptainer exec \
                        --bind "$db_dir/$db_name":/mnt/ \
                        "$straincascade_taxonomic_functional_analysis" \
                        /usr/bin/bash -c "
                            source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate dbcan_env && \
                            cd /mnt && \
                            dbcan_build --cpus $max_cpus --db-dir db --clean"; then
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
    databases=("checkm2_db" "gtdbtk_db" "bakta_db" "prokka_db" "microbeannotator_db" "plasmidfinder_db" "amrfinderplus_db" "resfinder_db" "dbcan3_db" "virsorter2_db")
    
    for db in "${databases[@]}"; do
        read -p "Do you want to install $db? (y/n) " choice
        case "$choice" in
            y|Y ) install_individual_database "$db" "$db_location";;
            n|N ) echo "Skipping $db installation.";;
            * ) echo "Invalid choice. Skipping $db installation.";;
        esac
    done
}

# Function to install all databases without user interaction
install_all_databases() {
    local db_location=$1
    
    databases=("checkm2_db" "gtdbtk_db" "bakta_db" "prokka_db" "microbeannotator_db" "plasmidfinder_db" "amrfinderplus_db" "resfinder_db" "dbcan3_db" "virsorter2_db")
    
    for db in "${databases[@]}"; do
        echo "Installing $db..."
        install_individual_database "$db" "$db_location"
    done
}

# Function to get necessary Apptainer images for database installation
get_apptainer_images() {
    local APPTAINER_IMAGES_DIR="$1"
    local images=("straincascade_assembly_qc_refinement" "straincascade_genome_annotation" "straincascade_taxonomic_functional_analysis" "straincascade_crisprcas_phage_is_elements")

    for image in "${images[@]}"; do
        local image_file=("$APPTAINER_IMAGES_DIR/${image}"*.sif)
        if [ "${#image_file[@]}" -eq 0 ]; then
            echo "No matching .sif file for $image database installation found in $APPTAINER_IMAGES_DIR. Continuing with the next installation step."
            return 1
        elif [ "${#image_file[@]}" -gt 1 ]; then
            echo "Warning: Multiple matching .sif files for $image database installation found. Using the first match: ${image_file[0]}"
            declare -g "${image}=${image_file[0]}"
        else
            declare -g "${image}=${image_file[0]}"
        fi
    done

    return 0
}

# Function to check for existing installations
check_existing_installations() {
    local existing_installations=false

    # Check for existing Apptainer images
    if [ "$(ls -A "$APPTAINER_IMAGES_DIR")" ]; then
        existing_installations=true
    fi

    # Check for existing databases
    if [ -d "$DEFAULT_DB_LOCATION" ] && [ "$(ls -A "$DEFAULT_DB_LOCATION")" ]; then
        existing_installations=true
    fi

    if $existing_installations; then
        echo "Warning: Existing Apptainer images or databases detected."
        echo "A full installation will overwrite these existing installations."
    fi
}