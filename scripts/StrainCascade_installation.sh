#!/bin/bash

# StrainCascade_installation.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Installation script for StrainCascade

# Change to the directory containing this script
cd "$(dirname "$0")" || { echo "Failed to change to the scripts directory. Please launch the installation from within the scripts directory"; exit 1; }

# Main script commands
echo "Running StrainCascade installation..."

# Define required directories
required_dirs=("scripts" "databases" "apptainer_images")

# Define minimum disk space in GB (last measurement: 534G)
min_disk_space=550

# Default database location
default_db_location="../databases"

# Default Apptainer images directory
apptainer_images_dir="../apptainer_images"

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
        echo "Recommended disk space: ${min_disk_space}GB"
    else
        echo "Sufficient disk space available: ${available_space}GB"
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

# Function to pull Apptainer images
pull_apptainer_images() {
    mkdir -p "$apptainer_images_dir"
    # Add the list of docker images to be pulled
    docker_images=(
        "sbujordi/straincascade_genome_assembly:v1"
        "sbujordi/straincascade_assembly_qc_refinement:v1"
        "sbujordi/straincascade_genome_annotation:v1"
        "sbujordi/straincascade_taxonomic_functional_analysis:v1"
        "sbujordi/straincascade_phage_detection:v1"
        "sbujordi/python_3.12.4:v1"
        "sbujordi/r_4.4.1:v1"
        "sbujordi/cisa_1.3:v1"
    )

    for image in "${docker_images[@]}"; do
        local image_file="$apptainer_images_dir/$(basename $image).sif"
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
    script_dir=$(pwd)
    declare -a shell_configs=("$HOME/.bashrc" "$HOME/.zshrc" "$HOME/.config/fish/config.fish")

    for config in "${shell_configs[@]}"; do
        if [[ -f "$config" ]]; then
            if grep -q "$script_dir" "$config"; then
                echo "Scripts directory already in PATH in $config."
            else
                echo "Adding scripts directory to PATH in $config."
                echo "# Added by StrainCascade_installation.sh" >> "$config"
                echo "export PATH=\$PATH:$script_dir" >> "$config"
            fi
        fi
    done

    echo "Please restart your terminal or source your shell configuration file to update your PATH."
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
            if curl -L "$db_url" | tar -xz -C "$db_dir/$db_name"; then
                mv "$db_dir/$db_name/CheckM2_database/"* "$db_dir/$db_name/"
                rmdir "$db_dir/$db_name/CheckM2_database"
                echo "$db_name database installed successfully in $db_dir/$db_name"
                return 0
            else
                echo "Error occurred during $db_name database installation."
                return 1
            fi
        ;;
        "gtdbtk_db")
            mkdir -p "$db_dir/$db_name"
            if wget -P "$db_dir/$db_name" https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz; then
                if tar -xzf "$db_dir/$db_name/gtdbtk_r220_data.tar.gz" -C "$db_dir/$db_name"; then
                    rm "$db_dir/$db_name/gtdbtk_r220_data.tar.gz"
                    echo "$db_name database installed successfully in $db_dir/$db_name"
                    return 0
                else
                    echo "Error occurred during $db_name database extraction."
                    return 1
                fi
            else
                echo "Error occurred during $db_name database download."
                return 1
            fi
        ;;
        "bakta_db")
            original_dir=$(pwd)
            mkdir -p "$db_dir/$db_name"
            cd "$db_dir/$db_name"

            echo "Downloading Bakta database..."
            
            if ! curl -L -o db-versions.json "https://zenodo.org/records/10522951/files/db-versions.json?download=1"; then
                echo "Error downloading db-versions.json" >&2
                cd "$original_dir"
                return 1
            fi

            if ! echo "3c1c6746afd96e10447fa30a44f3d630 db-versions.json" | md5sum -c --; then
                echo "db-versions.json MD5 check failed." >&2
                cd "$original_dir"
                return 1
            fi

            db_version=$(jq -r '.db' db-versions.json)
            echo "Bakta database version: $db_version"

            if ! curl -L -o db.tar.gz "https://zenodo.org/records/10522951/files/db.tar.gz?download=1"; then
                echo "Error downloading main database file" >&2
                cd "$original_dir"
                return 1
            fi

            if ! echo "f8823533b789dd315025fdcc46f1a8c1 db.tar.gz" | md5sum -c --; then
                echo "Main database file MD5 check failed." >&2
                cd "$original_dir"
                return 1
            fi

            echo "Extracting database..."
            if ! tar -xzf db.tar.gz; then
                echo "Error extracting database" >&2
                cd "$original_dir"
                return 1
            fi

            rm db.tar.gz

            echo "Bakta database version $db_version downloaded and set up successfully."
            cd "$original_dir"
            return 0
        ;;
        "prokka_db")
            if output=$(apptainer exec --bind "$db_dir":/opt/conda/envs/prokka_env/db "$straincascade_genome_annotation" prokka --listdb 2>&1); then
                echo "$db_name databases listed successfully:"
                echo "$output"
                return 0
            else
                echo "Error occurred during $db_name database listing. Make sure the straincascade_genome_annotation image is available"
                echo "Error output: $output"
                return 1
            fi
        ;;
        "microbeannotator_db")
            local threads=$(( $(getconf _NPROCESSORS_ONLN) / 2 ))
            threads=$(( threads < 1 ? 1 : threads ))

            if apptainer exec --bind "$db_dir":/data "$straincascade_genome_annotation" \
               microbeannotator_db_builder -d /data/microbeannotator_db -m blast -t $threads; then
                echo "$db_name database installed successfully in $db_dir/microbeannotator_db"
                return 0
            else
                echo "Error occurred during $db_name database installation. Make sure the straincascade_genome_annotation image is available"
                return 1
            fi
        ;;
        "plasmidfinder_db")
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
                echo "Error occurred during $db_name database installation. Make sure the straincascade_assembly_qc_refinement image is available"
                return 1
            fi
        ;;
        "rgi_db")
            mkdir -p "$db_dir/$db_name"
            if wget -P "$db_dir/$db_name" https://card.mcmaster.ca/latest/data; then
                if tar -xvf "$db_dir/$db_name/data" -C "$db_dir/$db_name"; then
                    rm "$db_dir/$db_name/data"
                    echo "$db_name database installed successfully in $db_dir/$db_name"
                    return 0
                else
                    echo "Error occurred during $db_name database extraction."
                    return 1
                fi
            else
                echo "Error occurred during $db_name database download."
                return 1
            fi
        ;;
        "resfinder_db")
            mkdir -p "$db_dir/resfinder_db"
            mkdir -p "$db_dir/pointfinder_db"
            mkdir -p "$db_dir/disinfinder_db"
            # ResFinder database
            if git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git "$db_dir/resfinder_db"; then
                # PointFinder database
                if git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git "$db_dir/pointfinder_db"; then
                    # DisinFinder database
                    if git clone https://bitbucket.org/genomicepidemiology/disinfinder_db.git "$db_dir/disinfinder_db"; then
                        echo "resfinder_db, pointfinder_db and disinfinder_db databases installed successfully in $db_dir"
                        return 0
                    else
                        echo "Error occurred during disinfinder_db database download."
                        return 1
                    fi
                else
                    echo "Error occurred during pointfinder_db database download."
                    return 1
                fi
            else
                echo "Error occurred during resfinder_db database download."
                return 1
            fi
        ;;
        "dbcan3_db")
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
                echo "Error occurred during $db_name database installation. Make sure the straincascade_taxonomic_functional_analysis image is available"
                return 1
            fi
        ;;
        "virsorter2_db")
            max_cpus=$(($(nproc) - 1))
            if apptainer exec \
                        --bind "$db_dir":/mnt/ \
                        "$straincascade_phage_detection" \
                        /usr/bin/bash -c "
                            source /opt/conda/etc/profile.d/conda.sh && \
                            conda activate virsorter2_env && \
                            cd /mnt && \
                            virsorter setup -d virsorter2_db -j $max_cpus"; then
                echo "$db_name database installed successfully in $db_dir/virsorter2_db"
                return 0
            else
                echo "Error occurred during $db_name database installation. Make sure the straincascade_phage_detection image is available"
                return 1
            fi
        ;;
        *)
            echo "Unknown database: $db_name"
            return 1
        ;;
    esac
}

# Function to install databases
install_databases() {
    local db_location=$1
    
    # List of available databases
    databases=("checkm2_db" "gtdbtk_db" "bakta_db" "prokka_db" "microbeannotator_db" "plasmidfinder_db" "rgi_db" "resfinder_db" "dbcan3_db" "virsorter2_db")
    
    for db in "${databases[@]}"; do
        read -p "Do you want to install $db? (y/n) " choice
        case "$choice" in
            y|Y ) install_individual_database "$db" "$db_location";;
            n|N ) echo "Skipping $db installation.";;
            * ) echo "Invalid choice. Skipping $db installation.";;
        esac
    done
}

# Function to get necessary Apptainer images for database installation
get_apptainer_images() {
    local apptainer_images_dir="$1"
    local images=("straincascade_assembly_qc_refinement" "straincascade_genome_annotation" "straincascade_taxonomic_functional_analysis" "straincascade_phage_detection")

    for image in "${images[@]}"; do
        local image_file=("$apptainer_images_dir/${image}"*.sif)
        if [ "${#image_file[@]}" -eq 0 ]; then
            echo "No matching .sif file for $image database installation found in $apptainer_images_dir. Continuing with the next installation step."
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

# Main script execution
check_main_directory
check_disk_space
check_apptainer_installation

# Ask for full installation or individual components
read -p "Do you want to do a full installation of StrainCascade? (y/n) " full_install_choice
case "$full_install_choice" in
    y|Y)
        pull_apptainer_images
        add_scripts_to_path

        if ! get_apptainer_images "$apptainer_images_dir"; then
            echo "Some Apptainer images were not found. Continuing without database installation."
        else
            # Database installation
            read -p "Do you want to install databases at the default location? (y/n) " install_db_choice
            case "$install_db_choice" in
                y|Y) 
                    db_location="$default_db_location"
                    ;;
                n|N) 
                    read -p "Enter the alternative directory for database installation: " alt_db_location
                    db_location="$alt_db_location"
                    ;;
                *) 
                    echo "Invalid choice. Skipping database installation."
                    db_location=""
                    ;;
            esac

            if [ -n "$db_location" ]; then
                mkdir -p "$db_location"
                install_databases "$db_location"
            fi
        fi
        ;;
    n|N)
        # Ask if user wants to install Apptainer images
        read -p "Do you want to install any Apptainer images? (y/n) " install_apptainer_choice
        case "$install_apptainer_choice" in
            y|Y)
                pull_apptainer_images
                add_scripts_to_path
                ;;
            n|N)
                # Ask if user wants to install databases
                read -p "Do you want to install any databases? (y/n) " install_db_choice
                case "$install_db_choice" in
                    y|Y)
                        if ! get_apptainer_images "$apptainer_images_dir"; then
                            echo "Some Apptainer images were not found. Continuing without database installation."
                        else
                            # Database installation
                            read -p "Do you want to install databases at the default location? (y/n) " install_db_location_choice
                            case "$install_db_location_choice" in
                                y|Y)
                                    db_location="$default_db_location"
                                    ;;
                                n|N)
                                    read -p "Enter the alternative directory for database installation: " alt_db_location
                                    db_location="$alt_db_location"
                                    ;;
                                *)
                                    echo "Invalid choice. Skipping database installation."
                                    db_location=""
                                    ;;
                            esac

                            if [ -n "$db_location" ]; then
                                mkdir -p "$db_location"
                                install_databases "$db_location"
                            fi
                        fi
                        ;;
                    n|N)
                        echo "Skipping database installation."
                        ;;
                    *)
                        echo "Invalid choice. Skipping database installation."
                        ;;
                esac
                ;;
            *)
                echo "Invalid choice. Skipping Apptainer image and database installation."
                ;;
        esac
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo "Installation process completed"