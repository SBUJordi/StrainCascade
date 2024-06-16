
### apptainer_images/info.md

# Apptainer Images Directory

## Overview
This directory stores the Apptainer (formerly Singularity) images pulled from Docker. These images contain the necessary software and dependencies to run the pipeline.

## Usage
- The `install.sh` script will pull the required Docker images and convert them to Apptainer images, placing them in this directory.
- The pipeline scripts will reference these images to ensure the correct software environment is used.

## Apptainer image list
- 

## Notes
- This directory is initially empty and will be populated by the installation script.
- Ensure that Apptainer is installed and configured on your system before running the installation script.
- If you need to manually add or update images, place the `.sif` files in this directory and update any relevant configuration files or scripts.
