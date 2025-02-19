#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# Function to display usage information
usage() {
    echo "Usage: $0 <partition> <email>"
    echo "  <partition>: SLURM partition to use for the job"
    echo "  <email>: Email address for notifications"
    echo "  -----------------------------------------------------------------"
    echo "  Example usage: ./scripts/SBATCH_scripts/StrainCascade_SBATCH_full_installation.sh your_partition your.email@example.com"
    exit 1
}

# Check if partition and email are provided
if [ $# -ne 2 ]; then
    usage
fi

PARTITION="$1"
EMAIL="$2"

# Create temporary SBATCH script
cat << EOF > installation_job.sh
#!/bin/bash
#SBATCH --job-name=StrainCascade_installation
#SBATCH --output=StrainCascade_installation_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --partition=$PARTITION
#SBATCH --mail-user=$EMAIL
#SBATCH --mail-type=END,FAIL

# Ensure we're in the correct directory
if [[ ! -d "./scripts" ]]; then
    echo "Error: This script must be run from the parent directory of the 'scripts' folder (= StrainCascade folder)."
    exit 1
fi

# Run the installation script and automatically answer 'y' to prompts
yes | ./scripts/StrainCascade_installation.sh
EOF

# Submit the job
sbatch installation_job.sh

# Clean up
rm installation_job.sh