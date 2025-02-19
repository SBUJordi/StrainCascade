#!/bin/bash

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_high_throughput_SBATCH_sequencing_files.sh
# A script for batch-submitting StrainCascade sbatch jobs for sequencing data processing

set -uo pipefail

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_directory> -p <partition> [-s <sequencing_type>] [-n <notification>] [-f <file_list>] [-c <cpus>] [-m <memory>]"
    echo "Required arguments:"
    echo "  -i <input_directory>: Directory containing input files"
    echo "  -p <partition>: SLURM partition to use for the job"
    echo "Optional arguments:"
    echo "  -s <sequencing_type>: Sequencing type (default: pacbio-hifi)"
    echo "                        Valid options: pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, nano-hq"
    echo "  -n <notification>: Email address for notifications"
    echo "  -f <file_list>: File with specific filenames to process (one filename per line)"
    echo "  -c <cpus>: Number of CPUs per task (default: 32)"
    echo "  -m <memory>: Memory per CPU in GB (default: 3)"
    echo "  -----------------------------------------------------------------"
    echo "  Example usage 1: $0 -i /path/to/input/directory -p your_partition -s pacbio-raw"
    echo "  Example usage 2: $0 -i /path/to/input/directory -p your_partition -s pacbio-hifi -n your.email@example.com -f file_list.txt -c 10 -m 5"
    exit 1
}

# Initialize variables
INPUT_DIR=""
PARTITION=""
SEQUENCING_TYPE="pacbio-hifi"
EMAIL=""
FILE_LIST=""
CPUS=32
MEMORY=3

# Parse command line arguments
while getopts "i:p:s:n:f:c:m:h" opt; do
    case ${opt} in
        i )
            INPUT_DIR=$OPTARG
            ;;
        p )
            PARTITION=$OPTARG
            ;;
        s )
            SEQUENCING_TYPE=$OPTARG
            ;;
        n )
            EMAIL=$OPTARG
            ;;
        f )
            FILE_LIST=$OPTARG
            ;;
        c )
            CPUS=$OPTARG
            ;;
        m )
            MEMORY=$OPTARG
            ;;
        h )
            usage
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_DIR" ] || [ -z "$PARTITION" ]; then
    echo "Error: Input directory (-i) and partition (-p) are required arguments."
    usage
fi

# Validate input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Validate CPU and memory values
if ! [[ "$CPUS" =~ ^[0-9]+$ ]] || [ "$CPUS" -lt 1 ]; then
    echo "Error: Invalid number of CPUs specified: $CPUS"
    exit 1
fi

if ! [[ "$MEMORY" =~ ^[0-9]+$ ]] || [ "$MEMORY" -lt 1 ]; then
    echo "Error: Invalid memory per CPU specified: $MEMORY"
    exit 1
fi

# Validate sequencing type
valid_seq_types=("pacbio-raw" "pacbio-corr" "pacbio-hifi" "nano-raw" "nano-corr" "nano-hq")
is_valid_seq_type=0
for seq_type in "${valid_seq_types[@]}"; do
    if [ "$SEQUENCING_TYPE" = "$seq_type" ]; then
        is_valid_seq_type=1
        break
    fi
done

if [ $is_valid_seq_type -eq 0 ]; then
    echo "Error: Invalid sequencing type: $SEQUENCING_TYPE"
    echo "Valid options are: ${valid_seq_types[*]}"
    exit 1
fi

# Create directories for SBATCH reports and runs
REPORT_DIR="${INPUT_DIR}/StrainCascade_sbatch_reports"
RUNS_DIR="${INPUT_DIR}/StrainCascade_runs"
mkdir -p "$REPORT_DIR" "$RUNS_DIR"

# Define run list file path
RUN_LIST="${INPUT_DIR}/straincascade_run_list.txt"

# Create array to store previously processed samples
declare -A PROCESSED_SAMPLES

# Check if run list exists and load previously processed samples
if [ -f "$RUN_LIST" ]; then
    echo "Found existing run list. Loading previously processed samples..."
    while IFS=$'\t' read -r sample_name timestamp job_id; do
        # Skip comments and empty lines
        [[ $sample_name =~ ^#.*$ || -z $sample_name ]] && continue
        PROCESSED_SAMPLES["$sample_name"]=1
        echo "Excluding previously processed sample: $sample_name"
    done < "$RUN_LIST"
else
    # Create new run list with header
    echo "# StrainCascade run list - $(date '+%Y-%m-%d %H:%M:%S')" > "$RUN_LIST"
    echo "# Format: sample_name\tsubmission_time\tslurm_job_id" >> "$RUN_LIST"
    echo "Creating new run list file: $RUN_LIST"
fi

# Function to submit a single StrainCascade job
submit_job() {
    local input_file="$1"
    local base_name=$(basename "$input_file")

    # Check if sample was previously processed
    if [ "${PROCESSED_SAMPLES[$base_name]:-0}" -eq 1 ]; then
        echo "Skipping $base_name - already processed in a previous run"
        return
    fi

    echo "Submitting job for file: $input_file"

    sbatch_command="sbatch --job-name=straincascade_${base_name} \
                          --output=${REPORT_DIR}/StrainCascade_${base_name}_%j.out \
                          --cpus-per-task=$CPUS \
                          --mem-per-cpu=${MEMORY}G \
                          --time=48:00:00 \
                          --partition=$PARTITION"

    if [ -n "$EMAIL" ]; then
        sbatch_command+=" --mail-user=$EMAIL --mail-type=end,fail"
    fi

    # Capture the job submission output
    job_submission=$($sbatch_command << EOF
#!/bin/bash
set -uo pipefail
echo "Processing file: $input_file"
straincascade -i "$input_file" -r main -t \$SLURM_CPUS_PER_TASK -s "$SEQUENCING_TYPE" -o "$RUNS_DIR"
EOF
)

    # Extract the job ID from the submission output
    if [[ $job_submission =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
        job_id="${BASH_REMATCH[1]}"
        # Add entry to run list with timestamp
        echo -e "${base_name}\t$(date '+%Y-%m-%d %H:%M:%S')\t${job_id}" >> "$RUN_LIST"
        echo "Successfully submitted job ${job_id} for ${base_name}"
    else
        echo "Warning: Failed to submit job for ${base_name}"
    fi
}

# Determine which files to process
FILES=()

if [ -n "$FILE_LIST" ]; then
    # Read the file list into the array if a file list is provided
    if [ -f "$FILE_LIST" ]; then
        while IFS= read -r line; do
            FILES+=("$line")
        done < "$FILE_LIST"
    else
        echo "Error: File list not found: $FILE_LIST"
        exit 1
    fi
else
    # Default to all files with the correct extensions if no file list is provided
    while IFS= read -r file; do
        FILES+=("$(basename "$file")")
    done < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.bam" \))
fi

# Counter for new submissions
new_submissions=0

# Loop over the determined list of files and submit jobs for each one
for file in "${FILES[@]}"; do
    full_path="$INPUT_DIR/$file"
    if [ -f "$full_path" ]; then
        if [ "${PROCESSED_SAMPLES[$(basename "$file")]:-0}" -eq 0 ]; then
            submit_job "$full_path"
            ((new_submissions++))
        fi
    else
        echo "Warning: File not found: $full_path"
    fi
done

# Final summary
echo "Processing complete:"
echo " - Found ${#FILES[@]} total input files"
echo " - ${#PROCESSED_SAMPLES[@]} previously processed samples"
echo " - $new_submissions new samples submitted"
echo "SBATCH reports will be saved in $REPORT_DIR"
echo "All StrainCascade outputs will be saved in $RUNS_DIR"
echo "Run list has been updated in $RUN_LIST"