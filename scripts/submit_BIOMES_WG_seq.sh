#!/bin/bash

# Submitter script:
# submit_BIOMES_WG_seq.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# Get the directory of this script
script_dir="${SLURM_SUBMIT_DIR}"
base_dir="$(dirname "${SLURM_SUBMIT_DIR}")"

# List all input files ending with .fastq.gz in the input_files directory
input_files=("$base_dir"/input_files/*.fastq.gz)

# Count the number of input files
num_files=${#input_files[@]}

# Calculate the number of batches
num_batches=$((num_files / 10))
remainder=$((num_files % 10))

# Submit jobs in batches of 10. If the cpus-per-task parameter is changed, adaption in the following will be necessary: lja --threads, spades.py -t, canu corThreads corConcurrency, classify_wf --cpus, prokka --cpus, bakta --threads
for ((i=0; i<num_batches; i++)); do
  start=$((i * 10))
  sbatch --array=$start-$((start + 9)) --cpus-per-task=32 $script_dir/BIOMES_WG_seq_pipeline.sh
done

# Submit a final job for the remaining files, if any. If the cpus-per-task parameter is changed, adaption in the following will be necessary: lja --threads, spades.py -t, canu corThreads corConcurrency, classify_wf --cpus, prokka --cpus, bakta --threads
if ((remainder > 0)); then
  start=$((num_batches * 10))
  sbatch --array=$start-$((start + remainder - 1)) --cpus-per-task=32 $script_dir/BIOMES_WG_seq_pipeline.sh
fi