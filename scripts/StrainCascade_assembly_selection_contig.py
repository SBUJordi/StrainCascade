#!/usr/bin/env python3

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_assembly_selection_contig.py - Version 1.0.0

import pandas as pd
import os
import sys
import numpy as np
import logging

# Set global seeds for reproducibility
SEED = 42
np.random.seed(SEED)
pd.util.hash_pandas_object.seed = SEED
pd.options.mode.chained_assignment = None  # Suppress warnings
random_state = np.random.RandomState(SEED)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Read the TSV file into a DataFrame
tsv_file = sys.argv[1]  # Get the file path from the command-line argument
output_dir = sys.argv[2]  # Get the output directory from the command-line argument

# Check if the TSV file exists
if not os.path.exists(tsv_file):
    logger.error(f"Error: The file {tsv_file} does not exist.")
    exit(1)

try:
    df = pd.read_csv(tsv_file, sep='\t')
except Exception as e:
    logger.error(f"Error: Failed to read the file {tsv_file}. Error message: {e}")
    exit(1)

# Check if the DataFrame is empty
if df.empty:
    logger.error(f"Error: The file {tsv_file} is empty.")
    exit(1)

# Define the columns
length_columns = [
    'Total length (>= 0 bp)',
    'Total length (>= 1000 bp)',
    'Total length (>= 5000 bp)',
    'Total length (>= 10000 bp)',
    'Total length (>= 25000 bp)',
    'Total length (>= 50000 bp)'
]
contig_column = '# contigs'

# Convert all columns to numeric values
for col in length_columns + [contig_column]:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Log the data types and first few rows
logger.info("Column data types:")
logger.info(df.dtypes)
logger.info("First few rows of the DataFrame:")
logger.info(df.head())

# Check if all assemblies have identical metrics
if df.drop('Assembly', axis=1).nunique().eq(1).all():
    logger.warning("All assemblies have identical metrics.")
    chosen_assembly = df.iloc[0].to_dict()
else:
    # Initial plausibility checks
    large_assemblies = df[df['Total length (>= 0 bp)'] > 10000000]
    if not large_assemblies.empty:
        for _, assembly in large_assemblies.iterrows():
            logger.warning(f"Assembly {assembly['Assembly']} has a total length of {assembly['Total length (>= 0 bp)']} bp, which is potentially implausible for gut bacteria")
    
    # Remove implausibly large or small assemblies
    df = df[
        (df['Total length (>= 0 bp)'] <= 100000000) & 
        (df['Total length (>= 0 bp)'] >= 580000)
    ]

    if df.empty:
        logger.info("No assembly with size between 580'000 and 100'000'000 bps found")
        exit(1)

    # New selection algorithm focusing on contig number
    # Step 1: Filter out assemblies with extreme sizes
    median_length = df['Total length (>= 0 bp)'].median()
    mad_length = np.median(np.abs(df['Total length (>= 0 bp)'] - median_length))
    mad_factor = 1.96 * 1.4826

    length_filter = (
        (df['Total length (>= 0 bp)'] >= median_length - mad_factor * mad_length) &  # Minimum size threshold
        (df['Total length (>= 0 bp)'] <= median_length + mad_factor * mad_length)  # Maximum size threshold
    )
    df = df[length_filter]

    # Step 2: Sort by contig number (ascending) and select assemblies with lowest contig count
    min_contigs = df[contig_column].min()
    df = df[df[contig_column] == min_contigs]

    # Step 3: If multiple assemblies have the same contig count, look at length ratios
    if len(df) > 1:
        for col in reversed(length_columns):
            df['Ratio'] = df[col].div(df['Total length (>= 0 bp)']).fillna(0)
            max_ratio = df['Ratio'].max()
            df = df[df['Ratio'] == max_ratio]
            if len(df) == 1:
                break

    # Step 4: If still multiple assemblies, prefer circularised/best assemblies
    if len(df) > 1:
        circularised_assemblies = df[df['Assembly'].str.contains('_circularised')]
        if not circularised_assemblies.empty:
            df = circularised_assemblies
        elif not df[df['Assembly'].str.contains('best_ev2')].empty:
            df = df[df['Assembly'].str.contains('best_ev2')]
        elif not df[df['Assembly'].str.contains('best_ev1')].empty:
            df = df[df['Assembly'].str.contains('best_ev1')]

    # If there are still multiple equivalent assemblies, save them
    if len(df) > 1:
        df.to_csv(os.path.join(output_dir, 'equivalent_assemblies.tsv'), sep='\t', index=False)

    # Select the final assembly
    chosen_assembly = df.iloc[0].to_dict()

# Save the chosen assembly
chosen_assembly_df = pd.DataFrame([chosen_assembly])
chosen_assembly_df.to_csv(os.path.join(output_dir, 'chosen_assembly.tsv'), sep='\t', index=False)

# Output the assembly name
print(f"{chosen_assembly['Assembly']}")