#!/usr/bin/env python3

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_assembly_selection.py

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

# Define the length columns and contig column separately
length_columns = [
    'Total length (>= 0 bp)',
    'Total length (>= 1000 bp)',
    'Total length (>= 5000 bp)',
    'Total length (>= 10000 bp)',
    'Total length (>= 25000 bp)',
    'Total length (>= 50000 bp)'
]
contig_column = '# contigs'

# Convert all length columns and contig column to numeric values
for col in length_columns + [contig_column]:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Log the data types of the columns
logger.info("Column data types:")
logger.info(df.dtypes)

# Log the first few rows of the DataFrame
logger.info("First few rows of the DataFrame:")
logger.info(df.head())

# Check if all assemblies have identical metrics
if df.drop('Assembly', axis=1).nunique().eq(1).all():
    logger.warning("All assemblies have identical metrics.")
    chosen_assembly = df.iloc[0].to_dict()
else:
    # Warn about assemblies where 'Total length (>= 0 bp)' > 10,000,000
    large_assemblies = df[df['Total length (>= 0 bp)'] > 10000000]
    if not large_assemblies.empty:
        for _, assembly in large_assemblies.iterrows():
            logger.warning(f"Assembly {assembly['Assembly']} has a total length of {assembly['Total length (>= 0 bp)']} bp, which is potentially implausible for gut bacteria (https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a)")
    
    # Exclude lines where 'Total length (>= 0 bp)' > 100,000,000 (implausibly large by a factor of 10)
    df = df[df['Total length (>= 0 bp)'] <= 100000000]

    # If there are no lines left, print a message and exit
    if df.empty:
        logger.info("No assembly found that is smaller than 100,000,000 bp (implausible [by a factor of 10] for gut bacteria; https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a)")
        exit(1)

    # If there are more than two lines left, apply the new filtering criteria
    if len(df) > 2:
        median_assembly_length = df['Total length (>= 0 bp)'].median()
        median_contig_number = df[contig_column].median()
        mad_length = np.median(np.abs(df['Total length (>= 0 bp)'] - median_assembly_length))
        mad_contigs = np.median(np.abs(df[contig_column] - median_contig_number))

        # Correction factor for MAD
        mad_factor = 1.96 * 1.4826

        # Remove assemblies that are significantly larger AND have significantly more contigs
        df = df[~((df['Total length (>= 0 bp)'] > median_assembly_length + mad_factor * mad_length) & 
                  (df[contig_column] > median_contig_number + mad_factor * mad_contigs))]

        # Remove assemblies that are significantly smaller and below 580,000 bp
        df = df[~((df['Total length (>= 0 bp)'] < median_assembly_length - mad_factor * mad_length) & 
                  (df['Total length (>= 0 bp)'] < 580000))] # Science. 1995 Oct 20;270(5235):397-403. doi: 10.1126/science.270.5235.397

    # Initialize variables
    equivalent_assemblies = []
    chosen_assembly = None

    # Iterate through the length columns in reverse order
    for col in reversed(length_columns):
        # Calculate the ratio for each assembly, handling potential NaN values
        df['Ratio'] = df[col].div(df['Total length (>= 0 bp)']).fillna(0)

        # Find the assembly with the highest ratio
        max_ratio = df['Ratio'].max()
        selected_assemblies = df[df['Ratio'] == max_ratio]

        # If there is only one assembly, stop the process
        if len(selected_assemblies) == 1:
            chosen_assembly = selected_assemblies.iloc[0]
            break

        # If there are multiple assemblies, proceed to the next criteria
        df = selected_assemblies.copy()

    # If there is still no single winner, apply additional criteria
    if chosen_assembly is None:
        criteria_columns = [contig_column, 'Largest contig', 'Total length']
        for col in criteria_columns:
            # Find assemblies with min/max value in the current column
            best_value = df[col].min() if col == contig_column else df[col].max()
            df = df[df[col] == best_value]

            # If there is only one assembly, stop the process
            if len(df) == 1:
                chosen_assembly = df.iloc[0]
                break

        # If there is still no single winner, save the DataFrame as 'equivalent_assemblies.tsv'
        if len(df) > 1:
            df.to_csv(os.path.join(output_dir, 'equivalent_assemblies.tsv'), sep='\t', index=False)

    # Check for assemblies with the '_circularised' pattern first
    circularised_assemblies = df[df['Assembly'].str.contains('_circularised')]
    if not circularised_assemblies.empty:
        chosen_assembly = circularised_assemblies.iloc[0].to_dict()
    # Then check for assemblies with the 'best_ev2' pattern
    elif not df[df['Assembly'].str.contains('best_ev2')].empty:
        chosen_assembly = df[df['Assembly'].str.contains('best_ev2')].iloc[0].to_dict()
    # Then check for assemblies with the 'best_ev1' pattern
    elif not df[df['Assembly'].str.contains('best_ev1')].empty:
        chosen_assembly = df[df['Assembly'].str.contains('best_ev1')].iloc[0].to_dict()
    # If none of the patterns are found, choose the first assembly in the list
    else:
        chosen_assembly = df.iloc[0].to_dict()

# Ensure that a chosen assembly is always returned
if chosen_assembly is None:
    logger.error("No assembly could be selected. All assemblies might have been filtered out.")
    exit(1)

# Save the chosen assembly as 'chosen_assembly.tsv'
chosen_assembly_df = pd.DataFrame([chosen_assembly])
chosen_assembly_df.to_csv(os.path.join(output_dir, 'chosen_assembly.tsv'), sep='\t', index=False)

# Pass the value in the 'Assembly' column to the bash script
chosen_assembly_value = chosen_assembly['Assembly']
print(f"{chosen_assembly_value}")