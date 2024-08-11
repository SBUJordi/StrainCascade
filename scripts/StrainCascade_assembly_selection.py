# StrainCascade_assembly_selection.py - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

import pandas as pd
import random
import os
import sys
import numpy as np
import logging

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

# Define the length columns
length_columns = [
    'Total length (>= 0 bp)',
    'Total length (>= 1000 bp)',
    'Total length (>= 5000 bp)',
    'Total length (>= 10000 bp)',
    'Total length (>= 25000 bp)',
    'Total length (>= 50000 bp)'
]

# Convert all length columns to numeric values
for col in length_columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Log the data types of the columns
logger.info("Column data types:")
logger.info(df.dtypes)

# Log the first few rows of the DataFrame
logger.info("First few rows of the DataFrame:")
logger.info(df.head())

# Exclude lines where 'Total length (>= 0 bp)' > 10000000 (implausibly large for gut bacteria (https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a))
df = df[df['Total length (>= 0 bp)'] <= 10000000]

# If there are no lines left, print a message and exit
if df.empty:
    logger.info("No assembly below 10'000'000 bp (implausible for gut bacteria; https://doi.org/10.1038/s41467-023-37396-x, Fig. 1a)")
    exit(1)

# If there are more than two lines left, calculate the median and exclude lines with deviation >= 1.96 standard deviations
if len(df) > 2:
    median_assembly_length = df['Total length (>= 0 bp)'].median()
    df['deviation_from_median_assembly_length'] = df['Total length (>= 0 bp)'] - median_assembly_length
    std_dev = df['deviation_from_median_assembly_length'].std()
    df = df[np.abs(df['deviation_from_median_assembly_length']) < 1.96 * std_dev]

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
    criteria_columns = ['# contigs', 'Largest contig', 'Total length']
    for col in criteria_columns:
        # Find assemblies with min/max value in the current column
        best_value = df[col].min() if col == '# contigs' else df[col].max()
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

# Save the chosen assembly as 'chosen_assembly.tsv'
chosen_assembly_df = pd.DataFrame([chosen_assembly])
chosen_assembly_df.to_csv(os.path.join(output_dir, 'chosen_assembly.tsv'), sep='\t', index=False)

# Pass the value in the 'Assembly' column to the bash script
chosen_assembly_value = chosen_assembly['Assembly']
print(f"{chosen_assembly_value}")