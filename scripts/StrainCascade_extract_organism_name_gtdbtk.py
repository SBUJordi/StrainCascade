#!/usr/bin/env python3

# Copyright (c) 2024-2025, University of Bern, Department for BioMedical Research, Sebastian Bruno Ulrich JORDI
# This source code is licensed under the MIT No Attribution License (MIT-0)
# found in the LICENSE file in the root directory of this source tree.

# StrainCascade_extract_organism_name_gtdbtk.py

# Import required libraries
import pandas as pd
import re
import glob
import argparse
import logging
import os
from typing import List, Optional

# Define constants for file patterns and column names
TSV_EXTENSION = "*summary.tsv"
CLASSIFICATION_COLUMN = 'classification'
CLOSEST_PLACEMENT_TAXONOMY_COLUMN = 'closest_placement_taxonomy'
# Dictionary mapping taxonomy prefixes to their full names
TAXONOMY_LEVELS = {
    's__': 'species',
    'g__': 'genus',
    'd__': 'domain',
    'p__': 'phylum',
    'c__': 'class',
    'o__': 'order',
    'f__': 'family'
}

def setup_logging():
    # Configure logging settings for the script
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    return parser.parse_args()

def find_tsv_file(output_dir: str) -> str:
    # Find the TSV file in the specified output directory
    tsv_file_path = os.path.join(output_dir, TSV_EXTENSION)
    logging.info(f"Searching for TSV files matching: {tsv_file_path}")
    tsv_files = glob.glob(tsv_file_path)
    logging.info(f"Found TSV files: {tsv_files}")
    if not tsv_files:
        raise ValueError(f"No .tsv file found matching {tsv_file_path}")
    if len(tsv_files) > 1:
        logging.warning(f"Multiple .tsv files found. Using the first one: {tsv_files[0]}")
    return tsv_files[0]

def load_dataframe(file_path: str) -> pd.DataFrame:
    # Load the TSV file into a pandas DataFrame
    return pd.read_csv(file_path, sep='\t', na_values=['NaN', 'N/A', 'NA', 'na', 'nan', '', 'Unclassified Bacteria', 'Unclassified Archaea', 'Unclassified'])

def check_columns(df: pd.DataFrame):
    # Verify that the required columns are present in the DataFrame
    required_columns = [CLASSIFICATION_COLUMN, CLOSEST_PLACEMENT_TAXONOMY_COLUMN]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

def print_dataframe_info(df: pd.DataFrame):
    # Log information about the DataFrame for debugging purposes
    logging.info(f"DataFrame shape: {df.shape}")
    logging.info("DataFrame columns:")
    for col in df.columns:
        logging.info(f"  - {col}")
    logging.info("First few rows of the DataFrame:")
    logging.info(df.head().to_string())

def get_organism_name(df: pd.DataFrame) -> Optional[str]:
    # Extract the organism name from the DataFrame
    classification = df[CLASSIFICATION_COLUMN].iloc[0]
    closest_placement = df[CLOSEST_PLACEMENT_TAXONOMY_COLUMN].iloc[0]

    if pd.isna(classification) and pd.isna(closest_placement):
        logging.warning("Both classification and closest placement taxonomy are N/A.")
        return None

    return classification if not pd.isna(classification) else closest_placement

def extract_phylum(organism_name: Optional[str]) -> Optional[str]:
    # Extract the phylum name from the organism name
    if organism_name is None:
        return None
    match = re.search(r'(p__[^;]*)', organism_name)
    return match.group(1) if match else None

def process_organism_name(organism_name: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    # Process the organism name to extract the most specific taxonomy level
    if organism_name is None:
        return None, None
    level = None
    while organism_name:
        last_pattern = organism_name.split(';')[-1]
        if re.search(r'\d', last_pattern):
            organism_name = organism_name.rsplit(';', 1)[0]
        else:
            for prefix, taxonomy_level in TAXONOMY_LEVELS.items():
                if prefix in last_pattern:
                    organism_name = last_pattern.replace(prefix, '')
                    level = taxonomy_level
                    return organism_name, level
    return None, None

def write_to_file(filename: str, content: Optional[str]):
    # Write the given content to a file
    with open(filename, "w") as file:
        file.write(str(content) if content is not None else "")

def main():
    # Main function to orchestrate the script's workflow
    setup_logging()
    args = parse_arguments()

    try:
        # Find and load the TSV file
        tsv_file = find_tsv_file(args.output_dir)
        df = load_dataframe(tsv_file)
        
        # Verify and display DataFrame information
        check_columns(df)
        print_dataframe_info(df)
        
        # Extract and process organism information
        organism_name = get_organism_name(df)
        organism_phylum = extract_phylum(organism_name)
        organism_name, level = process_organism_name(organism_name)
        
        # Write results to files
        write_to_file(f"{args.output_dir}/{args.sample_name}_gtdbtk_organism_name.txt", organism_name)
        write_to_file(f"{args.output_dir}/level_organism_name.txt", level)
        write_to_file(f"{args.output_dir}/gtdbtk_phylum_name.txt", organism_phylum)

        logging.info("Processing completed successfully.")
    except ValueError as ve:
        logging.error(f"Value Error: {str(ve)}")
    except pd.errors.EmptyDataError:
        logging.error("The TSV file is empty.")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {str(e)}")
        logging.exception("Exception details:")

if __name__ == "__main__":
    main()