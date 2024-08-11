#!/usr/bin/env python3

# StrainCascade_extract_organism_name_gtdbtk.py - Version 1.0.1
# Author: Sebastian Bruno Ulrich Jordi

import pandas as pd
import re
import glob
import argparse
import logging
import os
from typing import List, Optional

# Constants
TSV_EXTENSION = "*summary.tsv"
CLASSIFICATION_COLUMN = 'classification'
CLOSEST_PLACEMENT_TAXONOMY_COLUMN = 'closest_placement_taxonomy'
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
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    return parser.parse_args()

def find_tsv_file(output_dir: str) -> str:
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
    return pd.read_csv(file_path, sep='\t')

def check_columns(df: pd.DataFrame):
    required_columns = [CLASSIFICATION_COLUMN, CLOSEST_PLACEMENT_TAXONOMY_COLUMN]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

def print_dataframe_info(df: pd.DataFrame):
    logging.info(f"DataFrame shape: {df.shape}")
    logging.info("DataFrame columns:")
    for col in df.columns:
        logging.info(f"  - {col}")
    logging.info("First few rows of the DataFrame:")
    logging.info(df.head().to_string())

def get_organism_name(df: pd.DataFrame) -> str:
    unique_values_classification = df[CLASSIFICATION_COLUMN].nunique()
    unique_values_closest_placement_taxonomy = df[CLOSEST_PLACEMENT_TAXONOMY_COLUMN].nunique()

    if unique_values_classification == 0 and unique_values_closest_placement_taxonomy == 0:
        raise ValueError("No classification available")
    elif unique_values_classification > 1 and unique_values_closest_placement_taxonomy > 1:
        raise ValueError("Multiple classifications available")
    elif (unique_values_classification == 0 and unique_values_closest_placement_taxonomy > 1) or (unique_values_classification > 1 and unique_values_closest_placement_taxonomy == 0):
        raise ValueError("Unclear classification results")

    if unique_values_classification == 1:
        return df[CLASSIFICATION_COLUMN].unique()[0]
    return df[CLOSEST_PLACEMENT_TAXONOMY_COLUMN].unique()[0]

def extract_phylum(organism_name: str) -> Optional[str]:
    match = re.search(r'(p__[^;]*)', organism_name)
    return match.group(1) if match else None

def process_organism_name(organism_name: str) -> tuple[str, Optional[str]]:
    level = None
    while True:
        last_pattern = organism_name.split(';')[-1]
        if re.search(r'\d', last_pattern):
            organism_name = organism_name.rsplit(';', 1)[0]
        else:
            for prefix, taxonomy_level in TAXONOMY_LEVELS.items():
                if prefix in last_pattern:
                    organism_name = last_pattern.replace(prefix, '')
                    level = taxonomy_level
                    return organism_name, level
    return organism_name, level

def write_to_file(filename: str, content: str):
    with open(filename, "w") as file:
        file.write(str(content))

def main():
    setup_logging()
    args = parse_arguments()

    try:
        tsv_file = find_tsv_file(args.output_dir)
        df = load_dataframe(tsv_file)
        check_columns(df)
        print_dataframe_info(df)
        organism_name = get_organism_name(df)
        organism_phylum = extract_phylum(organism_name)
        organism_name, level = process_organism_name(organism_name)

        write_to_file(f"{args.output_dir}/{args.sample_name}_gtdbtk_organism_name.txt", organism_name)
        if level:
            write_to_file(f"{args.output_dir}/level_organism_name.txt", level)
        write_to_file(f"{args.output_dir}/gtdbtk_phylum_name.txt", str(organism_phylum))

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