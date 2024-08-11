import pandas as pd
import numpy as np
import sys
import re

# Get the paths to the TSV files from the command line arguments
quast_assembly_evaluation_tsv = sys.argv[1]
taxonomy_gtdbtk_tsv = sys.argv[2]
annotation_prokka_tsv = sys.argv[3]
annotation_bakta_tsv = sys.argv[4]
annotation_bakta_hypotheticals_tsv = sys.argv[5]

# Get the parent directory from the command line arguments
parent_directory = sys.argv[6]
sample_name = sys.argv[7]

# Read the TSV files into pandas DataFrames
quast_assembly_evaluation_df = pd.read_csv(quast_assembly_evaluation_tsv, sep='\t', na_values=[''])
taxonomy_gtdbtk_df = pd.read_csv(taxonomy_gtdbtk_tsv, sep='\t', na_values=[''])
annotation_prokka_df = pd.read_csv(annotation_prokka_tsv, sep='\t', dtype=str, na_values=[''])

# Then use the header_row variable to skip the appropriate number of rows
annotation_bakta_df = pd.read_csv(annotation_bakta_tsv, sep='\t', skiprows=5, header=0, dtype=str, na_values=[''])
annotation_bakta_hypotheticals_df = pd.read_csv(annotation_bakta_hypotheticals_tsv, sep='\t', skiprows=2, header=0, dtype=str, na_values=[''])

# Process the QUAST assembly evaluation TSV
# Define the new column names
new_column_names = {
    "Assembly": "assembly",
    "# contigs": "number_contigs",
    "# contigs (>= 1000 bp)": "number_contigs_>=1000_bp",
    "# contigs (>= 5000 bp)": "number_contigs_>=5000_bp",
    "# contigs (>= 10000 bp)": "number_contigs_>=10000_bp",
    "# contigs (>= 25000 bp)": "number_contigs_>=25000_bp",
    "# contigs (>= 50000 bp)": "number_contigs_>=50000_bp",
    "Total length (>= 0 bp)": "total_length_>=0_bp",
    "Largest contig": "largest_contig"
}

# Rename the columns
quast_assembly_evaluation_df = quast_assembly_evaluation_df.rename(columns=new_column_names)

# Filter the DataFrame to keep only the renamed columns
quast_assembly_evaluation_df = quast_assembly_evaluation_df[list(new_column_names.values())]

# Filter the DataFrame to keep only the rows where the assembly column contains "circular"
quast_assembly_evaluation_df = quast_assembly_evaluation_df[quast_assembly_evaluation_df['assembly'].str.contains("circular")]


# Process the GTDB-Tk taxonomic classification TSV
# Define the new column names
new_column_names = {
    "classification": "gtdbtk_classification",
    "closest_placement_taxonomy": "gtdbtk_closest_placement_taxonomy",
    "closest_placement_reference": "gtdbtk_closest_placement_reference",
    "closest_placement_ani": "gtdbtk_closest_placement_ani",
    "note": "gtdbtk_note",
    "warnings": "gtdbtk_warnings"
}

# Rename the columns
taxonomy_gtdbtk_df = taxonomy_gtdbtk_df.rename(columns=new_column_names)

# Filter the DataFrame to keep only the renamed columns
taxonomy_gtdbtk_df = taxonomy_gtdbtk_df[list(new_column_names.values())]


# Process the Prokka annotation TSV
# Define the new column names
new_column_names = {
    "COG": "cog",
    "gene": "gene",
    "product": "product"
}

# Rename the columns
annotation_prokka_df = annotation_prokka_df.rename(columns=new_column_names)

# Filter the DataFrame to keep only the renamed columns
annotation_prokka_df = annotation_prokka_df[list(new_column_names.values())]

# Create a new DataFrame with the rows where the product column contains "hypothetical protein"
annotation_prokka_df['product'] = annotation_prokka_df['product'].astype(str)
annotation_prokka_hypotheticals_df = annotation_prokka_df[annotation_prokka_df['product'] == "hypothetical protein"]

# Remove these rows from the original DataFrame
annotation_prokka_df = annotation_prokka_df[annotation_prokka_df['product'] != "hypothetical protein"]

# Drop the product column from both DataFrames
annotation_prokka_df = annotation_prokka_df.drop(columns=['product'])
annotation_prokka_hypotheticals_df = annotation_prokka_hypotheticals_df.drop(columns=['product'])


# Process the Bakta annotation TSVs
# Keep only the required columns and rename them
# Define the new column names
new_column_names = {
    "Gene": "gene",
    "DbXrefs": "cog"
}

# Rename the columns
annotation_bakta_df = annotation_bakta_df.rename(columns=new_column_names)

# Convert the cog column to a string
annotation_bakta_df['cog'] = annotation_bakta_df['cog'].astype(str)

# Extract the COG code from the DbXrefs column
annotation_bakta_df['cog'] = annotation_bakta_df['cog'].apply(
    lambda x: re.search('COG:(.*?),', x).group(1) if re.search('COG:(.*?),', x) is not None else np.nan
)

# For the hypotheticals DataFrame, create the required columns and drop the rest
annotation_bakta_hypotheticals_df['gene'] = 'hypothetical protein'
annotation_bakta_hypotheticals_df['cog'] = np.nan
annotation_bakta_hypotheticals_df = annotation_bakta_hypotheticals_df[['gene', 'cog']]


# Summarise the annotation data
# Replace empty strings with NaN
annotation_prokka_df.replace('', np.nan, inplace=True)
annotation_bakta_df.replace('', np.nan, inplace=True)
annotation_prokka_hypotheticals_df.replace('', np.nan, inplace=True)
annotation_bakta_hypotheticals_df.replace('', np.nan, inplace=True)

# Calculate the values for each column
number_genes_names_prokka = annotation_prokka_df['gene'].count()
number_genes_prokka_with_COGref = annotation_prokka_df['cog'].count()
number_genes_names_Bakta = annotation_bakta_df['gene'].count()
number_genes_Bakta_with_COGref = annotation_bakta_df['cog'].count()
overlap_genes_name_prokka_Bakta = len(set(annotation_prokka_df['gene'].dropna()).intersection(annotation_bakta_df['gene'].dropna()))
overlap_genes_COGref_prokka_Bakta = len(set(annotation_prokka_df['cog'].dropna()).intersection(annotation_bakta_df['cog'].dropna()))
number_of_hypotheticals_prokka = annotation_prokka_hypotheticals_df['gene'].count()
number_of_hypotheticals_Bakta = annotation_bakta_hypotheticals_df['gene'].count()
percentage_of_hypotheticals_prokka = number_of_hypotheticals_prokka / (number_genes_names_prokka + number_of_hypotheticals_prokka) if (number_genes_names_prokka + number_of_hypotheticals_prokka) != 0 else np.nan
percentage_of_hypotheticals_Bakta = number_of_hypotheticals_Bakta / (number_genes_names_Bakta + number_of_hypotheticals_Bakta) if (number_genes_names_Bakta + number_of_hypotheticals_Bakta) != 0 else np.nan

# Create the gene_annotation_summary DataFrame
gene_annotation_summary = pd.DataFrame({
    'number_genes_names_prokka': [number_genes_names_prokka],
    'number_genes_prokka_with_COGref': [number_genes_prokka_with_COGref],
    'number_genes_names_Bakta': [number_genes_names_Bakta],
    'number_genes_Bakta_with_COGref': [number_genes_Bakta_with_COGref],
    'overlap_genes_name_prokka_Bakta': [overlap_genes_name_prokka_Bakta],
    'overlap_genes_COGref_prokka_Bakta': [overlap_genes_COGref_prokka_Bakta],
    'number_of_hypotheticals_prokka': [number_of_hypotheticals_prokka],
    'number_of_hypotheticals_Bakta': [number_of_hypotheticals_Bakta],
    'percentage_of_hypotheticals_prokka': [percentage_of_hypotheticals_prokka],
    'percentage_of_hypotheticals_Bakta': [percentage_of_hypotheticals_Bakta]
})

# Combine the dataframes
quast_assembly_evaluation_df = quast_assembly_evaluation_df.reset_index(drop=True)
taxonomy_gtdbtk_df = taxonomy_gtdbtk_df.reset_index(drop=True)
gene_annotation_summary = gene_annotation_summary.reset_index(drop=True)

# result_summary = pd.concat([quast_assembly_evaluation_df, taxonomy_gtdbtk_df, gene_annotation_summary], axis=1, join='outer') # This is the original line but the annotation summaries are not correct yet
result_summary = pd.concat([quast_assembly_evaluation_df, taxonomy_gtdbtk_df], axis=1, join='outer')

# Write the resulting dataframe to a TSV file
result_summary.to_csv(f'{parent_directory}/{sample_name}_BIOMES_WGseq_short_summary.tsv', sep='\t', index=False)