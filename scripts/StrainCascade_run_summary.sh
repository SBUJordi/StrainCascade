#!/bin/bash

# StrainCascade_run_summary.sh - Version 1.0.1
# Author: Sebastian Bruno Ulrich Jordi

# Assign the command line arguments to named variables
utils_file=$3
apptainer_images_dir=$4
input_file=$5
sample_name=$6
sequencing_type=$7
threads=$8
databases_dir=$9
main_results_dir_abs=${10}
selected_modules=${11}
version=${12}

# Convert the selected_modules string to an array
IFS=' ' read -r -a selected_modules <<< "$selected_modules"

source "$utils_file"

# Function to check if a module was selected
was_module_selected() {
    local sc_info="$1"
    for module in "${selected_modules[@]}"; do
        if [[ "$sc_info" == *"$module"* ]]; then
            return 0  # True, module was selected
        fi
    done
    return 1  # False, module was not selected
}

# Create an array to store selected SC info
selected_sc_info=()

# Find the necessary .sif files
straincascade_python=$(find_apptainer_sif_file "$apptainer_images_dir" 'python_3.12.4*.sif')
straincascade_r=$(find_apptainer_sif_file "$apptainer_images_dir" 'r_4.4.1*.sif')
straincascade_genome_assembly=$(find_apptainer_sif_file "$apptainer_images_dir" 'straincascade_genome_assembly*.sif')
straincascade_assembly_qc_refinement=$(find_apptainer_sif_file "$apptainer_images_dir" 'straincascade_assembly_qc_refinement*.sif')
straincascade_genome_annotation=$(find_apptainer_sif_file "$apptainer_images_dir" 'straincascade_genome_annotation*.sif')
straincascade_taxonomic_functional_analysis=$(find_apptainer_sif_file "$apptainer_images_dir" 'straincascade_taxonomic_functional_analysis*.sif')
straincascade_phage_detection=$(find_apptainer_sif_file "$apptainer_images_dir" 'straincascade_phage_detection*.sif')

# Define the version and get the current date-time
version="1.0.0"
current_date_time=$(date "+%Y-%m-%d %H:%M:%S")
os_name=$(uname -s)
os_version=$(cat /etc/os-release | grep "^PRETTY_NAME" | cut -d '=' -f 2 | tr -d '"')
kernel_version=$(uname -r)
architecture=$(uname -m)

# Define software versions
straincascade_python_version=$(apptainer exec "$straincascade_python" python --version 2>&1)

straincascade_r_version=$(apptainer exec "$straincascade_r" R --version 2>&1 | grep "^R version")

canu_version=$(apptainer exec "$straincascade_genome_assembly" bash -c "source activate genome_assembly_env && canu --version 2>&1 | grep '^canu' | awk '{print \$2}'")

lja_version="LJA 0.2"

spades_version=$(apptainer exec "$straincascade_genome_assembly" bash -c "source activate spades_env && spades.py --version 2>&1 | tail -n 1")

flye_version=$(apptainer exec "$straincascade_genome_assembly" bash -c "source activate genome_assembly_env && flye --version 2>&1")

quast_version=$(apptainer exec "$straincascade_assembly_qc_refinement" bash -c "source activate quast_env && quast.py --version 2>&1 | grep 'QUAST' | awk '{print $1, $2}'")

mac2_version="MAC2.0 v2.1 (with MUMmer 4.0.0rc1)"

circlator_version=$(apptainer exec "$straincascade_assembly_qc_refinement" bash -c "source activate circlator_env && circlator version 2>&1")

ngmlr_version=$(apptainer exec "$straincascade_assembly_qc_refinement" bash -c "source activate bbmap_ngmlr_env && ngmlr --version 2>&1 | grep 'ngmlr' | awk '{print $1, $2}'")

bbmap_version=$(apptainer exec "$straincascade_assembly_qc_refinement" bash -c "source activate bbmap_ngmlr_env && pileup.sh --version 2>&1 | grep 'BBMap' | awk '{print $1, "v"$3}'")

checkm2_version=$(apptainer exec "$straincascade_assembly_qc_refinement" bash -c "source activate checkm2_env && checkm2 --version 2>&1 | grep 'checkm2' | awk '{print $2}'")
checkm2_db_version=$(grep '"version"' "$databases_dir/checkm2_db/CONTENTS.json" | awk -F '"' '{print $10}')

gtdbtk_version=$(apptainer exec "$straincascade_taxonomic_functional_analysis" bash -c "source activate gtdbtk_env && gtdbtk --version 2>&1 | grep 'gtdbtk' | awk '{print "GTDB-Tk", "v"$3}'")
gtdbtk_db_version=$(ls "$databases_dir/gtdbtk_db/")

arrow_version="arrow (GenomicConsensus v2.3.3)"
medaka_version="medaka v1.11.3"

bakta_version=$(apptainer exec "$straincascade_genome_annotation" bash -c "source activate bakta_env && bakta --version 2>&1 | grep 'bakta' | awk '{print $1, "v"$2}'")
bakta_db_version_json_file="$databases_dir/bakta_db/db/version.json"
# Extract overall database information
bakta_db_date=$(grep -oP '"date": "\K[^"]+' "$bakta_db_version_json_file")
bakta_db_major=$(grep -oP '"major": \K\d+' "$bakta_db_version_json_file" | head -n 1)
bakta_db_minor=$(grep -oP '"minor": \K\d+' "$bakta_db_version_json_file" | head -n 1)
bakta_db_type=$(grep -oP '"type": "\K[^"]+' "$bakta_db_version_json_file")
bakta_db_doi=$(grep -oP '"doi": "\K[^"]+' "$bakta_db_version_json_file")
# Initialize the final string with overall database information
bakta_db_version="Overall Bakta database:\n"
bakta_db_version+="Date: $bakta_db_date"
bakta_db_version+="Major: $bakta_db_major;"
bakta_db_version+="Minor: $bakta_db_minor;"
bakta_db_version+="Type: $bakta_db_type;"
bakta_db_version+="DOI: $bakta_db_doi;"
bakta_db_version+="Bakta database dependencies:\n"
# Extract and append dependencies information
grep -oP '"name": "\K[^"]+' "$bakta_db_version_json_file" | while read -r name; do
    release=$(grep -A 1 "\"name\": \"$name\"" "$bakta_db_version_json_file" | grep -oP '"release": "\K[^"]+')
    bakta_db_version+="$name: $release\n"
done

prokka_version=$(apptainer exec "$straincascade_genome_annotation" bash -c "source activate prokka_env && prokka --version 2>&1 | grep 'prokka' | awk '{print $1, $2}'")

microbeannotator_version=$(apptainer exec "$straincascade_genome_annotation" bash -c "source activate microbeannotator_env && microbeannotator --version 2>&1 | grep 'MicrobeAnnotator' | awk '{print $1, $2}'")
microbeannotator_kofam_date="Kofam data accessed on: $(cat "${databases_dir}/microbeannotator_db/kofam_data/date_accessed.txt")"
microbeannotator_uniprot_release="UniProt release version: $(cat "${databases_dir}/microbeannotator_db/protein_db/uniprot_release.txt")"

plasmidfinder_version="PlasmidFinder v2.1.6"
plasmidfinder_db_version=$(cat "$databases_dir/plasmidfinder_db/VERSION.txt")

rgi_version=$(apptainer exec "$straincascade_taxonomic_functional_analysis" bash -c "source activate rgi_env && rgi main --version 2>&1 | awk '{print "RGI", "v"$1}'")
rgi_db_version=$(awk '
    BEGIN {
        version = "";
        timestamp = "";
    }

    /"_version"/ {
        match($0, /"_version":[ ]*"[^"]+"/);
        version = substr($0, RSTART + 11, RLENGTH - 11);  # Adjust to remove _version:
        gsub(/"/, "", version);  # Remove remaining quotes
    }

    /"_timestamp"/ {
        match($0, /"_timestamp":[ ]*"[^"]+"/);
        timestamp = substr($0, RSTART + 13, RLENGTH - 13);  # Adjust to remove _timestamp:
        gsub(/"/, "", timestamp);  # Remove remaining quotes
    }

    END {
        print "v" version " released on " timestamp;
    }
' "$databases_dir/rgi_db/card.json")

resfinder_version=$(apptainer exec "$straincascade_taxonomic_functional_analysis" bash -c "source activate resfinder_env && python /opt/conda/envs/resfinder_env/lib/python3.8/site-packages/resfinder/run_resfinder.py -v 2>&1 | awk '{print "ResFinder", "v"$1}'")
resfinder_db_version=$(cat "$databases_dir/resfinder_db/VERSION")

dbcan_version="v4.1.4"
dbcan_dbCAN_sub_hmm_version=$(awk '/^DATE/ {print $2, $3, $4, $5, $6; exit}' "$databases_dir/dbcan3_db/db/dbCAN_sub.hmm")
dbcan_dbCAN_txt_version=$(awk '/^DATE/ {print $2, $3, $4, $5, $6; exit}' "$databases_dir/dbcan3_db/db/dbCAN.txt")

islandpath_version="v1.0.6"

virsorter2_version=$(apptainer exec "$straincascade_phage_detection" bash -c 'source activate virsorter2_env && virsorter --version 2>&1 | grep -oP "VirSorter \K[0-9]+\.[0-9]+\.[0-9]+"')
virsorter2_combined_hmm_version=$(awk '/^DATE/ {print $2, $3, $4, $5, $6; exit}' "$databases_dir/virsorter2_db/hmm/viral/combined.hmm")

deepvirfinder_version="1.0"

# Print the introductory message
echo "$sample_name was processed StrainCascade version $version"
echo "Current date (YYYY-MM-DD) and time (HH:MM:SS): $current_date_time"
echo ""
echo "###################"
echo ""
echo "System Information:"
echo "-------------------"
echo "Operating System: $os_name"
echo "OS Version: $os_version"
echo "Kernel Version: $kernel_version"
echo "Architecture: $architecture"
echo "Number of utilised CPUs: $threads"
echo ""
echo "###################"
echo ""
echo "Run and Tool Information:"
echo "-------------------"
echo "You ran the following modules of StrainCascade for the sample: $input_file with user-declared sequencing type: $sequencing_type"
echo "--- Format: Module SC ID: module name (module_script.sh). Main tools and databases used in the module ---"

SC1="SC1: Canu Correction and Trimming (StrainCascade_Canu_correct_trim.sh). Using Canu v$canu_version."
SC2="SC2: LJA Assembly (StrainCascade_LJA_assembly.sh). Using $lja_version."
SC3="SC3: SPAdes Assembly (StrainCascade_SPAdes_assembly.sh). Using $spades_version."
SC4="SC4: Canu Assembly (StrainCascade_Canu_assembly.sh). Using Canu v$canu_version."
SC5="SC5: Flye Assembly (StrainCascade_Flye_assembly.sh). Using Flye $flye_version."
SC6="SC6: Assembly Evaluation 1 (StrainCascade_assembly_evaluation1.sh). Using $quast_version for assembly assessment and StrainCascade's assembly selection algorithm for assembly selection with $straincascade_python_version."
SC7="SC7: MAC2.0 Assembly Merging (StrainCascade_MAC2_assembly_merging.sh). Using $mac2_version."
SC8="SC8: Assembly Evaluation 2 (StrainCascade_assembly_evaluation2.sh). Using $quast_version for assembly assessment and StrainCascade's assembly selection algorithm for assembly selection with $straincascade_python_version."
SC9="SC9: Circlator Circularisation (StrainCascade_Circlator_circularisation.sh). Using Circlator v$circlator_version."
SC10="SC10: Assembly Evaluation 3 (StrainCascade_assembly_evaluation3.sh). Using $quast_version for assembly assessment and StrainCascade's assembly selection algorithm for assembly selection with $straincascade_python_version."
SC11="SC11: Arrow Medaka Polishing (StrainCascade_arrow_medaka_polishing.sh). Using either $arrow_version or $medaka_version."
SC12="SC12: NGMLR BBMap Coverage (StrainCascade_NGMLR_BBMap_coverage.sh). Using either $ngmlr_version or $bbmap_version."
SC13="SC13: CheckM2 QC (StrainCascade_CheckM2_QC.sh). Using CheckM2 v$checkm2_version with CheckM2 database v$checkm2_db_version."
SC14="SC14: GTDB-Tk Taxonomy (StrainCascade_GTDB-Tk_taxonomy.sh). Using $gtdbtk_version with GTDB-Tk database $gtdbtk_db_version."
SC15="SC15: GTDB-Tk De Novo Tree (StrainCascade_GTDB-Tk_de_novo_tree.sh). Using $gtdbtk_version with GTDB-Tk database $gtdbtk_db_version."
SC16="SC16: Bakta Annotation (StrainCascade_Bakta_annotation.sh). Using $bakta_version with Bakta database: $bakta_db_version and ARMFinder update during run (amrfinder_update --force_update)."
SC17="SC17: Prokka Annotation (StrainCascade_Prokka_annotation.sh). Using $prokka_version."
SC18="SC18: MicrobeAnnotator Annotation (StrainCascade_MicrobeAnnotator_annotation.sh). Using $microbeannotator_version with $microbeannotator_kofam_date and $microbeannotator_uniprot_release."
SC19="SC19: PlasmidFinder Identification (StrainCascade_PlasmidFinder_identification.sh). Using $plasmidfinder_version with PlasmidFinder database v$plasmidfinder_db_version."
SC20="SC20: RGI Antimicrobial Resistance Identification (StrainCascade_RGI_antimicrobial_resistance_identification.sh). Using $rgi_version with CARD database $rgi_db_version."
SC21="SC21: ResFinder Antimicrobial Resistance Identification (StrainCascade_ResFinder_antimicrobial_resistance_identification.sh). Using $resfinder_version with ResFinder database v$resfinder_db_version."
SC22="SC22: dbCAN3 CAZymes Identification (StrainCascade_dbCAN3_CAZymes_identification.sh). Using dbCAN3 $dbcan_version with dbCAN.txt version $dbcan_dbCAN_txt_version and dbCAN_sub.hmm version $dbcan_dbCAN_sub_hmm_version."
SC23="SC23: IslandPath Genomic Islands Identification (StrainCascade_IslandPath_genomic_islands_identification.sh). Using IslandPath-DIMOB $islandpath_version."
SC24="SC24: VirSorter2 Phage Identification (StrainCascade_VirSorter2_phage_identification.sh). Using VirSorter2 v$virsorter2_version with combined.hmm version $virsorter2_combined_hmm_version."
SC25="SC25: DeepVirFinder Phage Identification (StrainCascade_DeepVirFinder_phage_identification.sh). Using DeepVirFinder v$deepvirfinder_version."
SC26="SC26: Data Integration (StrainCascade_data_integration.sh). Using $straincascade_r_version."

# Check each SC variable and add to selected_sc_info if it was run
for i in {1..26}; do
    sc_var="SC$i"
    if was_module_selected "${!sc_var}"; then
        selected_sc_info+=("${!sc_var}")
    fi
done

# Generate the StrainCascade_run_documentation.txt file
{
    echo "$sample_name was processed StrainCascade version $version"
    echo "Current date (YYYY-MM-DD) and time (HH:MM:SS): $current_date_time"
    echo ""
    echo "###################"
    echo ""
    echo "System Information:"
    echo "-------------------"
    echo "Operating System: $os_name"
    echo "OS Version: $os_version"
    echo "Kernel Version: $kernel_version"
    echo "Architecture: $architecture"
    echo ""
    echo "###################"
    echo ""
    echo "Run and Tool Information:"
    echo "-------------------"
    echo "You ran the following modules of StrainCascade:"
    for info in "${selected_sc_info[@]}"; do
        echo -e "$info"
    done
} > "$main_results_dir_abs/StrainCascade_run_documentation.txt"

echo "StrainCascade run documentation has been saved to ${main_results_dir_abs}/StrainCascade_run_documentation.txt"