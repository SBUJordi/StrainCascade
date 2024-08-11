# Central configuration file of the BIOMES WGseq pipeline
# config.sh - Version 1.0.0
# Author: Sebastian Bruno Ulrich Jordi

# The default is that the pipeline submits an array job with 1 node per sample. Settings can be changed in submit_jobs_in_batches.sh

# Set this to "yes" if you only want the key results (intermediate or additional results will be deleted). If you choose "no", some other value or let this variable undefined, all results will be kept.
only_main_results="no"

# Set number of threads used for parallelisation where possible
threads=32

# Modify these paths once here according to your setup
apptainer_images_path="/storage/workspaces/dbmr_mg/nbmisi/sebastian/isoBIOMES/container_images"
straincascade_genome_assembly_path="${apptainer_images_path}/straincascade_genome_assembly_v1.sif"
straincascade_assembly_qc_refinement_path="${apptainer_images_path}/straincascade_assembly_qc_refinement_v1.sif"
straincascade_genome_annotation_path="${apptainer_images_path}/straincascade_genome_annotation_v1.sif"

# Modify these paths once here according to your setup
miniconda_path="/storage/workspaces/dbmr_mg/nbmisi/sebastian/miniconda3"
lja_path="${miniconda_path}/envs/LJA/LJA/bin/lja"
spades_path="${miniconda_path}/envs/SPAdes/bin/spades.py"
cisa_merge_py_path="${miniconda_path}/envs/cisa/CISA1.3/Merge.py"
cisa_cisa_py_path="${miniconda_path}/envs/cisa/CISA1.3/CISA.py"
bakta_db_path="/storage/homefs/sj21o643/Bakta_db/db"
plasmidfinder_path="${miniconda_path}/envs/PlasmidFinder/plasmidfinder/plasmidfinder.py"
plasmidfinder_db_path="/storage/homefs/sj21o643/plasmidfinder_db"
dbcan_db_path="/storage/homefs/sj21o643/dbCAN_db"
dbcan_path="${miniconda_path}/envs/dbCAN3/lib/python3.8/site-packages/dbcan/cli/run_dbcan.py"
resfinder_path="${miniconda_path}/envs/ResFinder/lib/python3.12/site-packages/resfinder/run_resfinder.py"
resfinder_db_path="${miniconda_path}/envs/ResFinder/db_resfinder" # ResFinder for acquired resistance genes
pointfinder_db_path="${miniconda_path}/envs/ResFinder/db_pointfinder" # PointFinder for chromosomal mutations
islandpath_path="${miniconda_path}/envs/IslandPath-DIMOB/opt/islandpath/Dimob.pl"

# Read type for SPAdes
# Set this for SPAdes to either 'PacBio CLR' when  or 'single' (for PacBio CCS use single); further options are possible, see SPAdes manual (pipeline adaptions will be necessary)
read_type="single"

# Sequencing type for Canu; options [-pacbio|-nanopore|-pacbio-hifi]
sequencing_type_canu="-pacbio-hifi"

# Sequencing type for Flye; options (--pacbio-raw | --pacbio-corr | --pacbio-hifi | --nano-raw | --nano-corr | --nano-hq )
sequencing_type_flye="pacbio-hifi"

# Sequencing type for ngmlr; options (-x pacbio | -x ont)
sequencing_type_ngmlr="pacbio"

# Minimum coverage for ResFinder. At least 60% of a resistance gene must be covered by the sequencing reads for it to be reported.
# This helps to avoid false positives from partial gene matches. The specific value of 0.6 is chosen as a balance between sensitivity (detecting all true resistance genes) and specificity (avoiding false positives).
resfinder_coverage=0.6

# Identity threshold for ResFinder. The sequencing reads must have at least 80% identity to a resistance gene for it to be reported.
# This helps to ensure that the reported genes are indeed the resistance genes and not some other genes that are slightly similar. The specific value of 0.8 is chosen for similar reasons as the -l parameter. It's a balance between sensitivity and specificity.
resfinder_identity=0.8

# List of selected modules
# Uncomment (= remove "#" before a line) the modules you want to include in the pipeline
# Always run BIOMES_assembly_evaluation.sh if you plan any downstream modules (classification, annotation etc.)
selected_modules=(
  "BIOMES_LJA_assembly.sh"
  "BIOMES_SPAdes_assembly.sh"
  "BIOMES_Canu_assembly.sh"
  "BIOMES_Flye_assembly.sh"
  "BIOMES_CISA_merging.sh"
  "BIOMES_assembly_evaluation.sh"
  "BIOMES_Circlator_circularisation.sh"
  "BIOMES_assembly_evaluation2.sh"
  "BIOMES_BBMap_coverage.sh"
  "BIOMES_gtdbtk_taxonomy.sh"
  "BIOMES_gtdbtk_de_novo_tree.sh"
  "BIOMES_prokka_annotation.sh"
  "BIOMES_Bakta_annotation.sh"
  "BIOMES_CheckM2_QC.sh"
  "BIOMES_PlasmidFinder_identification.sh"
  "BIOMES_MicrobeAnnotator_annotation.sh"
  "BIOMES_RGI_antimicrobial_resistance_identification.sh"
  "BIOMES_ResFinder_antimicrobial_resistance_identification.sh"
  "BIOMES_dbCAN_CAZymes_identification.sh"
  "BIOMES_IslandPath_genomic_islands_identification.sh"
  "BIOMES_results_summary.sh"
  #"BIOMES_result_processing.sh"
)