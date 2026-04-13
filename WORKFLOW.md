# StrainCascade Workflow Documentation

> This document describes the module execution order, inputs, outputs, and dependencies as implemented in `StrainCascade_pipeline_wrapper.sh` and the individual module scripts.

------------------------------------------------------------------------

## Overview

StrainCascade contains **30 sequential modules (SC1–SC30)** plus a pre-processing input handler and a post-processing run summary. Modules are grouped into five stages:

| Stage | Modules | Description |
|------------------|----------------------|--------------------------------|
| **I. Genome Assembly & QC** | SC1–SC14 | Read correction, multi-assembler strategy, iterative evaluation, polishing, coverage, QC |
| **II. Taxonomic Classification** | SC15–SC16 | GTDB-Tk classification and de novo phylogenetic tree |
| **III. Genome Annotation** | SC17–SC20 | Bakta, Prokka, DeepFRI, and MicrobeAnnotator annotation |
| **IV. Functional Profiling** | SC21–SC30 | AMR, CAZymes, mobile genetic elements, phage, CRISPR-Cas, IS elements, data integration |
| **Summary** | — | Run summary report generation (MD + PDF) |

### Numbering Changes After Reviewer Feedback

The original manuscript (v1) described 27 modules (SC1–SC27). After reviewer feedback, the following changes were made, expanding the pipeline to 30 modules (SC1–SC30):

| Change | Details | Reviewer Comment |
|------------------|-------------------|-----------------------------------|
| **SC6 added** | Unicycler hybrid/long-read assembly | Reviewer requested hybrid short+long read support |
| **SC12 updated** | Assembly polishing now includes short-read polishing (Polypolish) in addition to long-read polishing (Medaka/Arrow/Racon) | Reviewer suggested optional short-read polishing step |
| **SC13 changed** | NGMLR replaced by minimap2 for read mapping; NGMLR fully removed from codebase and container | Stability improvement |
| **SC17 updated** | Bakta updated to v1.11.4 with db schema v6.0 | Reviewer noted Prokka is outdated since 2019 |
| **SC19 added** | DeepFRI deep-learning functional annotation | Reviewer specifically requested deep-learning comparison |
| **SC27 added** | geNomad phage and plasmid identification (replaces DeepVirFinder) | Reviewer recommended geNomad for mobile genetic elements |
| **SC30 added** | Data integration as explicit numbered module | Was implicit in v1 |

**Old → New numbering mapping (shifted due to SC6 insertion and SC19 insertion):**

| Manuscript v1 | Current (post-review) | Module |
|------------------------|------------------------|------------------------|
| SC1–SC5 | SC1–SC5 | Unchanged |
| — | SC6 | **New:** Unicycler Assembly |
| SC6 | SC7 | Assembly Evaluation 1 |
| SC7 | SC8 | MAC2 Assembly Merging |
| SC8 | SC9 | Assembly Evaluation 2 |
| SC9 | SC10 | Circlator Circularisation |
| SC10 | SC11 | Assembly Evaluation 3 |
| SC11 | SC12 | Assembly Polishing (updated: +Polypolish) |
| SC12 | SC13 | minimap2/BBMap Coverage (changed: NGMLR→minimap2) |
| SC13 | SC14 | CheckM2 QC |
| SC14 | SC15 | GTDB-Tk Taxonomy |
| SC15 | SC16 | GTDB-Tk De Novo Tree |
| SC16 | SC17 | Bakta Annotation (updated: v1.11.4, db v6.0) |
| SC17 | SC18 | Prokka Annotation |
| — | SC19 | **New:** DeepFRI Annotation |
| SC18 | SC20 | MicrobeAnnotator Annotation |
| SC19 | SC21 | PlasmidFinder Identification |
| SC20 | SC22 | AMRFinderPlus AMR Identification |
| SC21 | SC23 | ResFinder AMR Identification |
| SC22 | SC24 | dbCAN3 CAZymes Identification |
| SC23 | SC25 | IslandPath Genomic Islands Identification |
| SC25 | SC26 | VirSorter2 Phage Identification |
| SC26 (DeepVirFinder) | SC27 (geNomad) | Replaced: DeepVirFinder → geNomad |
| SC24 | SC28 | CRISPRCasFinder Identification |
| SC27 | SC29 | ISEScan IS Elements Identification |
| — | SC30 | **New:** Data Integration |

------------------------------------------------------------------------

## Directory Structure

All output is placed under `StrainCascade_output_<sample_name>/` with the following structure:

```         
StrainCascade_output_<sample>/
└── <sample>_StrainCascade_main_results/
    ├── 00_sequencing_reads/
    ├── 01_genome_assembly/
    ├── 02_taxonomic_classification/
    ├── 03_genome_annotation/
    ├── 04_functional_analysis/
    ├── 05_results_integration/
    └── pipeline_logs/
```

------------------------------------------------------------------------

## Pre-Processing: Input File Handler

**Script:** `StrainCascade_input_file_handler.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Routes and converts input files to the appropriate directory based on input type. Converts BAM to FASTQ if needed. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SEQUENCING_READS_DIR`, `GENOME_ASSEMBLY_DIR`, `INPUT_TYPE`, `SHORT_READS_R1`, `SHORT_READS_R2`, `THREADS` |
| **Inputs** | Raw sequencing file (FASTQ/FASTQ.GZ/FASTA/BAM) or pre-assembled genome (FASTA/FA/FNA) |
| **Outputs** | Processed reads in `00_sequencing_reads/` or assembly in `01_genome_assembly/`; returns tab-separated: `processed_input`, `bam_file`, `processed_short_r1`, `processed_short_r2` |
| **Container** | `straincascade_genome_assembly*.sif` |

------------------------------------------------------------------------

## Stage I: Genome Assembly & Quality Control (SC1–SC14)

### SC1 — Canu Correction and Trimming

**Script:** `StrainCascade_Canu_correct_trim.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Performs read correction and/or trimming using Canu. Skipped for HiFi data. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Raw long reads (from input handler) |
| **Outputs** | Corrected/trimmed reads in `Canu_correct_trim_results/` |
| **Supported types** | `pacbio-raw`, `pacbio-corr`, `nano-raw`, `nano-corr`, `nano-hq` (skips `pacbio-hifi`) |
| **Container** | `straincascade_genome_assembly*.sif` |
| **Deterministic mode** | Threads set to 1; entropy file bound to `/dev/random` and `/dev/urandom` |

------------------------------------------------------------------------

### SC2 — LJA Assembly

**Script:** `StrainCascade_LJA_assembly.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome assembly using La Jolla Assembler (LJA). Optimized for PacBio HiFi. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Long reads (raw or corrected/trimmed from SC1) |
| **Outputs** | `LJA_assembly_results/`; final assembly copied as `<sample>_assembly_lja.fasta` → `01_genome_assembly/` |
| **Validation** | Genome size check (1.3–10 Mbp range) |
| **Container** | `straincascade_lja_genome_assembly*.sif` |
| **Deterministic mode** | Threads=1, deterministic entropy file |

------------------------------------------------------------------------

### SC3 — SPAdes Assembly

**Script:** `StrainCascade_SPAdes_assembly.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome assembly using SPAdes with `-s` flag for long reads. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Long reads |
| **Outputs** | `SPAdes_assembly_results/`; final assembly `<sample>_assembly_spades.fasta` → `01_genome_assembly/` |
| **Container** | `straincascade_genome_assembly*.sif` |
| **Deterministic mode** | Threads=1, deterministic entropy file |

------------------------------------------------------------------------

### SC4 — Canu Assembly

**Script:** `StrainCascade_Canu_assembly.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome assembly using Canu. Uses corrected/trimmed reads from SC1 (`-trimmed -corrected` flags) for non-HiFi data. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Corrected/trimmed reads from SC1 (or raw HiFi reads) |
| **Outputs** | `Canu_assembly_results/`; final assembly `<sample>_assembly_canu.fasta` → `01_genome_assembly/` |
| **Container** | `straincascade_genome_assembly*.sif` |
| **Deterministic mode** | Threads=1, deterministic entropy file |

------------------------------------------------------------------------

### SC5 — Flye Assembly

**Script:** `StrainCascade_Flye_assembly.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome assembly using Flye. Uses informed genome size estimate if available, otherwise defaults to 4.5m. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Long reads; optionally `informed_genome_size_estimation.txt` |
| **Outputs** | `Flye_assembly_results/`; final assembly `<sample>_assembly_flye.fasta` → `01_genome_assembly/` |
| **Container** | `straincascade_genome_assembly*.sif` |
| **Deterministic mode** | Threads=1, deterministic entropy file |

------------------------------------------------------------------------

### SC6 — Unicycler Assembly *(added post-review)*

**Script:** `StrainCascade_Unicycler_assembly.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome assembly using Unicycler in **hybrid mode** (long + short reads) or **long-read-only mode**. Added in response to reviewer feedback requesting hybrid short+long read assembly support. |
| **Parameters (14)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE`, `SHORT_READS_R1`, `SHORT_READS_R2` |
| **Inputs** | Long reads; optionally paired-end short reads (R1, R2) |
| **Outputs** | `Unicycler_assembly_results/`; final assembly `<sample>_assembly_uc.fasta` → `01_genome_assembly/` |
| **Container** | `straincascade_genome_assembly*.sif` |
| **Note** | Cleans existing output directory before running |

------------------------------------------------------------------------

### SC7 — Assembly Evaluation 1

**Script:** `StrainCascade_assembly_evaluation1.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | First iterative evaluation. Runs QUAST on all assemblies, selects the best one. |
| **Parameters (10)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `SELECTION_ALGORITHM` |
| **Inputs** | All `*_assembly_*.fasta` files in `01_genome_assembly/` (excludes `*_best_ev*`, `*_circularised*`, `*_mac2.0_merged*`) |
| **Outputs** | `StrainCascade_assembly_evaluation/`; QUAST report `<sample>_quast_assembly_evaluation.tsv`; best assembly tagged as `*_best_ev1.fasta` |
| **Selection** | Python script (`StrainCascade_assembly_selection_contig.py` or `StrainCascade_assembly_selection_continuity.py`) based on `SELECTION_ALGORITHM` |
| **Containers** | `straincascade_assembly_qc_refinement*.sif` (QUAST), `python_3.12.4*.sif` (selection) |

------------------------------------------------------------------------

### SC8 — MAC2 Assembly Merging

**Script:** `StrainCascade_MAC2_assembly_merging.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Merges multiple assemblies using MAC2 (Meta-Assembly Consensus). Uses `best_ev1` as reference. Requires minimum 5 assemblies. |
| **Parameters (8)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `GENOME_ASSEMBLY_DIR` |
| **Inputs** | All assemblies in `01_genome_assembly/` + `best_ev1` as reference |
| **Outputs** | `MAC2_optimization_results/`; merged assembly `<sample>_mac2.0_merged_assembly.fasta` → `01_genome_assembly/` |
| **Constraints** | Max runtime: 200 hours; filter threshold: 24 hours |
| **Container** | `straincascade_assembly_qc_refinement*.sif` |

------------------------------------------------------------------------

### SC9 — Assembly Evaluation 2

**Script:** `StrainCascade_assembly_evaluation2.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Second iterative evaluation. Re-runs QUAST on all assemblies including the MAC2 merged assembly. |
| **Parameters (10)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `SELECTION_ALGORITHM` |
| **Inputs** | All assemblies (excludes `*_best_ev2*`, `*_best_ev3*`, `*_circularised*`) |
| **Outputs** | QUAST report TSV; best assembly tagged as `*_best_ev2.fasta` |
| **Containers** | `straincascade_assembly_qc_refinement*.sif`, `python_3.12.4*.sif` |

------------------------------------------------------------------------

### SC10 — Circlator Circularisation

**Script:** `StrainCascade_Circlator_circularisation.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome circularization using Circlator. Prioritizes: `best_ev2` \> `best_ev1` \> any `.fasta`. |
| **Parameters (11)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `GENOME_ASSEMBLY_DIR`, `RESULTS_INTEGRATION_DIR`, `VERSION` |
| **Inputs** | Best assembly from SC9 (or SC7 fallback); raw reads for remapping |
| **Outputs** | `Circlator_circularisation_results/`; circularised assembly `*_circularised.fasta` → `01_genome_assembly/`; QS files → `05_results_integration/qs_files/` |
| **Constraints** | Max runtime: 200 hours |
| **Containers** | `straincascade_assembly_qc_refinement*.sif` (Circlator), `r_4.4.1*.sif` (R processing) |

------------------------------------------------------------------------

### SC11 — Assembly Evaluation 3

**Script:** `StrainCascade_assembly_evaluation3.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Final (third) iterative evaluation. Evaluates all assemblies including circularised version. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `RESULTS_INTEGRATION_DIR`, `VERSION`, `SELECTION_ALGORITHM` |
| **Inputs** | All assemblies (excludes `*_best_ev3*`) |
| **Outputs** | QUAST report TSV; best assembly tagged as `*_best_ev3.fasta`; QS files → `05_results_integration/qs_files/` |
| **Containers** | `straincascade_assembly_qc_refinement*.sif`, `python_3.12.4*.sif`, `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC12 — Assembly Polishing *(updated post-review)*

**Script:** `StrainCascade_assembly_polishing.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Unified polishing module with both long-read and short-read polishing (2 rounds each). Updated in response to reviewer feedback to include optional short-read polishing with Polypolish. |
| **Parameters (17)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `BAM_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `REPRODUCIBILITY_MODE`, `RESULTS_INTEGRATION_DIR`, `VERSION`, `SHORT_READS_R1`, `SHORT_READS_R2` |
| **Inputs** | Best assembly (analysis assembly); raw long reads; optionally BAM file; optionally paired-end short reads |
| **Long-read polishing** | Arrow (PacBio with BAM), Racon (PacBio without BAM), or Medaka (Nanopore) |
| **Short-read polishing** | Polypolish (if paired-end short reads provided) |
| **Outputs** | `assembly_polishing_results/longread_polishing/`; `assembly_polishing_results/shortread_polishing/` |
| **Container** | `straincascade_assembly_qc_refinement*.sif` |

------------------------------------------------------------------------

### SC13 — minimap2/BBMap Coverage *(changed post-review)*

**Script:** `StrainCascade_minimap2_BBMap_coverage.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Read mapping and coverage analysis. minimap2 replaces the previously used NGMLR for long-read mapping. BBMap is used for short-read mapping and coverage statistics (pileup.sh). |
| **Parameters (11)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `INPUT_FILE`, `OUTPUT_DIR`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `THREADS`, `GENOME_ASSEMBLY_DIR` |
| **Inputs** | Raw reads; analysis assembly from `01_genome_assembly/` |
| **Outputs** | `minimap2_BBMap_coverage_results/` (SAM/BAM files, coverage statistics) |
| **Mapping presets** | `map-hifi` (PacBio HiFi), `map-pb` (PacBio), `map-ont` (Nanopore); BBMap for short reads (≤600 bp) |
| **Container** | `straincascade_assembly_qc_refinement*.sif` |

------------------------------------------------------------------------

### SC14 — CheckM2 Quality Control

**Script:** `StrainCascade_CheckM2_QC.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome completeness and contamination assessment using CheckM2. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `CheckM2_QC_results/`; quality report `<sample>_checkm2_quality_report.tsv` → `01_genome_assembly/`; QS files → `05_results_integration/qs_files/` |
| **Database** | `uniref100.KO.1.dmnd` |
| **Containers** | `straincascade_assembly_qc_refinement*.sif` (CheckM2), `r_4.4.1*.sif` (R processing) |

------------------------------------------------------------------------

## Stage II: Taxonomic Classification (SC15–SC16)

### SC15 — GTDB-Tk Taxonomy

**Script:** `StrainCascade_GTDB-Tk_taxonomy.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Taxonomic classification using GTDB-Tk `classify_wf`. |
| **Parameters (13)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `TAXONOMIC_CLASS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `GTDB-Tk_taxonomy_results/`; phylum name file; QS files → `05_results_integration/qs_files/` |
| **Database** | GTDB reference database |
| **Containers** | `straincascade_taxonomic_functional_analysis*.sif` (conda env `gtdbtk_env`), `python_3.12.4*.sif`, `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC16 — GTDB-Tk De Novo Tree

**Script:** `StrainCascade_GTDB-Tk_de_novo_tree.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | De novo phylogenetic tree construction. Places assembled genome with optional external assemblies. Uses phylum from SC15 as outgroup. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `EXTERNAL_ASSEMBLY_DIR`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `TAXONOMIC_CLASS_DIR`, `DATABASES_DIR` |
| **Inputs** | Analysis assembly; GTDB-Tk taxonomy results from SC15; optional external assemblies |
| **Outputs** | `GTDB-Tk_de_novo_tree_results/` (phylogenetic tree files) |
| **Container** | `straincascade_taxonomic_functional_analysis*.sif` (conda env `gtdbtk_env`) |

------------------------------------------------------------------------

## Stage III: Genome Annotation (SC17–SC20)

### SC17 — Bakta Annotation *(updated post-review)*

**Script:** `StrainCascade_Bakta_annotation.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome annotation with Bakta v1.11.4 (db schema v6.0). Updated in response to reviewer comment that Prokka is outdated. Uses taxonomic info from SC15 for organism/genus context. |
| **Parameters (17)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `TAXONOMIC_CLASS_DIR`, `GENOME_ANNOTATION_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `LOCUS_TAG`, `FORCE_OVERWRITE`, `VERSION`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Analysis assembly; taxonomic classification from SC15 |
| **Outputs** | `Bakta_annotation_results/`; annotated files (`*_annotation_bakta.tsv`, `.ffn`, `.faa`, `.gbff`, `.png`, `.svg`) → `03_genome_annotation/`; QS files |
| **Containers** | `straincascade_genome_annotation*.sif`, `r_4.4.1*.sif` |
| **Deterministic mode** | Threads=1, entropy file |

------------------------------------------------------------------------

### SC18 — Prokka Annotation

**Script:** `StrainCascade_Prokka_annotation.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genome annotation with Prokka. Retained alongside Bakta for cross-validation; reviewer noted Prokka is outdated since 2019. |
| **Parameters (15)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `INITIAL_THREADS`, `GENOME_ASSEMBLY_DIR`, `GENOME_ANNOTATION_DIR`, `RESULTS_INTEGRATION_DIR`, `LOCUS_TAG`, `FORCE_OVERWRITE`, `VERSION`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `Prokka_annotation_results/`; annotated files (`*_annotation_prokka.gff`, `.faa`, `.ffn`, `.tsv`) → `03_genome_annotation/`; QS files |
| **Containers** | `straincascade_genome_annotation*.sif`, `r_4.4.1*.sif` |
| **Deterministic mode** | Threads=1, entropy file |

------------------------------------------------------------------------

### SC19 — DeepFRI Annotation *(added post-review)*

**Script:** `StrainCascade_DeepFRI_annotation.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Deep-learning-based functional annotation using DeepFRI across 4 Gene Ontology categories. Added in direct response to reviewer recommendation for deep-learning model comparison. |
| **Parameters (11)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `GENOME_ANNOTATION_DIR`, `RESULTS_INTEGRATION_DIR`, `VERSION`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Protein sequences (`.faa`) from Bakta (SC17) or Prokka (SC18) in `03_genome_annotation/` |
| **Outputs** | `DeepFRI_annotation_results/`; per-ontology predictions (`*_deepfri_mf*`, `*_deepfri_bp*`, `*_deepfri_cc*`, `*_deepfri_ec*`); QS files |
| **Ontologies** | MF (Molecular Function), BP (Biological Process), CC (Cellular Component), EC (Enzyme Commission) |
| **Containers** | `straincascade_genome_annotation*.sif` (conda env `deepfri_env`), `r_4.4.1*.sif` |
| **Deterministic mode** | `PYTHONHASHSEED=42` |

------------------------------------------------------------------------

### SC20 — MicrobeAnnotator Annotation

**Script:** `StrainCascade_MicrobeAnnotator_annotation.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Functional annotation with MicrobeAnnotator for KEGG metabolic modules. |
| **Parameters (14)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `INITIAL_THREADS`, `GENOME_ANNOTATION_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION`, `REPRODUCIBILITY_MODE` |
| **Inputs** | Protein sequences (`.faa`) — priority: Bakta → Prokka → any `.faa` fallback |
| **Outputs** | `MicrobeAnnotator_annotation_results/`; QS files |
| **Containers** | `straincascade_genome_annotation*.sif`, `r_4.4.1*.sif` |
| **Deterministic mode** | Supported |

------------------------------------------------------------------------

## Stage IV: Functional Profiling (SC21–SC30)

### SC21 — PlasmidFinder Identification

**Script:** `StrainCascade_PlasmidFinder_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Plasmid identification using PlasmidFinder with KMA-indexed database. |
| **Parameters (11)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `GENOME_ASSEMBLY_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `PlasmidFinder_plasmid_identification_results/`; `*_plasmidfinder_*.txt`, `*.json`; QS files |
| **Containers** | `straincascade_assembly_qc_refinement*.sif` (conda env `plasmidfinder_env`, `blastn`), `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC22 — AMRFinderPlus Antimicrobial Resistance Identification

**Script:** `StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | AMR gene and point-mutation detection using AMRFinderPlus. Determines organism from GTDB-Tk results and checks against 28 supported organisms for species-specific point-mutation detection. |
| **Parameters (14)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `TAXONOMIC_CLASS_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly; taxonomic classification from SC15 |
| **Outputs** | `AMRFinderPlus_AMR_identification_results/`; `*_identification_amrfinderplus.*`, `*_mutations_amrfinderplus.*`; QS files |
| **Containers** | `straincascade_taxonomic_functional_analysis*.sif`, `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC23 — ResFinder Antimicrobial Resistance Identification

**Script:** `StrainCascade_ResFinder_antimicrobial_resistance_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | AMR identification using ResFinder. Detects sequencing type for Nanopore-specific flags. |
| **Parameters (14)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `SEQUENCING_TYPE`, `GENOME_ASSEMBLY_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `ResFinder_AMR_identification_results/` (alignments, acquired/point mutations); QS files |
| **Databases** | `resfinder_db`, `pointfinder_db`, `disinfinder_db` |
| **Containers** | `straincascade_taxonomic_functional_analysis*.sif` (conda env `resfinder_env`), `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC24 — dbCAN3 CAZymes Identification

**Script:** `StrainCascade_dbCAN3_CAZymes_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Carbohydrate-active enzyme (CAZyme) identification using dbCAN3 `run_dbcan`. Uses three integrated tools: HMMER, DIAMOND, and dbCAN_sub (note: Hotpep was replaced by dbCAN-sub in 2022). |
| **Parameters (13)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `dbCAN3_CAZymes_identification_results/` (HMMER, DIAMOND, dbCAN_sub, CGC results); QS files |
| **Thresholds** | HMMER: e-value \<1e-15, coverage \>0.35; dbCAN_sub/DIAMOND: e-value \<1e-102 |
| **Containers** | `straincascade_taxonomic_functional_analysis*.sif`, `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC25 — IslandPath Genomic Islands Identification

**Script:** `StrainCascade_IslandPath_genomic_islands_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Genomic island identification using IslandPath-DIMOB. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `GENOME_ASSEMBLY_DIR`, `GENOME_ANNOTATION_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `VERSION` |
| **Inputs** | GenBank file from Bakta (`.gbff`) → Prokka (`.gbk`) → any `.gbk`/`.gbff`/`.embl` fallback |
| **Outputs** | `IslandPath_genomic_island_identification_results/`; `*_results_genomic_island_identification_islandpath.gff3`; QS files |
| **Containers** | `straincascade_taxonomic_functional_analysis*.sif` (conda env `islandpath_env`, `Dimob.pl`), `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC26 — VirSorter2 Phage Identification

**Script:** `StrainCascade_VirSorter2_phage_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Viral/phage identification using VirSorter2 `run ... all` with `--min-length 200`. |
| **Parameters (13)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `VirSorter2_phage_identification_results/`; `final-viral-combined.fa`, `final-viral-score.tsv`, `final-viral-boundary.tsv` (copied to `04_functional_analysis/` with sample prefix) |
| **Container** | `straincascade_crisprcas_phage_is_elements*.sif` (conda env `virsorter2_env`), `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC27 — geNomad Phage and Plasmid Identification *(added post-review)*

**Script:** `StrainCascade_geNomad_phage_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Phage and plasmid identification using geNomad `end-to-end`. Added in response to reviewer recommendation to include geNomad for comprehensive mobile genetic element detection (replaces DeepVirFinder from manuscript v1). |
| **Parameters (13)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `DATABASES_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `geNomad_phage_identification_results/`; virus summary/genes/sequences TSVs and FNA/FAA files → `04_functional_analysis/`; QS files |
| **Settings** | Min score: 0.7; splits: 8 |
| **Container** | `straincascade_crisprcas_phage_is_elements*.sif` (conda env `genomad_env`), `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC28 — CRISPRCasFinder Identification

**Script:** `StrainCascade_CRISPRCasFinder_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | CRISPR-Cas system identification using `CRISPRCasFinder.pl`. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `CRISPRCasFinder_identification_results/`; `Result_*/` directory, `result.json` (also copied to `04_functional_analysis/`) |
| **Flags** | `-cas -keep` |
| **Container** | `straincascade_crisprcas_phage_is_elements*.sif` (conda env `crisprcasfinder`) |

------------------------------------------------------------------------

### SC29 — ISEScan IS Elements Identification

**Script:** `StrainCascade_ISEScan_IS_elements_identification.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Insertion Sequence (IS) element identification using ISEScan. |
| **Parameters (12)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `OUTPUT_DIR`, `SAMPLE_NAME`, `THREADS`, `GENOME_ASSEMBLY_DIR`, `FUNCTIONAL_ANALYSIS_DIR`, `RESULTS_INTEGRATION_DIR`, `VERSION` |
| **Inputs** | Analysis assembly from `01_genome_assembly/` |
| **Outputs** | `ISEScan_identification_results/`; `.gff`, `.tsv`, `.is.fna`, `.orf.fna`, `.orf.faa` (copied to `04_functional_analysis/`); QS files |
| **Container** | `straincascade_crisprcas_phage_is_elements*.sif` (conda env `isescan_env`, `isescan.py`), `r_4.4.1*.sif` |

------------------------------------------------------------------------

### SC30 — Data Integration

**Script:** `StrainCascade_data_integration.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Integrates results from all upstream modules into comprehensive R-based reports and datasets. Handles contig-name normalization. |
| **Parameters (9)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_DIR`, `GENOME_ASSEMBLY_DIR`, `RESULTS_DIR`, `VERSION`, `SAMPLE_NAME` |
| **Inputs** | All results from `05_results_integration/qs_files/` (QS files from all modules) |
| **Outputs** | Contig name mapping; integrated QS data files; `StrainCascade_Results*.RData`; `StrainCascade_Analysis*.html` in `05_results_integration/` |
| **Container** | `r_4.4.1*.sif` |

------------------------------------------------------------------------

## Post-Processing: Run Summary

**Script:** `StrainCascade_run_summary.sh`

|  |  |
|------------------------------------|------------------------------------|
| **Purpose** | Generates a detailed run summary document (Markdown + secure PDF). Collects system info, tool versions, runtime duration, and container hashes. |
| **Parameters (20)** | `SCRIPT_DIR`, `LOGS_DIR`, `LOG_NAME`, `UTILS_FILE`, `APPTAINER_IMAGES_DIR`, `INPUT_FILE`, `INPUT_TYPE`, `SAMPLE_NAME`, `SEQUENCING_TYPE`, `THREADS`, `DATABASES_DIR`, `MAIN_RESULTS_DIR`, `GENOME_ASSEMBLY_DIR`, `RESULTS_INTEGRATION_DIR`, `SELECTED_MODULES`, `VERSION`, `SELECTION_ALGORITHM`, `REPRODUCIBILITY_MODE`, `START_TIME`, `STOP_TIME` |
| **Inputs** | All pipeline metadata, container image directory |
| **Outputs** | Run summary (`.md` + `.pdf`) with system info, tool versions, runtime stats |
| **Containers** | All container images are queried for version extraction |

------------------------------------------------------------------------

## Module Execution Flow Diagram

```         
Input File Handler
       │
       ▼
┌─── Stage I: Genome Assembly & QC ───────────────────────────────┐
│                                                                  │
│  SC1 Canu Correct/Trim                                          │
│       │                                                          │
│       ▼                                                          │
│  SC2 LJA ──┐                                                    │
│  SC3 SPAdes─┤                                                    │
│  SC4 Canu ──┼── (parallel assemblers)                            │
│  SC5 Flye ──┤                                                    │
│  SC6 Unicycler┘                                                  │
│       │                                                          │
│       ▼                                                          │
│  SC7 Assembly Evaluation 1 ──→ best_ev1.fasta                   │
│       │                                                          │
│       ▼                                                          │
│  SC8 MAC2 Assembly Merging (needs ≥5 assemblies)                │
│       │                                                          │
│       ▼                                                          │
│  SC9 Assembly Evaluation 2 ──→ best_ev2.fasta                   │
│       │                                                          │
│       ▼                                                          │
│  SC10 Circlator Circularisation                                  │
│       │                                                          │
│       ▼                                                          │
│  SC11 Assembly Evaluation 3 ──→ best_ev3.fasta (analysis asm)   │
│       │                                                          │
│       ▼                                                          │
│  SC12 Assembly Polishing (long-read + optional short-read)       │
│       │                                                          │
│       ▼                                                          │
│  SC13 minimap2/BBMap Coverage                                    │
│       │                                                          │
│       ▼                                                          │
│  SC14 CheckM2 QC                                                │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
       │
       ▼
┌─── Stage II: Taxonomic Classification ──────────────────────────┐
│  SC15 GTDB-Tk Taxonomy                                          │
│       │                                                          │
│       ▼                                                          │
│  SC16 GTDB-Tk De Novo Tree                                      │
└──────────────────────────────────────────────────────────────────┘
       │
       ▼
┌─── Stage III: Genome Annotation ────────────────────────────────┐
│  SC17 Bakta Annotation                                          │
│       │                                                          │
│       ▼                                                          │
│  SC18 Prokka Annotation                                         │
│       │                                                          │
│       ▼                                                          │
│  SC19 DeepFRI Annotation (uses .faa from SC17/SC18)             │
│       │                                                          │
│       ▼                                                          │
│  SC20 MicrobeAnnotator Annotation (uses .faa from SC17/SC18)    │
└──────────────────────────────────────────────────────────────────┘
       │
       ▼
┌─── Stage IV: Functional Profiling ──────────────────────────────┐
│  SC21 PlasmidFinder                                             │
│  SC22 AMRFinderPlus (uses taxonomy from SC15)                   │
│  SC23 ResFinder                                                 │
│  SC24 dbCAN3 CAZymes                                            │
│  SC25 IslandPath (uses .gbff from SC17 or .gbk from SC18)      │
│  SC26 VirSorter2                                                │
│  SC27 geNomad                                                   │
│  SC28 CRISPRCasFinder                                           │
│  SC29 ISEScan                                                   │
│       │                                                          │
│       ▼                                                          │
│  SC30 Data Integration                                          │
└──────────────────────────────────────────────────────────────────┘
       │
       ▼
   Run Summary (MD + PDF)
```

------------------------------------------------------------------------

## Assembly Input Mode

When `INPUT_TYPE=assembly`, the following modules are automatically **skipped** (regex match in wrapper):

-   All assemblers: SC1 (Canu correct/trim), SC2–SC6 (LJA, SPAdes, Canu, Flye, Unicycler)
-   SC7, SC9 (Assembly Evaluation 1 & 2)
-   SC8 (MAC2)
-   SC10 (Circlator)
-   SC12 (Assembly Polishing)
-   SC13 (Coverage)

The pipeline skips to SC11 (Assembly Evaluation 3), then proceeds to SC14 (CheckM2 QC) and onward.

------------------------------------------------------------------------

## Execution Modes

| Mode | Modules Included |
|--------------------|----------------------------------------------------|
| **minimal** | SC3, SC7, SC14, SC15, SC17, SC30 |
| **efficient** | SC1, SC2, SC3, SC7, SC8, SC9, SC10, SC11, SC12, SC14, SC15, SC17, SC21, SC30 |
| **standard** | SC1–SC15, SC17–SC30 (all except SC16 GTDB-Tk De Novo Tree) |
| **comprehensive** | SC1–SC30 (all modules including SC16) |

### Module Bundles

| Bundle         | Modules                      |
|----------------|------------------------------|
| **assembly**   | SC1–SC6, SC7–SC14            |
| **annotation** | SC17, SC18, SC19, SC20       |
| **functional** | SC15, SC22, SC23, SC24, SC25 |
| **phage**      | SC26, SC27, SC28, SC29       |

------------------------------------------------------------------------

## Container Images

| Container | Modules |
|---------------------------------------|---------------------------------|
| `straincascade_genome_assembly*.sif` | SC1–SC6, Input Handler |
| `straincascade_lja_genome_assembly*.sif` | SC2 (LJA only) |
| `straincascade_assembly_qc_refinement*.sif` | SC7–SC14, SC21 |
| `straincascade_genome_annotation*.sif` | SC17–SC20 |
| `straincascade_taxonomic_functional_analysis*.sif` | SC15–SC16, SC22–SC25 |
| `straincascade_crisprcas_phage_is_elements*.sif` | SC26–SC29 |
| `r_4.4.1*.sif` | R processing across most modules |
| `python_3.12.4*.sif` | Assembly selection (SC7, SC9, SC11), SC15 |
| `straincascade_document_processing*.sif` | Run Summary |

------------------------------------------------------------------------

## Key Dependencies Between Modules

| Downstream Module | Depends On |
|------------------------------------|------------------------------------|
| SC4 (Canu Assembly) | SC1 (corrected/trimmed reads) |
| SC7 (Evaluation 1) | SC2–SC6 (assemblies) |
| SC8 (MAC2) | SC7 (best_ev1) + ≥5 assemblies |
| SC9 (Evaluation 2) | SC8 (merged assembly) |
| SC10 (Circlator) | SC9 (best_ev2) or SC7 (best_ev1) fallback |
| SC11 (Evaluation 3) | SC10 (circularised assembly) |
| SC12 (Polishing) | Analysis assembly (best_ev3 or best_ev2 or best_ev1) |
| SC15 (GTDB-Tk) | Analysis assembly |
| SC16 (De Novo Tree) | SC15 (phylum/taxonomy results) |
| SC17 (Bakta) | Analysis assembly + SC15 (taxonomy, optional) |
| SC18 (Prokka) | Analysis assembly |
| SC19 (DeepFRI) | SC17 or SC18 (`.faa` protein file) |
| SC20 (MicrobeAnnotator) | SC17 or SC18 (`.faa` protein file) |
| SC22 (AMRFinderPlus) | SC15 (taxonomy for organism matching) |
| SC25 (IslandPath) | SC17 (`.gbff`) or SC18 (`.gbk`) |

------------------------------------------------------------------------