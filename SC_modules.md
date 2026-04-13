# StrainCascade Module Reference

Complete reference of all 30 SC modules, derived from `StrainCascade_pipeline_wrapper.sh`.

| Module | Script | Description | Key Inputs |
|----------------|----------------|---------------------|-------------------|
| **SC1** | `StrainCascade_Canu_correct_trim.sh` | Long-read correction & trimming (Canu) | Raw reads, sequencing type, threads, reproducibility mode |
| **SC2** | `StrainCascade_LJA_assembly.sh` | La Jolla Assembler — multiplex de Bruijn graph assembly optimized for HiFi reads | Processed reads, sequencing type, threads, reproducibility mode |
| **SC3** | `StrainCascade_SPAdes_assembly.sh` | SPAdes — multi-*k*-mer de Bruijn graph assembly | Processed reads, sequencing type, threads, reproducibility mode |
| **SC4** | `StrainCascade_Canu_assembly.sh` | Canu — overlap-layout-consensus assembly with adaptive *k*-mer weighting | Processed reads, sequencing type, threads, reproducibility mode |
| **SC5** | `StrainCascade_Flye_assembly.sh` | Flye — repeat graph assembly for long error-prone reads | Processed reads, sequencing type, threads, reproducibility mode |
| **SC6** | `StrainCascade_Unicycler_assembly.sh` | Unicycler — hybrid (long+short) or long-read-only assembly (miniasm+Racon) | Processed reads, sequencing type, threads, reproducibility mode, **short reads R1/R2** (optional) |
| **SC7** | `StrainCascade_assembly_evaluation1.sh` | Assembly quality evaluation round 1 (pre-merge) | Assembly outputs from SC2–SC6, threads, selection algorithm |
| **SC8** | `StrainCascade_MAC2_assembly_merging.sh` | MAC2 consensus assembly merging from multiple assemblers | Assembly dir (output of SC2–SC7) |
| **SC9** | `StrainCascade_assembly_evaluation2.sh` | Assembly quality evaluation round 2 (post-merge) | Merged assembly from SC8, threads, selection algorithm |
| **SC10** | `StrainCascade_Circlator_circularisation.sh` | Contig circularization (Circlator) | Selected assembly, original reads, results integration dir, version |
| **SC11** | `StrainCascade_assembly_evaluation3.sh` | Assembly quality evaluation round 3 (post-circularization) | Circularized assembly from SC10, threads, results integration dir, version, selection algorithm |
| **SC12** | `StrainCascade_assembly_polishing.sh` | Assembly error correction (Arrow/Racon/Medaka + optional Polypolish with short reads) | Selected assembly, original reads, BAM file, sequencing type, threads, reproducibility mode, results integration dir, version, **short reads R1/R2** (optional) |
| **SC13** | `StrainCascade_minimap2_BBMap_coverage.sh` | Read mapping & coverage statistics (minimap2 + BBMap) | Polished assembly, original reads, sequencing type, threads |
| **SC14** | `StrainCascade_CheckM2_QC.sh` | Genome completeness & contamination QC (CheckM2) | Final assembly, threads, results integration dir, CheckM2 database, version |
| **SC15** | `StrainCascade_GTDB-Tk_taxonomy.sh` | Taxonomic classification (GTDB-Tk) | Final assembly, threads, taxonomic classification dir, results integration dir, GTDB database, version |
| **SC16** | `StrainCascade_GTDB-Tk_de_novo_tree.sh` | De novo phylogenetic tree construction (GTDB-Tk) | Final assembly, external assemblies dir, threads, taxonomic classification dir, GTDB database |
| **SC17** | `StrainCascade_Bakta_annotation.sh` | Genome annotation (Bakta v2) | Final assembly, taxonomy from SC15, threads, genome annotation dir, results integration dir, Bakta database, locus tag, force overwrite, version, reproducibility mode |
| **SC18** | `StrainCascade_Prokka_annotation.sh` | Prokaryotic genome annotation (Prokka) | Final assembly, threads, genome annotation dir, results integration dir, locus tag, force overwrite, version, reproducibility mode |
| **SC19** | `StrainCascade_DeepFRI_annotation.sh` | Deep-learning protein function prediction (DeepFRI) across GO categories + EC numbers | Protein sequences from SC17 (Bakta) or SC18 (Prokka) via genome annotation dir, results integration dir, version, reproducibility mode |
| **SC20** | `StrainCascade_MicrobeAnnotator_annotation.sh` | Metabolic pathway annotation (MicrobeAnnotator) | Annotation output from SC17/SC18 via genome annotation dir, threads, functional analysis dir, results integration dir, MicrobeAnnotator database, version, reproducibility mode |
| **SC21** | `StrainCascade_PlasmidFinder_identification.sh` | Plasmid replicon identification (PlasmidFinder) | Final assembly, results integration dir, PlasmidFinder database, version |
| **SC22** | `StrainCascade_AMRFinderPlus_antimicrobial_resistance_identification.sh` | Antimicrobial resistance gene detection (AMRFinderPlus) | Final assembly, taxonomy from SC15, threads, functional analysis dir, results integration dir, AMRFinderPlus database, version |
| **SC23** | `StrainCascade_ResFinder_antimicrobial_resistance_identification.sh` | Resistance gene identification (ResFinder) | Final assembly, threads, sequencing type, functional analysis dir, results integration dir, ResFinder database, version |
| **SC24** | `StrainCascade_dbCAN3_CAZymes_identification.sh` | Carbohydrate-active enzyme identification (dbCAN3) | Final assembly, threads, functional analysis dir, results integration dir, dbCAN database, version |
| **SC25** | `StrainCascade_IslandPath_genomic_islands_identification.sh` | Genomic island prediction (IslandPath-DIMOB) | Final assembly, annotation from SC17/SC18 (GBK file) via genome annotation dir, functional analysis dir, results integration dir, version |
| **SC26** | `StrainCascade_VirSorter2_phage_identification.sh` | Viral/phage sequence detection (VirSorter2, HMM multi-classifier) | Final assembly, threads, functional analysis dir, results integration dir, VirSorter2 database, version |
| **SC27** | `StrainCascade_geNomad_phage_identification.sh` | Viral & plasmid identification (geNomad, marker-based classification) | Final assembly, threads, functional analysis dir, results integration dir, geNomad database, version |
| **SC28** | `StrainCascade_CRISPRCasFinder_identification.sh` | CRISPR-Cas system detection (CRISPRCasFinder) | Final assembly, threads, functional analysis dir, results integration dir, version |
| **SC29** | `StrainCascade_ISEScan_IS_elements_identification.sh` | Insertion sequence element detection (ISEScan) | Final assembly, threads, functional analysis dir, results integration dir, version |
| **SC30** | `StrainCascade_data_integration.sh` | Consolidated result integration & interactive HTML report | Assembly dir, results integration dir, version, sample name |

**Notes:**

- **Shared inputs:** Every module receives `SCRIPT_DIR` and `LOGS_DIR` (injected by `execute_module()`), plus `LOG_NAME`, `UTILS_FILE`, and `APPTAINER_IMAGES_DIR` as the first three explicit parameters. Most modules also receive `OUTPUT_DIR` and `sample_name`, but not all (e.g., SC30 does not receive `OUTPUT_DIR`). `VERSION` is passed to SC10–SC12, SC14–SC15, SC17–SC30, but not to SC1–SC9, SC13, or SC16.
- **"Final assembly"** refers to the analysis assembly file found by `find_analysis_assembly_file()` within the genome assembly directory — the best assembly selected through the evaluation/merge/polish pipeline.
- **Assembly-input skip logic:** When the input type is `assembly` (user provides a pre-assembled genome via `-a`), the wrapper skips modules matching the regex `(Canu|LJA|SPAdes|Flye|Unicycler|evaluation[12]|MAC2|Circlator|assembly_polishing|BBMap)`. This skips **SC1–SC10, SC12, and SC13**. Notably, **SC11** (assembly_evaluation3) is **not** skipped because `evaluation[12]` only matches `evaluation1` and `evaluation2`, not `evaluation3`. SC14 and above always run.
- **Run summary:** `StrainCascade_run_summary.sh` is listed in `all_modules` but has no SC number — it runs automatically after all selected modules complete, outside the module selection loop.