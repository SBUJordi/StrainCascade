# BIOMES_WGseq

BIOMES_WGseq is a modular pipeline to analyse whole genome sequencing data of bacterial isolates. The following tasks are performed 
- Assembly (different assemblers)
- Quality assessment of assembly
- Selection of best assembly
- Taxonomic classification
- Gene annotation (different tools)

The pipeline can be highly parallelised and handle as many samples as desired at once. Hoewver, each sample will be submitted as an individual task needing 32 CPUs with 10GB memory per CPU (320GB total). Thus, the system's resources can become a constraint. The 32 CPUs are used to parallelise the processes of many of the modules.

The processing time for a PacBio HiFi sequencing file from a bacterial isolate (5Mbps), including execution of all modules, is approximately 2.5h.


## Pipeline Structure

- **`submit_BIOMES_WG_seq.sh`** (submitter script to initiate the pipeline)
    - **[`config.sh`](#configurations-in-configsh)** (configuration script)
    - **`utils.sh`** (defines re-occurring functions)
    - **[`BIOMES_WG_seq_pipeline.sh`](#biomes-whole-genome-sequencing-wgseq-analysis-pipeline-overview)** (wrapper script initiating and orchestrating the pipeline modules)
        1. **[`BIOMES_LJA_assembly.sh`](#la-jolla-assembler-lja-module)** (La Jolla Assembler)
        1. **[`BIOMES_SPAdes_assembly.sh`](#st-petersburg-genome-assembler-spades-module)** (SPAdes Assembler)
        1. **[`BIOMES_Canu_assembly.sh`](#canu-assembly-module)** (Canu Assembler)
        2. **[`BIOMES_CISA_merging.sh`](#cisa-contig-integrator-for-sequence-assembly-assembly-merging-module)** (Contig Integrator for Sequence Assembly; merging of assemblies)
        3. **[`BIOMES_assembly_evaluation.sh`](#assembly-evaluation-quast-and-selection-module)** (QUality ASsessment Tool [QUAST] and BIOMES_assembly_selection.py; evaluation and selection fo best assembly)
            - **[`BIOMES_assembly_selection.py`](#biomes_assembly_selectionpy)** (home-made assembly selection algorithm)
        4. **[`BIOMES_Circlator_circularisation.sh`](#circlator-circularization-module)** (Circlator; assembly circularisation)
        5. **[`BIOMES_assembly_evaluation2.sh`](#assembly-evaluation-quast-and-selection-module)** (QUality ASsessment Tool [QUAST] only; quality control of the circularised assembly)
            1. **[`BIOMES_gtdbtk_taxonomy.sh`](#gtdb-tk-taxonomic-classification-module)** (Genome Taxonomy Database [GTDB] Toolkit [Tk]; taxonomic classification)
                - **`BIOMES_extract_organism_name_gtdbtk.py`** (script to extract key taxonomic info)
                    1. **[`BIOMES_prokka_annotation.sh`](#prokka-genome-annotation-module)** (Prokka; genome annotation)
                    1. **[`BIOMES_Bakta_annotation.sh`](#bakta-genome-annotation-module)** (Bakta; genome annotation)
                    1. **[`BIOMES_MicrobeAnnotator_annotation.sh`](#microbeannotator-genome-annotation-and-pathway-analysis-modules)** (MicrobeAnnotator; genome annotation and pathway analysis)
                    
                        1. **[`BIOMES_results_summary.sh`](#results-summary-module)** (generates a concise summary)
                        1. **[`BIOMES_result_processing.sh`](#results-summary-module)** (processes results and creates an Rproject for downstream analysis in R)
                            - **[`BIOMES_result_processing_for_R.r`](#results-summary-module)** (generates an .RData file with the results)
                            - **[`BIOMES_WGseq_shiny_app.r`](#results-summary-module)** (code to run a shiny app with the provided .RData file)

                        
                        
            
    
## General Instructions 

### Necessary Directory Structure for Run

- `run_directory` (= parent dir; can have any name), 
    - `input_files` (1. order child directory; needs to be named "input_files"),
        - `read_file1.fasta`
        - `read_file2.fasta`
        - `read_file3.fasta`
        - ...
    - `script_directory` (1. order child directory; can have any name)
        - all pipeline scripts
        - no subfolders


### Necessary Manual Configurations

- Adaptions in the `config.sh`
    - dd
- Adaptions in the `submit_BIOMES_WG_seq.sh`
    - dd


### Running the Job

- Before you start you need to place all sequencing files (.fasta) in the `input_files` directory. Chose file name carefully as they will be used for naming the results
- To start, you need to launch the job as a SLURM job (`sbatch submit_BIOMES_WG_seq.sh`) from the `script_directory` (script_directory = working directory). 

    > Explanation: The `SLURM_SUBMIT_DIR` (used in many of the scripts) environment variable is set to the directory from which the SLURM job was submitted, not the location of the script itself. This means it's the directory you were in when you ran the `sbatch` command to submit the job.


## Instructions for our group @unibe

Follow general instructions. Before submitting (`sbatch submit_BIOMES_WG_seq.sh`) activate your workspace (e.g., `HPC_WORKSPACE=nbmisi module load Workspace` for the NBMISI workspace). 


## Detailed instructions

### [Configurations in `config.sh`](#configurations-in-configsh)

You need to modify certain variables in the `config.sh` file:

- `only_main_results`: Set this to "yes" if you only want the key results. Otherwise, all results will be kept.
- `miniconda_path`: Set this to the path where miniconda is installed.
- `bakta_db_path`: Set this to the path where the Bakta database is located.
- `read_type`: Set this to the type of reads for SPAdes.
- `genome_size_Canu`: Set this to an estimate of the genome size for Canu.
- `genome_size_CISA`: Set this to an estimate of the genome size for CISA.
- `sequencing_type`: Set this to the sequencing type for Canu.
- `selected_modules`: Uncomment the modules you want to include in the pipeline.


### [BIOMES Whole Genome Sequencing (WGseq) Analysis Pipeline Overview](#biomes-whole-genome-sequencing-wgseq-analysis-pipeline-overview)

The pipeline is divided into two main scripts: 

- **Submitter Script (`submit_BIOMES_WG_seq.sh`):**
  - Manages the batch processing of input files.
  - Identifies all `.fastq.gz` files in the `input_files` directory.
  - Submits them in batches of 10 to the SLURM job scheduler for processing by the wrapper script.
  - If there are remaining files, an additional job is submitted to handle them.
  - Sets the `--cpus-per-task` argument which is important for the parallelisation of many of the modules; default is 32 CPUs per task (sample). 
  - Each sample is run as an individual task sequentially through the modules.

- **Wrapper Script (`BIOMES_WG_seq_pipeline.sh`):**
  - Handles the actual processing of each input file.
  - Loads the configuration and utility scripts, creates directories and defines variables.
  - For each input file, it creates an output directory and a series of subdirectories for storing the results of each processing step.

- **Modularity:**
  - The pipeline is modular, allowing for the execution of selected modules based on the configuration.
  - New modules can easily be incorporated.
  - Modules include various assembly methods (La Jolla Assembler, SPAdes, Canu), assembly merging (CISA), assembly evaluation (home-made), circularisation (Circlator), taxonomic classification (GTDB-Tk), and genome annotation (Prokka, Bakta).
  - Each module is executed with its specific parameters, and the results are stored in the appropriate subdirectory.
  - If the `only_main_results` configuration is set to "yes", only the main results are kept and the rest is deleted.

- **Flexibility:**
  - This pipeline structure allows for a flexible, customizable, and efficient processing of whole genome sequencing data.


### [La Jolla Assembler (LJA) Module](#la-jolla-assembler-lja-module)

The LJA module is implemented in the `BIOMES_LJA_assembly.sh` script.


#### Conda Environment

The script activates the *LJA* conda environment where the `lja` function is installed. The environment needs to be called "LJA".


#### Running LJA

The script runs the La Jolla Assembler utilising 32 CPU threads.


### [St. Petersburg genome Assembler (SPAdes) Module](#st-petersburg-genome-assembler-spades-module)

The SPAdes assembly module is implemented in the `BIOMES_SPAdes_assembly.sh` script.


#### Conda Environment

The script activates the *SPAdes* conda environment where the `spades.py` function is installed. The environment needs to be called "SPAdes".


#### Running SPAdes

The script runs SPAdes utilising 32 CPU threads. The specific SPAdes command used depends on the `read_type` argument, which can be either "PacBio_CLR" or "single".


### [Canu Assembly Module](#canu-assembly-module)

The Canu assembly module is implemented in the `BIOMES_Canu_assembly.sh` script.


#### Conda Environment

The script activates the *Canu* conda environment where the `canu` function is installed. The environment needs to be called "Canu".


#### Running Canu

The script runs Canu on the input file with specified genome size, utilising 32 CPU threads with 10GB memory each.


### [CISA (Contig Integrator for Sequence Assembly) Assembly Merging Module](#cisa-contig-integrator-for-sequence-assembly-assembly-merging-module)

The BIOMES_CISA_merging module is implemented in the `BIOMES_CISA_merging.sh` script.


#### Input files

The script uses all previously created assembly files.


#### Conda Environment

The script activates the *cisa* conda environment where the `Merge.py` and `CISA.py` functions are installed. The environment needs to be called "cisa".


#### Running CISA

If there is more than one assembly file, the script will use all of them and merge them using CISA.


### [Assembly Evaluation (QUAST) and Selection Module](#assembly-evaluation-quast-and-selection-module)

#### BIOMES_assembly_evaluation(2).sh

This bash script is used for evaluating genome assemblies. It processes and compares all assemblies.

1. Finds all assembly files including the merged assembly from CISA.
2. The script activates the *quast* conda environment where the `quast.py` function is installed. The environment needs to be called "quast".
3. Runs QUAST iteratively on the copied assemblies.
4. Creates an assembly evaluation file with the key quality indices for all assemblies.
5. Deactivates the quast environment.
6. Activates the *default_python_env* environment where the `python` and python's `pandas` are installed. The environment needs to be called "default_python_env".
7. Runs the `BIOMES_assembly_selection.py` script for assembly selection (not for `BIOMES_assembly_evaluation2.sh`).
8. Deactivates the default Python environment.


#### [BIOMES_assembly_selection.py](#biomes_assembly_selectionpy)

This Python script is used to select the best assembly from the available set of assemblies. It takes as input the TSV file containing assembly statistics generated with QUAST. The script performs the following steps:


##### Contig Size Ranges and Ratios

The script calculates ratios for each assembly across different contig size ranges, expressing the contribution of each contig size range to the total assembly length in base pairs (bp). The following contig size ranges are used.

1. **Total assembly length (contigs >= 0 bp):** Overall assessment of assembly length.
2. **Assembly length (bp) contributed by contigs >= 1000 bp:** Excludes very short contigs.
3. **Assembly length (bp) contributed by contigs >= 5000 bp:** Prioritizes assemblies with a higher fraction of somewhat longer contigs.
4. **Assembly length (bp) contributed by contigs >= 10'000 bp:** Prioritizes assemblies with a higher fraction of even longer contigs.
5. **Assembly length (bp) contributed by contigs >= 25'000 bp:** Prioritizes assemblies with a higher fraction of relatively long contigs.
6. **Assembly length (bp) contributed by contigs >= 50'000 bp:** Prioritizes assemblies with a higher fraction of very long contigs, suggesting high-quality genomes.

Ratios are calculated by dividing the length contributions of each contig size range by the total length of all contigs (assembly length contributed by contigs >= 0 bp). These ratios facilitate quantitative comparisons, aiding the selection of assemblies based on both overall length and the distribution of contig lengths in specific size categories.


##### Selection Process

The selection of assemblies is based on the following process:

1. **Sanity Check - Genome Size:**
   - Excludes assemblies from selection process if their total lenths is >10Mbps (Fig. 1a https://doi.org/10.1038/s41467-023-37396-x).

2. **Sanity Check - Genome Size:**
   - If there are more than two assemblies left, it calculates the median of their total length. Assemblies that deviate >=1.96 standard deviations from the median are excluded.

3. **Ratio Calculation:**
   - Ratios are calculated for each assembly by dividing the contribution of different contig size range by the total length of the assemblies (length of contigs >= 0 bp).

4. **Selection Criteria:**
   - The assembly with the highest ratio is considered the most favorable for within each contig size range.

5. **Iteration through Size Ranges:**
   - The script iterates through these ratios starting with the largest size range (contribution of contigs >= 50'000 bp to total assembly length).
   - At each step the script aims to find a single winning assembly based on the highest ratio value.
   - If there is  one single assembly with the highest ratio among all assemblies at a step, that assembly is immediately chosen as the final selected assembly, and the process stops.

6. **Priority Order of Size Ranges:**
   - The process of comparing ratios goes in descending order of the contig size range. 
   - This determines the priority, starting with the largest contig size range as most relevant and proceeding to smaller size ranges.

7. **Additional Criteria if Needed:**
   - If, after considering all contig size ranges, there is still no single winning assembly, the script applies additional criteria (number of contigs, largest contig size, and total length) to break ties.

8. **Random Selection if Still No Single Winner:**
   - If, after applying all criteria, there is still no single winning assembly, the script randomly chooses one from the remaining assemblies.


##### Summary

In summary, the script systematically selects a single assembly based on the highest ratio in specific contig size ranges. **Ratios are calculated by expressing each contig size range as a proportion of the total contig length**. Thereby **prefering assemblies with higher contributions of large contigs** to the total assembly length. 


### [Circlator Circularization Module](#circlator-circularization-module)

The BIOMES_Circlator_circularisation module is implemented in the `BIOMES_Circlator_circularisation.sh` script.


#### Input files

The script uses the previously chosen assembly file (compare section *Assembly Evaluation (QUAST) and Selection Module*) as well as the original sequencing reads file.


#### Conda Environment

The script activates the *Circlator* conda environment where the `circlator` function is installed. The environment needs to be called "Circlator".


#### Running Circlator

Circlator circularises the chosen assembly file. It needs the original sequencing reads as additional input. Circlator is implemented that it depends on SPAdes; Canu could be used alternatively.


### [GTDB-Tk Taxonomic Classification Module](#gtdb-tk-taxonomic-classification-module)

The BIOMES_gtdbtk_taxonomy module is implemented in the `BIOMES_gtdbtk_taxonomy.sh` script.


#### Input file

The script uses the previously chosen and circularised assembly file (compare sections *Assembly Evaluation (QUAST) and Selection Module* and *Circlator Circularization Module*).


#### Conda Environment 1

The script activates the *gtdbtk-2.3.2* conda environment where the `classify_wf` function is installed. The environment needs to be called "gtdbtk-2.3.2" even if another version is installed; alternatively the module script needs to be adapted.


#### Running GTDB-Tk

GTDB-Tk classifies the assembly using the `classify_wf` function. The function utilises 32 CPU threads.


#### Conda Environment 2

The script activates the *default_python_env* environment where the `python` and python's `pandas` are installed. The environment needs to be called "default_python_env".


#### Running Taxonomic Information Extraction

The the Python script `BIOMES_extract_organism_name_gtdbtk.py` extracts taxonomic information from the GTDB-Tk results. These information are used by downstream modules (e.g. the *Bakta Genome Annotation Module*).


### [Prokka Genome Annotation Module](#prokka-genome-annotation-module)

The BIOMES_prokka_annotation module is implemented in the `BIOMES_prokka_annotation.sh` script.


#### Input files

The script uses the previously chosen and circularised assembly file (compare sections *Assembly Evaluation (QUAST) and Selection Module* and *Circlator Circularization Module*).


#### Conda Environment

The script activates the *prokka* environment where the `prokka` function is installed. The environment needs to be called "prokka".


#### Running Prokka

Prokka performs genome annotation on the chosen and circularised assembly file. The annotation process utilises 32 CPU threads.


### [Bakta Genome Annotation Module](#bakta-genome-annotation-module)

The BIOMES_Bakta_annotation module is implemented in the `BIOMES_Bakta_annotation.sh` script.


#### Input files

The script uses the previously chosen and circularised assembly file (compare sections *Assembly Evaluation (QUAST) and Selection Module* and *Circlator Circularization Module*).


#### Conda Environment

The script activates the *Bakta* environment where the `bakta` function is installed. The environment needs to be called "Bakta".


#### Running Bakta

Bakta performs genome annotation on the chosen and circularised assembly file. The annotation process utilises 32 CPU threads and the extracted taxonomic info from the *GTDB-Tk Taxonomic Classification Module* if available.


### [MicrobeAnnotator Genome Annotation and Pathway Analysis Modules](#microbeannotator_genome-annotation-and-pathway-analysis-module)


### [BIOMES Results Summary Module](#results-summary-module)

The module to summarise the assembly, classification and genome annotation is implemented in the `BIOMES_results_summary.sh` script.


#### Conda Environment

The script activates the *default_python_env* conda environment where the `python` and python's `pandas` are installed. The environment needs to be called "default_python_env".


#### Running the Summary Script

The script executes a Python script that generates a concise summary of BIOMES_WGseq results, covering aspects such as assembly size, contig count, taxonomic identification, annotated gene count, etc.


## Debugging Comments

### Circlator

#### Installation

**Error** (during installation): 
```
Found canu but couldn't get version.
```

**Possible solution**:
Find the script "external_progs.py" and change "Canu" to "canu" in line 22.


#### Problem with nucmer version

**Error** (while running): 
```
circlator.external_progs.Error: Version of nucmer too low. I found ., but must be at least 3.1.
```

**Possible solution**: 
If you have perl warnings when running the prevent circlator to recognise the nucmer version (string-based recognision). Example of warnings:
```
nucmer --version
perl: warning: Setting locale failed.
perl: warning: Please check that your locale settings:
	LANGUAGE = (unset),
	LC_ALL = (unset),
	LC_CTYPE = "UTF-8",
	LANG = "en_US.UTF-8"
    are supported and installed on your system.
perl: warning: Falling back to a fallback locale ("en_US.UTF-8").
nucmer
NUCmer (NUCleotide MUMmer) version 3.1
```

To solve this you can add `export LC_ALL=C` to your bash profile file:

```
nano ~/.bashrc
```

add `export LC_ALL=C`


```
source ~/.bashrc
```

### MicrobeAnnotator
#### Database download
The pipeline was tested with the full Database. The database needs to be called MicrobeAnnotator_DB and should be located in the ${miniconda_path}/envs/microbeannotator/MicrobeAnnotator_DB location

### Example of Working Bash Profile

```
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/storage/homefs/sj21o643/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/storage/homefs/sj21o643/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/storage/homefs/sj21o643/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/storage/homefs/sj21o643/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# Add the following line to include SPAdes in the PATH
export PATH="/storage/homefs/sj21o643/miniconda3/envs/SPAdes/SPAdes-3.15.5-Linux/bin:$PATH"

# Add the following line to include Canu in the PATH
export PATH=/storage/homefs/sj21o643/miniconda3/envs/Canu/canu-2.2/bin:$PATH

export PATH=$PATH:/storage/homefs/sj21o643/miniconda3/envs/cisa/CISA1.3

# The following (export LC_ALL=C) is to  prevent perl warnings that interfered with the nucmer version recognision by circlator
export LC_ALL=C
```


## Authorship and Ownership Statement

This pipeline, BIOMES_WGseq, is a proprietary tool developed by Sebastian Bruno Ulrich Jordi. The pipeline is not open-source but might circulate within specific groups for collaborative purposes. Sebastian Bruno Ulrich Jordi retains full ownership and intellectual property rights to the BIOMES_WGseq pipeline.


### Author
- Sebastian Bruno Ulrich JORDI


### Ownership
The BIOMES_WGseq pipeline is owned and maintained by Sebastian Bruno Ulrich JORDI. Any distribution, modification, or use of this pipeline beyond the authorized group should be with explicit permission from the author only.


### Disclaimer
This README file serves as documentation for understanding and using the BIOMES_WGseq pipeline. Any unauthorized distribution of the BIOMES_WGseq pipeline, parts of it, modifications, or use is not permitted without the explicit consent of Sebastian Bruno Ulrich Jordi.


### Third-Party Software
The BIOMES_WGseq pipeline utilises third-party software, including but not limited to La Jolla Assembler (LJA), St. Petersburg genome Assembler (SPAdes), Canu, Contig Integrator for Sequence Assembly (CISA), QUality ASsessment Tool (QUAST), Circlator, Genome Taxonomy Database Toolkit (GTDB-Tk), Prokka, and Bakta. Users are required to adhere to the respective licenses of these tools and give appropriate credit as specified in the individual licenses.


### Contact Information
For inquiries or collaboration requests, please contact the author:
- Sebastian Bruno Ulrich JORDI: [sbuj_research.87w0y@passfwd.com](mailto:sbuj_research.87w0y@passfwd.com)

