![StrainCascade Logo](assets/straincascade_logo_colour_whitebackground.png)

## Overview

StrainCascade is a modular bioinformatics pipeline designed to comprehensively process genomic data of bacterial isolates, supporting both long-read sequencing data (PacBio or ONT) and pre-assembled genomes as input. Its automatic and customisable workflow spans genome assembly (5 assemblers with consensus merging), taxonomic classification, genome annotation (including deep-learning-based DeepFRI), functional analysis, plasmid detection, antimicrobial resistance screening, CAZyme identification, phage detection (VirSorter2 + geNomad), CRISPR-Cas detection, and more — across 30 analysis modules (SC1–SC30). Version 2.0.0 adds hybrid assembly support (Unicycler + Polypolish short-read polishing) and a deterministic reproducibility mode. For further information on usage visit the [StrainCascade documentation page](https://sbujordi.github.io/StrainCascade/). Below you will find the minimum necessary installation information in case you already know your way around and do not need any further information.

## System requirements

-   Linux OS (e.g., CentOS, Ubuntu, Debian, etc.)

-   \~ 1 TB disk space for software, Apptainer images, and reference databases

-   For optimal performance, we recommend running StrainCascade with 32 CPU cores, with 3GB of RAM allocated per core (= total of 96GB of RAM). If you run SC16 for *de novo* tree generation more than 100GB of RAM might be required.

## Installation dependencies

-   [Git](https://git-scm.com/)

-   [Apptainer](https://apptainer.org/)

-   ([Bash](https://tiswww.case.edu/php/chet/bash/bashtop.html) as command-line shell)

-   All further software dependencies are bundled in Apptainer images, which are automatically downloaded and used when installing/running StrainCascade.

## Installation

For standard installation with the command-line interface, follow these steps:

Clone the StrainCascade GitHub repository to your current working directory

``` bash
git clone https://github.com/SBUJordi/StrainCascade.git
```

Navigate to the new StrainCascade directory

``` bash
cd StrainCascade
```

Make all scripts (including the ones in subdirectories) of StrainCascade executable

``` bash
find scripts/ -type f -exec chmod +x {} \;
```

Execute the installation script, which pulls Apptainer images and downloads databases:

``` bash
# Interactive (recommended for first-time users)
./scripts/StrainCascade_installation.sh

# Non-interactive, idempotent (safe to re-run after a transient failure)
./scripts/StrainCascade_installation.sh --full --yes
```

For all installation options (custom installation, resuming a failed installation, verbosity, advanced multi-install setup with Conda) see the [installation documentation](https://sbujordi.github.io/StrainCascade/docs/installation.html).

## Citation

If you use StrainCascade in your research, please cite:

> Jordi SBU et al. StrainCascade: An automated, modular workflow for high-throughput long-read bacterial genome reconstruction and characterization. *iScience* 2026; 29:116189. https://doi.org/10.1016/j.isci.2026.116189