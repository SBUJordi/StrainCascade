![StrainCascade Logo](assets/straincascade_logo_colour_whitebackground.png)

## Overview

StrainCascade is a modular bioinformatics tool designed to comprehensively process genomic data of bacterial isolates, supporting both long-read sequencing data (PacBio or ONT) and pre-assembled genomes as input. Its automatic and customisable workflow includes everything from genome assembly, taxonomic identification, genome annotation, functional analysis, plasmid detection, as well as screens for antimicriobial resistance genes, CAZymes, phages and more. For further information on usage visit the [StrainCascade documentation page](https://sbujordi.github.io/StrainCascade_documentation/). Below you will find the minimum necessary installation information in case you already know your way around and do not need any further information.

## System requirements

-   Linux OS (e.g., CentOS, Ubuntu, Debian, etc.)

-   \~ 590GB disk space for software and reference databases

-   For optimal performance, we recommend running StrainCascade with 32 CPU cores, with 3GB of RAM allocated per core (= total of 96GB of RAM). If you run SC15 for *de novo* tree generation more than 100GB of RAM might be required.

## Installation dependencies

-   [Git](https://git-scm.com/)

-   [Apptainer](https://apptainer.org/)

-   ([Bash](https://tiswww.case.edu/php/chet/bash/bashtop.html) as command-line shell)

-   (Optional: [Conda/Miniconda](https://docs.anaconda.com/miniconda/))

-   All further software dependencies are bundled in Apptainer images, which are automatically downloaded and used when installing/running StrainCascade.

## Installation

For standard installation with the command-line interface, follow these steps:

Optional - Create and activate a conda environment for StrainCascade:

``` bash
conda create -n StrainCascade_env
conda activate StrainCascade_env
```

Optional - Navigate to your StrainCascade conda environment directory:

``` bash
cd $CONDA_PREFIX
```

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
./scripts/StrainCascade_installation.sh
```