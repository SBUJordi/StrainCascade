# StrainCascade

## Overview

StrainCascade is a modular bioinformatics pipeline designed to comprehensively process genomic data of bacterial isolates. Its automatic and customisable workflow includes everything from genome assembly, taxonomic identification, genome annotation, functional analysis, plasmid detection, as well as screens for antimicriobial resistance genes, CAZymes, phages and more. For further information on usage visit the [StrainCascade documentation page](https://sbujordi.github.io/StrainCascade/). Below you will find the minimum necessary installation information in case you already know your way around and do not need any further information.

## System requirements

-   Linux OS (e.g., CentOS, Ubuntu, Debian, etc.)

-   \~ 500GB disk space software and reference databases

## Installation dependencies

-   git

-   apptainer

## Installation steps

```         
# Clone the git repository
git clone https://github.com/SBUJordi/StrainCascade.git

# Change to the StrainCascade directory
cd StrainCascade

# Make the installation script executable
chmod +x scripts/*

# Execute the installation script (downloads databases & pulls docker images)
./scripts/StrainCascade_installation.sh
```
