# StrainCascade

## Overview
StrainCascade is a modular bioinformatics pipeline designed to process genomic data. For further information visit the [StrainCascade documentation page](https://sbujordi.github.io/StrainCascade/). 

## Dependencies
- git
- apptainer

## Installation
1. Clone the repository:
   ```bash
   # Clone the git repository
   git clone https://github.com/SBUJordi/StrainCascade.git

   # Change to the StrainCascade directory
   cd StrainCascade

   # Make the installation script executable
   chmod +x scripts/*
   
   # Execute the installation script (downloads databases & pulls docker images)
   ./scripts/StrainCascade_installation.sh
