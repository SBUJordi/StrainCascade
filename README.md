# StrainCascade

## Overview
StrainCascade is a modular bioinformatics pipeline designed to process genomic data of bacterial isolates. For further information visit the [StrainCascade documentation page](https://sbujordi.github.io/StrainCascade/). 
Below you will find the minimum necessary installation information in case you already know your way around and do not need any further information. 

## Installation dependencies
- git
- apptainer

## Installation steps
   ```bash
   # Clone the git repository
   git clone https://github.com/SBUJordi/StrainCascade.git

   # Change to the StrainCascade directory
   cd StrainCascade

   # Make the installation script executable
   chmod +x scripts/*
   
   # Execute the installation script (downloads databases & pulls docker images)
   ./scripts/StrainCascade_installation.sh
