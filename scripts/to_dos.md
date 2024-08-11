# To dos

## Additions

- IN DEVELOPMENT (Sebastian): Flye (genome assembler) (https://github.com/fenderglass/Flye)
- IN DEVELOPMENT (Sebastian): (Resistance Gene Identifier) module (https://github.com/arpcard/rgi)
- IN DEVELOPMENT (Sebastian): ResFinder (Resistence Gene Finder) module (https://github.com/cadms/resfinder)
- IN DEVELOPMENT (Sebastian): PlasmidFinder (https://github.com/kcri-tz/plasmidfinder)
- IN DEVELOPMENT (Sebastian): run_dbcan (Cazymes) (https://dbcan.readthedocs.io/en/latest/)


## Improvements
- Add another assembly_evaluation step directly after LJA to take the genome size from the LJA assembly and use it to overwrite the default settings (config file) for expected genome size for Canu, Flye and CISA
- Simplify structure of the results
    - rename pathway_analysis -> functional_analyis
    - move pathway results from microbeannotator to functional_analysis
- Improve selection for assembly: BIOMES_assembly_selection.py
- Test BIOMES_result_processing_for_R.r in April (genbankr package was taken from bioconductor but should go back in april - otherwise find replacement)
- Add script that retrieves all relevant system info and versions of used software. Should be written as an output
- Improve/expand the shiny app according to new read outs

## Set-up
- Docker
    - check how to handle miniconda envs in docker (or alternatives)
    - check if docker or singularity
    - request docker installation
