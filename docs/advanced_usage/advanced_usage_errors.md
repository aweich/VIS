### Error handling

#### Snakemake
General debugging ressources for everything related to snakemake can be found in the snakemake [FAQ](https://snakemake.readthedocs.io/en/v6.15.5/project_info/faq.html).

#### Log files
The pipeline is designed with rule-specific `log` files, which are stored in the `intermediate` output directory. These logs serve as the primary resource for identifying and addressing any rule-specific issues that arise during execution. If you encounter errors or unexpected behavior, these files should be your first point of reference for debugging. 
