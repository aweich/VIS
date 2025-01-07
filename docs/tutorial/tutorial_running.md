## Running the pipeline
---
### Expected jobs 

```bash
> snakemake --cores 2 -n

...
Job stats:
job                                      count
-------------------------------------  -------
BAM_to_BED                                   2
Non_insertion_mapping                        2
all                                          1
basic_insertion_plots                        1
build_insertion_reference                    1
calc_distance_to_elements                    2
calculate_exact_insertion_coordinates        2
clean_postcut_by_maping_quality              2
collect_outputs                              2
copy_config_version                          1
detailed_fragmentation_length_plot           2
extract_by_length                            2
extract_fastq_insertions                     2
extract_mapping_quality                      2
finalize_mapping_quality                     2
find_insertion_BLASTn                        2
find_insertion_BLASTn_in_Ref                 2
fragmentation_distribution_plots             2
generate_mapq_heatmap                        2
get_coordinates_for_fasta                    2
hardcode_blast_header                        2
insertion_fragmentation                      1
insertion_mapping                            2
make_blastn_DB                               1
make_fasta_without_tags                      2
minimap_index                                1
multiqc                                      1
nanoplot                                     2
prepare_insertion                            1
read_level_fastqc                            2
sort_insertion_file                          2
split_fasta                                  2
total                                       55

```

All of these jobs will be executed in the correct order by the workflow. So let's finally run it.

<br>

### Execution
```bash
    snakemake --cores 2
```
!!! info 
    
    Depending on the number of cores specified and whether the environments need to be built for the first time, this process may take a while. However, since the simulated data is very small, the expected runtime should not exceed 5-10 minutes.

If you see this message, the workflow has executed successfully and completed. If not, refer to the [error handling](../advanced_usage/advanced_usage_errors.md) section. 

```bash
    Finished job 0.
    55 of 55 steps (100%) done
``` 

Now, run snakemake again with its built-in `--report` functionality to get a comprehensive overview of the workflow's runtime and output. 

```bash
    snakemake --report
```
<br>