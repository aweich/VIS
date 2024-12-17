# Tutorial

Let's now look at a minimal example of what the pipeline can do for you. For this, we have simulated some sequencing data and randomly added insertions to some of them. For more on this, check [here]().

Before we are getting started, make sure that you have all at hand that is needed. You can easily check this by typing: 

```bash
    snakemake -n
```
If you now see a list of jobs that wait for execution, you are fully equipped for what comes next.

!!! info For this tutorial, we execute the workflow based on the alternative setup. Just add `--use-conda` if this is not what you chose. 

## Collect the reference files

To let the pipeline know where our samples and other dependencies are located, we now need to open and edit the configuration file. Open the `config.yml` file and fill the missing dependencies as follows:

```yaml
# tutorial config
experiment: "tutorial"
samples:
    S1: "tutorial/simulated/S1.bam",
    S2: "tutorial/simulated/S2.bam"
processing_dir: "tutorial/out"
vector_fasta: "tutorial/references/vectorseq.fa"
splitmode: "Buffer"
fragment_size: 100
MinLength: 1
MAPQ: 10
MinInsertionLength: 500
ref_genome_ctrl: "tutorial/references/chr1region.fa"
annotation_1: "tutorial/references/UCSC_GENCODEV44_chr1region.bed"
detection: "rules/detection.smk"
quality_control: "rules/qc.smk"
functional_genomics: "rules/functional_genomics.smk"
```

## Run the pipeline

Let's run one more dry-run of the workflow to see if the `config` has been put together correctly. 

```bash
> snakemake -n

job                                      count
-------------------------------------  -------
BAM_to_BED                                   1
Non_insertion_mapping                        1
all                                          1
basic_insertion_plots                        1
blast_to_gff                                 1
build_insertion_reference                    1
calc_distance_to_elements                    1
calculate_exact_insertion_coordinates        1
clean_postcut_by_maping_quality              1
collect_outputs                              1
copy_config_version                          1
detailed_fragmentation_length_plot           1
extract_by_length                            1
extract_fastq_insertions                     1
extract_mapping_quality                      1
finalize_mapping_quality                     1
find_vector_BLASTn                           1
find_vector_BLASTn_in_humanRef               1
fragmentation_distribution_plots             1
generate_mapq_heatmap                        1
get_cleavage_sites_for_fasta                 1
hardcode_blast_header                        1
insertion_mapping                            1
make_blastn_DB                               1
make_fasta_without_tags                      1
minimap_index                                1
multiqc                                      1
nanoplot                                     1
plot_distance_to_elements                    1
plot_scoring                                 1
prepare_vector                               1
read_level_fastqc                            1
sort_insertion_file                          1
split_fasta                                  1
vector_fragmentation                         1
total                                       35
```

All of these jobs will be executed by the workflow. So let's finally run it.

```bash
    snakemake --cores 2
```
!!! info Depending on the amount of cores specified, this might take a while. However, since the simulated data is reasonably small, the expected runtime should not exceed more than 5 minutes.  

If you see this, the workflow executed successfully and ran through completely. If not, jump to the [error handling](#error-handling) section. 

```bash
    Finished job 0.
    35 of 35 steps (100%) done
``` 

Now run snakemake again with its inbuilt `--report` functionality to get a comprehensive overview about the workflow's runtime and output. 

```bash
    snakemake --report
```

### Inspecting the output

#### Report

Since we now have also automatically generated a general report for the workflow stored in the working directory of the pipeline. 
Take a look at the Statistics in `report.html`. Some rules took significantly longer to finish than others.

<img src="images/report_statistics.png" alt="Report_Statistics" width="400">

Throughout the pipeline, some simple plots are also generated to give us a glimpse of what the insertions look like in terms of their length or chromosomal specificity. Navigate to the results tab and take a look at the detected lengths of your insertions. It looks like some of the reads only contained parts of the insertion.

<img src="images/Insertion_length.png" alt="images/Insertion_length.png" width="800">

If you further want to get an idea about quality control metrics, navigate to the `multiqc.html` report in the results tab. 

#### Output files

Let's now have a look at the directly generated output files. Navigate to the output folder as defined in the `config`. If you want to have an overview about the file structure of this directory, run `tree "tutorial/out/tutorial"`. An overview can also be found [here](other.md/#output-directory-structure). 

####

To do:
# Run two sample simulated data with in total less than 25mb?
# Add their output to the tutorial documentation
# Add their data to the github
# finish the documentation

Next to this overview, there are also other files generated that can subsequently be used for further investigations. The directory structure provided will pre-sort the outputs into intermediate and final output. Navigate to the final output/final/localization first.

Here, the most interesting one, of course, is the sample-specific BED file in containing the genomic positions of the insertions. This is also used for the annotation of each insertion (see functional genomics).

In output/final/qc, the subfolder fragmentation contains detailed information about each sample's fragmented insertion coverage, i.e. which parts of the insertion sequence was detected in the reads, what are the consecutive intervals for each read's insertion, and if there were also parts of the insertions detected, that match the human reference genome. In summary, this output gives you an idea about the effectiveness of the fragmentation of your inserted sequence for the detection of the read. Depending on the size of the insertion, you might want to run the analysis pipeline again and change the fragment-length in the config to further increase your level of detail.   

In the output/intermediate folder, you can find various subdirectories with intermediate files created during the pipeline run. Most of them are self-explanatory when looking at the workflow of the pipeline. To make things easier, an illustration about which files are created were in the DAG, please see below.  

Congratulations, you have finished the quick start using simulated data! If you want to unleash the full power of the pipeline, feel free to continue with the advanced usage tutorial below.


### functional genomics


## Error handling

You can now follow the different rules in your terminal window. If you encounter errors, make sure to double-check your initial input. If your error is for a specific rule, check the detailed documentation for this step in ./log/rulename

Other general debugging ressources for everything related to snakemake can be found [here] () or [here] ().