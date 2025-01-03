# Tutorial


Now, let's explore a minimal example of what the pipeline can accomplish. To demonstrate this, we’ve simulated some sequencing data and randomly introduced insertions in some of the samples. For more details on the data, refer to [this section](other.md/#simulate_data_for_tutorial).

## Before running the pipeline
---
### Prepare config.yml

To inform the pipeline about the location of our samples and other dependencies, we need to open and edit the configuration file. Open the `config.yml` file and fill in the missing dependencies as outlined below:

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
<br>

### Check setup

Before we begin, make sure you have all the necessary files ready. You can quickly verify this by running: 

```bash
    snakemake -n
```

If you see a list of jobs waiting for [execution](#expected-jobs), you're all set for the next steps.

!!! info 
    
    For this tutorial, we will run the workflow using the alternative setup (i.e., one virtual environment). If this isn't the setup you've chosen, simply add `--use-conda` to the snakemake commands. 

<br>

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

If you see this message, the workflow has executed successfully and completed. If not, refer to the [error handling](#error-handling) section. 

```bash
    Finished job 0.
    55 of 55 steps (100%) done
``` 

Now, run snakemake again with its built-in `--report` functionality to get a comprehensive overview of the workflow's runtime and output. 

```bash
    snakemake --report
```
<br>

## After running the pipeline
---
### Snakemake report

We have also automatically generated a general report for the workflow, which is stored in the working directory of the pipeline. Take a look at the statistics in `report.html`. Some rules took way longer to complete than others.

![report_statistics.png](images/tutorial/report_statistics.png)

Throughout the pipeline, several simple plots are generated to give insights into the insertions' characteristics, such as their length and chromosomal specificity. Navigate to the results tab to explore the detected insertion lengths. It appears that some reads only contain parts of the insertion.

![Insertion_Length.png](images/tutorial/Insertion_length.png)

If you would like to explore quality control metrics, check out the `multiqc.html` report in the results tab. 

<br>

### Output directory structure

Now, let's examine the output files directly generated by the pipeline. Navigate to the output folder as specified in the `config`. To get an overview of the file structure in this directory, run `tree "tutorial/out/tutorial"`. 

<details><summary> Output directory structure </summary>

```plaintext

├── config_settings.yml
├── final
│   ├── functional_genomics
│   │   ├── Functional_distances_to_Insertions_S1.bed
│   │   └── Functional_distances_to_Insertions_S2.bed
│   ├── localization
│   │   ├── ExactInsertions_S1.bed
│   │   ├── ExactInsertions_S2.bed
│   │   ├── Heatmap_Insertion_Chr.png
│   │   └── Insertion_length.png
│   └── qc
│       ├── Fragmentation
│       │   ├── Insertions
│       │   │   ├── insertions_100_S1
│       │   │   │   ├── 100_fragmentation_distribution.png
│       │   │   │   └── 100_read_match_fragmentation_distribution.png
│       │   │   └── insertions_100_S2
│       │   │       ├── 100_fragmentation_distribution.png
│       │   │       └── 100_read_match_fragmentation_distribution.png
│       │   ├── Longest_Interval
│       │   │   ├── S1
│       │   │   │   ├── Longest_interval_Read-343.png
│       │   │   │   ├── Longest_interval_Read-555.png
│       │   │   │   ├── Longest_interval_Read-561.png
│       │   │   │   ├── Longest_interval_Read-745.png
│       │   │   │   └── Longest_interval_Read-902.png
│       │   │   └── S2
│       │   │       ├── Longest_interval_Read-262.png
│       │   │       ├── Longest_interval_Read-417.png
│       │   │       ├── Longest_interval_Read-522.png
│       │   │       ├── Longest_interval_Read-682.png
│       │   │       └── Longest_interval_Read-824.png
│       │   └── Reference
│       │       ├── reference_100_S1
│       │       │   └── 100_fragmentation_distribution.png
│       │       └── reference_100_S2
│       │           └── 100_fragmentation_distribution.png
│       ├── mapq
│       │   ├── Insertions_S1_mapq.txt
│       │   ├── Insertions_S2_mapq.txt
│       │   ├── S1_mapq_heatmap_image.png
│       │   └── S2_mapq_heatmap_image.png
│       └── multiqc_report.html
└── intermediate
    ├── blastn
    │   ├── Coordinates_100_InsertionMatches_S1.blastn
    │   ├── Coordinates_100_InsertionMatches_S2.blastn
    │   ├── Filtered_Annotated_100_InsertionMatches_S1.blastn
    │   ├── Filtered_Annotated_100_InsertionMatches_S2.blastn
    │   ├── Readnames_100_InsertionMatches_S1.txt
    │   ├── Readnames_100_InsertionMatches_S2.txt
    │   └── ref
    │       ├── Filtered_Annotated_100_InsertionMatches_S1.blastn
    │       └── Filtered_Annotated_100_InsertionMatches_S2.blastn
    ├── fasta
    │   ├── fragments
    │   │   ├── 100_Insertion_fragments.fa
    │   │   ├── 100_Insertion_fragments.fa.ndb
    │   │   ├── 100_Insertion_fragments.fa.nhr
    │   │   ├── 100_Insertion_fragments.fa.nin
    │   │   ├── 100_Insertion_fragments.fa.not
    │   │   ├── 100_Insertion_fragments.fa.nsq
    │   │   ├── 100_Insertion_fragments.fa.ntf
    │   │   ├── 100_Insertion_fragments.fa.nto
    │   │   └── Forward_Backward_Insertion.fa
    │   ├── Full_S1.fa
    │   ├── Full_S2.fa
    │   ├── Insertion_S1.fa
    │   ├── Insertion_S2.fa
    │   ├── Modified_S1.fa
    │   └── Modified_S2.fa
    ├── localization
    │   ├── ExactInsertions_S1.bed
    │   └── ExactInsertions_S2.bed
    ├── log
    │   ├── detection
    │   │   ├── BAM_to_BED
    │   │   │   ├── Postcut_S1.log
    │   │   │   ├── Postcut_S2.log
    │   │   │   ├── Precut_S1.log
    │   │   │   └── Precut_S2.log
    │   │   ├── basic_insertion_plots
    │   │   │   ├── heat.log
    │   │   │   └── length.log
    │   │   ├── build_insertion_reference
    │   │   │   └── out.log
    │   │   ├── calculate_exact_insertion_coordinates
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── clean_postcut_by_maping_quality
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── collect_outputs
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── copy_config_version
    │   │   │   └── out.log
    │   │   ├── extract_by_length
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── find_insertion_BLASTn
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── find_insertion_BLASTn_in_Ref
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── get_coordinates_for_fasta
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── hardcode_blast_header
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── insertion_fragmentation
    │   │   │   └── out.log
    │   │   ├── insertion_mapping
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── make_blastn_DB
    │   │   │   └── out.log
    │   │   ├── make_fasta_without_tags
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── minimap_index
    │   │   │   └── out.log
    │   │   ├── Non_insertion_mapping
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   ├── prepare_insertion
    │   │   │   └── out.log
    │   │   └── split_fasta_by_borders
    │   │       ├── S1.log
    │   │       └── S2.log
    │   ├── functional_genomics
    │   │   ├── calc_distance_to_elements
    │   │   │   ├── S1.log
    │   │   │   └── S2.log
    │   │   └── sort_insertion_file
    │   │       ├── S1.log
    │   │       └── S2.log
    │   └── qc
    │       ├── detailed_fragmentation_length_plot
    │       │   ├── S1.log
    │       │   └── S2.log
    │       ├── extract_fastq_insertions
    │       │   ├── S1.log
    │       │   └── S2.log
    │       ├── extract_mapping_quality
    │       │   ├── S1.log
    │       │   └── S2.log
    │       ├── finalize_mapping_quality
    │       │   ├── S1.log
    │       │   └── S2.log
    │       ├── fragmentation_distribution_plots
    │       │   ├── fragmentation_match_distribution_S1.log
    │       │   ├── fragmentation_match_distribution_S2.log
    │       │   ├── fragmentation_read_match_distribution_S1.log
    │       │   └── fragmentation_read_match_distribution_S2.log
    │       ├── generate_mapq_heatmap
    │       │   ├── S1.log
    │       │   └── S2.log
    │       ├── multiqc
    │       │   └── out.log
    │       ├── nanoplot
    │       │   ├── S1.log
    │       │   └── S2.log
    │       └── read_level_fastqc
    │           ├── S1.log
    │           └── S2.log
    ├── mapping
    │   ├── Postcut_S1.bed
    │   ├── Postcut_S1_sorted.bam
    │   ├── Postcut_S1_sorted.bam.bai
    │   ├── Postcut_S1_unfiltered_sorted.bam
    │   ├── Postcut_S1_unfiltered_sorted.bam.bai
    │   ├── Postcut_S2.bed
    │   ├── Postcut_S2_sorted.bam
    │   ├── Postcut_S2_sorted.bam.bai
    │   ├── Postcut_S2_unfiltered_sorted.bam
    │   ├── Postcut_S2_unfiltered_sorted.bam.bai
    │   ├── Precut_S1.bed
    │   ├── Precut_S1_sorted.bam
    │   ├── Precut_S1_sorted.bam.bai
    │   ├── Precut_S2.bed
    │   ├── Precut_S2_sorted.bam
    │   └── Precut_S2_sorted.bam.bai
    └── qc
        ├── fastqc
        │   ├── readlevel_S1
        │   │   ├── S1_read_Read-343.fastq
        │   │   ├── S1_read_Read-343_fastqc.html
        │   │   ├── S1_read_Read-343_fastqc.zip
        │   │   ├── S1_read_Read-555.fastq
        │   │   ├── S1_read_Read-555_fastqc.html
        │   │   ├── S1_read_Read-555_fastqc.zip
        │   │   ├── S1_read_Read-561.fastq
        │   │   ├── S1_read_Read-561_fastqc.html
        │   │   ├── S1_read_Read-561_fastqc.zip
        │   │   ├── S1_read_Read-745.fastq
        │   │   ├── S1_read_Read-745_fastqc.html
        │   │   ├── S1_read_Read-745_fastqc.zip
        │   │   ├── S1_read_Read-902.fastq
        │   │   ├── S1_read_Read-902_fastqc.html
        │   │   └── S1_read_Read-902_fastqc.zip
        │   ├── readlevel_S2
        │   │   ├── S2_read_Read-262.fastq
        │   │   ├── S2_read_Read-262_fastqc.html
        │   │   ├── S2_read_Read-262_fastqc.zip
        │   │   ├── S2_read_Read-417.fastq
        │   │   ├── S2_read_Read-417_fastqc.html
        │   │   ├── S2_read_Read-417_fastqc.zip
        │   │   ├── S2_read_Read-522.fastq
        │   │   ├── S2_read_Read-522_fastqc.html
        │   │   ├── S2_read_Read-522_fastqc.zip
        │   │   ├── S2_read_Read-682.fastq
        │   │   ├── S2_read_Read-682_fastqc.html
        │   │   ├── S2_read_Read-682_fastqc.zip
        │   │   ├── S2_read_Read-824.fastq
        │   │   ├── S2_read_Read-824_fastqc.html
        │   │   └── S2_read_Read-824_fastqc.zip
        │   ├── S1_filtered.fastq
        │   └── S2_filtered.fastq
        ├── multiqc_data
        │   ├── multiqc_citations.txt
        │   ├── multiqc_data.json
        │   ├── multiqc_fastqc.txt
        │   ├── multiqc_general_stats.txt
        │   ├── multiqc.log
        │   ├── multiqc_nanostat.txt
        │   ├── multiqc_software_versions.txt
        │   └── multiqc_sources.txt
        ├── multiqc_report.html
        └── nanoplot
            ├── S1
            │   ├── AlignedReadlengthvsSequencedReadLength_dot.html
            │   ├── AlignedReadlengthvsSequencedReadLength_dot.png
            │   ├── AlignedReadlengthvsSequencedReadLength_kde.html
            │   ├── AlignedReadlengthvsSequencedReadLength_kde.png
            │   ├── MappingQualityvsReadLength_dot.html
            │   ├── MappingQualityvsReadLength_dot.png
            │   ├── MappingQualityvsReadLength_kde.html
            │   ├── MappingQualityvsReadLength_kde.png
            │   ├── NanoPlot_20250103_1228.log
            │   ├── NanoPlot-report.html
            │   ├── NanoStats.txt
            │   ├── Non_weightedHistogramReadlength.html
            │   ├── Non_weightedHistogramReadlength.png
            │   ├── Non_weightedLogTransformed_HistogramReadlength.html
            │   ├── Non_weightedLogTransformed_HistogramReadlength.png
            │   ├── PercentIdentityHistogramDynamic_Histogram_percent_identity.html
            │   ├── PercentIdentityHistogramDynamic_Histogram_percent_identity.png
            │   ├── PercentIdentityvsAlignedReadLength_dot.html
            │   ├── PercentIdentityvsAlignedReadLength_dot.png
            │   ├── PercentIdentityvsAlignedReadLength_kde.html
            │   ├── PercentIdentityvsAlignedReadLength_kde.png
            │   ├── WeightedHistogramReadlength.html
            │   ├── WeightedHistogramReadlength.png
            │   ├── WeightedLogTransformed_HistogramReadlength.html
            │   ├── WeightedLogTransformed_HistogramReadlength.png
            │   ├── Yield_By_Length.html
            │   └── Yield_By_Length.png
            └── S2
                ├── AlignedReadlengthvsSequencedReadLength_dot.html
                ├── AlignedReadlengthvsSequencedReadLength_dot.png
                ├── AlignedReadlengthvsSequencedReadLength_kde.html
                ├── AlignedReadlengthvsSequencedReadLength_kde.png
                ├── MappingQualityvsReadLength_dot.html
                ├── MappingQualityvsReadLength_dot.png
                ├── MappingQualityvsReadLength_kde.html
                ├── MappingQualityvsReadLength_kde.png
                ├── NanoPlot_20250103_1228.log
                ├── NanoPlot-report.html
                ├── NanoStats.txt
                ├── Non_weightedHistogramReadlength.html
                ├── Non_weightedHistogramReadlength.png
                ├── Non_weightedLogTransformed_HistogramReadlength.html
                ├── Non_weightedLogTransformed_HistogramReadlength.png
                ├── PercentIdentityHistogramDynamic_Histogram_percent_identity.html
                ├── PercentIdentityHistogramDynamic_Histogram_percent_identity.png
                ├── WeightedHistogramReadlength.html
                ├── WeightedHistogramReadlength.png
                ├── WeightedLogTransformed_HistogramReadlength.html
                ├── WeightedLogTransformed_HistogramReadlength.png
                ├── Yield_By_Length.html
                └── Yield_By_Length.png

66 directories, 219 files


    
```

</details>

<br>

### Output files

#### 1. Localization

The sequence-guided detection of insertions is the core of the workflow. In addition to simply identifying the insertions, several other interesting parameters are automatically evaluated during the execution of the pipeline.

##### Genomic location 


File: `../final/localization/ExactInsertions_{sample}.bed`

**Simulated S1:**
```plaintext

    chr1	270204	272451	Read-561	[257666, 291832]	+
    chr1	314899	323644	Read-343	[296872, 297968]	+
    chr1	432141	440886	Read-902	[428005, 432140]	+
        
```

!!! warning

    The `strand` column in `ExactInsertions_{sample}.bed` refers to the alignment of the read, not the insertion itself.

!!! info 
    
    This file is the primary output and shows the positions of the detected insertions, which are dependent on the reference. It follows the standard [BED6](https://samtools.github.io/hts-specs/BEDv1.pdf) format with the columns: `Chromosome - Start - End - Read - Original Read Start/End - Strand`.
 

##### Orientation and structure

In addition to the main output, it can be useful to examine the orientation of the insertion and the exact structure of the inserted sequence within the read.
    
File: `../final/qc/Fragmentation/Longest_Interval/{sample}/Longest_interval_{read}.bed`

**S1 Read-343:**

![Longest_interval_Read-343](images/tutorial/Longest_interval_Read-343.png)

The small numbers displayed above the line represent the matching vector fragments, while the x-axis indicates the actual length in base pairs (bp) of the longest consecutive interval.

The longest detected interval of this read contained all possible 100 bp vector fragments from 0 to 87, with ambiguous 100 bp matches in the region around positions 6/7 and 55/56 of the insertion sequence. This ambigous region of the insertion corresponds to the long-terminal reapeats (LTRs) of the [vector construct](other.md/#vector-map). 

!!! info

    Since the underlying vector sequence FASTA is in the 5'-3' orientation, and this order is maintained in the longest-matching interval of the fragmented sequence, the insertion and the read share the same `+` orientation. 

**S2 Read-262:**

![Longest_interval_Read-536](images/tutorial/Longest_interval_Read-262.png)

The small numbers displayed above the line represent the borders of the matching vector fragments, while the x-axis indicates the actual length in base pairs (bp) of the interval.

The longest consecutively detected interval of this read included only a subset of all 100 bp vector fragments, resulting in a shorter insertion of approximately 2500 bp. Additionally, the fragment numbers appear to be detected in descending order.

!!! info 
    
    Since the insertion sequence FASTA is oriented in the 5'-3' direction, and this order is **not** preserved in the longest-matching interval of the fragmented sequence, the insertion in the read has a `-` orientation. This indicates that the vector sequence is located in the `-` orientation on a `+` directional read.


#### 2. Quality control

The workflow automatically assesses the quality of the input sequencing data, the alignments performed with and without fragmentation, and the fragmentation itself. This allows not only for detecting insertions but also for evaluating the likelihood of true positives and the overall effectiveness of the search strategy employed by the pipeline. 

##### Input data quality

The pipeline integrates basic quality assessment tools from widely established resources, including [FastQC](), [MultiQC](), and [NanoPlot](). An overview of the results can be accessed via Snakemake's workflow report, which is generated using `snakemake --report` or directly in the output directory.

File: `../final/qc/multiqc_report.html`

!!! info 
    
    The pipeline uses fastqc by processing the FASTQ of each read with a detected insertion individually. 

!!! Hint "Further Details"
    
    For detailed explanations of the plots provided in the report, consult the documentation of each quality control tool. To access the individual quality control results, navigate to the following directories within the output folder:

    fastqc: `../intermediate/qc/fastqc/`<br>
    multiqc: `../intermediate/qc/multiqc/`<br>
    nanoplot: `../intermediate/qc/nanoplot/`<br>
 
##### Mapping quality

The pipeline incorporates two mapping steps to improve the quality of mapping by modifying reads that contain insertions. These steps are essential for accurately localizing the insertions, making it crucial to track the mapping quality of the affected reads at each key alignment stage. 

 File: `../intermediate/qc/mapq/Insertions_{sample}_mapq.txt`

**S1:** 

```plaintext
Read	    PrecutChr	        PrecutMAPQ	PostcutChr	PostcutMAPQ	FilteredChr	FilteredMAPQ
Read-343	pSLCAR-CD19-CD3z	60	        chr1	    44	        chr1	    44.0
Read-555	pSLCAR-CD19-CD3z	60	        *	        0		
Read-561	chr1	            60	        chr1	    60	        chr1	    60.0
Read-745	pSLCAR-CD19-CD3z	60	        *	        0		
Read-902	pSLCAR-CD19-CD3z	60	        chr1	    60	        chr1	    60.0
```
    
The table illustrates changes in mapping quality and chromosome alignment for each read with an insertion across three stages: **Precut** mapping before any modifications, **Postcut** mapping after the reads were modified, and **Filtered** mapping after filtering based on mapping quality. 
    
!!! info 
    
    During the initial mapping of the unaltered reads, four out of the five reads containing detected insertions predominantly aligned with high quality to the vector reference. However, after the modification (`Buffer`), where every base of the insertion was replaced with `N`, two additional reads successfully mapped to a region in the reference genome, while the other two reads became unmappable.

!!! info

    The scores from the table are automatically visualized in the plot. However, due to overlapping quality scores, some reads may be obscured by others with identical values. In the example data, this occurs with `Read-902` and `Read-561`, as well as for `Read-555` and `Read-745`.

    **S1:**

    ![S1_mapq_heatmap_image.png](images/tutorial/S1_mapq_heatmap_image.png)


##### Fragmentation
The fragmentation process is a crucial step not only for detecting insertions but also for gaining a detailed understanding of the exact composition and orientation of the inserted sequence. Some aspects of fragmentation quality control align closely with the analysis of the [orientation and structure](#orientation-and-structure) of the detected insertions.

However, the analysis of the previously mentioned output files overlooks another critical factor: The existence of fragments with significant sequence similarity to other "normal" sequences in the reference FASTA.

The pipeline includes functionality to perform a BLASTN search of the fragmented insertion sequence against a pre-built version of your reference's BLAST database. To enable this feature, simply specify the `blastn_db` argument in the `config.yml`. 

!!! Danger 
    
    The potential similarity of the insertion sequence to other sequences in your reference is particularly important when using the pipeline in conjunction with complex vector expression systems. For example, CAR T cell vector constructs (like our example vector [construct](other.md/#insertion-sequence) ) often insert sequences partially derived from human genes.

As this option is not configured for the tutorial, we can instead rely on two other automatically generated plots to gain insights into potential false-positive matches for the insertion sequence. 

Directory: `../final/qc/Fragmentation/Insertions_{fragmentsize}_{sample}/`
    
![Combined_Barplots.png](images/tutorial/Combined_Barplots.png)
    
These two plots illustrate the distributions of all insertion fragments (left) and the number of fragment matches "contributed" by each read (right).

!!! info 
    
    The `Combined distribution of all 100 bp fragments` plot reveals that every fragment of the vector is represented at least four times. However, fragments `6`,`7`,`55`, and `56` are noticeably overrepresented in the reads. As mentioned in the [orientation and structure](#orientation-and-structure) section, these fragments correspond to the vector's LTRs, making their alignment ambiguous. The slight plateau observed between fragments `57` and `78` is better understood in conjunction with the second plot. 
    
    The `Contribution of reads to the toal count of 100 bp fragments` plot clarifies this plateau by showing the read-specific contributions of fragments. Four reads contribute the maximum number of vector fragments, whereas `Read-561` includes only about 21 vector fragments. This leads to the slight overrepresentation of fragments `57` to `78` in the `Combined distribution of all 100 bp fragments` plot. 

!!! Attention 
    
    Observations like these are crucial for determining the most accurate `MinInsertionLength` threshold in the `config.yml`. 

!!! Hint "Further Details"
    
    For the example data, we selected only a very small portion of the reference genome to generate reads. This is why there are no additional "off-target" fragment matches within our reads. Since the vector construct contains several human-derived components in its architecture, a real sequencing dataset would likely result in a more complex barplot.

    As mentioned before, the safest way to identify potential misleading fragment matches in advance is to provide a human BLASTn reference database to the pipeline. The vector fragments are then automatically aligned against this reference, and the resulting plots offer an overview of the vector regions that are highly likely to appear, even in the absence of an actual insertion.

    <details><summary> S1 Barplots when provided a BLASTn reference: </summary>
        ![Combined_BlastN_Barplots.png](images/tutorial/Combined_BlastN_Barplots.png)

    The bar plots now illustrate which vector fragments are likely to produce false positives. When comparing these fragments with the structure of the [construct](other.md/#insertion-sequence), you can identify three main regions of fragment matches: fragments `22`–`25` correspond to the EF-1a core promoter, fragments `35`–`39` align with CD28 and CD247, and fragments `56`–`57` represent the 3'LTR. These are all human components in the vector architecture that we can also anticipate detecting with the pipeline by using the vector genome as the target sequence.
   
    </details>
    
<br>

#### 3. Functional annotation 
Typically, identifying the genomic localization of an insertion is just the starting point, with the proximity to the insertion site being the next area of interest. A basic yet essential functionality for annotating the detected insertion sites is included in the pipeline through the `functional_genomics.smk` rule collection. The pipeline can work with up to four different user-defined BED annotation files that can be provided in the `config.yml` as `annotation_1`, `annotation_2`, `annotation_3`, and `annotation_4`.

##### Genes in proximity

For the tutorial, we have defined only one annotation file in the `config.yml`, which simply contains the known genes located in our specified reference FASTA. For details on generating this file, refer to [this](other.md/#simulate-data-for-tutorial).The pipeline compares the locations of the insertions with the entries in the provided annotation file and reports the closest match, producing the file shown below.

File: `../final/functional_genomics/Functional_distances_to_Insertions_{sample}.bed`

**S1:**

| InsertionChromosome | InsertionStart | InsertionEnd | InsertionRead | InsertionOrig      | InsertionStrand | AnnotationChromosome | AnnotationStart | AnnotationEnd | AnnotationID    | AnnotationScore | AnnotationStrand | AnnotationSource                             | Distance |
|---------------------|----------------|--------------|---------------|--------------------|-----------------|----------------------|-----------------|---------------|-----------------|-----------------|------------------|---------------------------------------------|----------|
| chr1                | 270204         | 272451       | Read-561      | [257666, 291832]   | +               | chr1                 | 266854          | 268655        | ENSG00000286448 | .               | +                | UCSC_genes_chr1_0_500000_processed          | -1550    |
| chr1                | 314899         | 323644       | Read-343      | [296872, 297968]   | +               | chr1                 | 360056          | 366052        | ENSG00000236601 | .               | +                | UCSC_genes_chr1_0_500000_processed          | 36413    |
| chr1                | 432141         | 440886       | Read-902      | [428005, 432140]   | +               | chr1                 | 450739          | 451678        | OR4F29          | .               | -                | UCSC_genes_chr1_0_500000_processed          | 9854     |

 
A good starting point to get familiar with the personalization of the pipeline tailored to your specific research question can be including a rule for the visualisation of this table. Check out the [advanced usage](advanced_usage.md/#custom-insertion-annotations) for more on this. 

!!! Hint "Further Details"   

    The reads for this tutorial are artificially generated based on the first `50kb` of sequence from human chromosome 1. The regions at the beginning of chromosomes (near the centromeres and telomeres) are typically less gene-dense compared to the more gene-rich areas toward the middle of the chromosomes. This relative scarcity of coding genes also makes these regions less accessible for the integration of lentiviral-based vector systems, thus reducing the biological plausibility of our simulated data.


#### 4. Intermediate files

The workflow generates numerous additional files beyond those listed above. Most of these files are straightforward to understand once you are familiar with the pipeline's functionality. They are typically not critical for most use cases unless [debugging](#error-handling) is required or you integrate [custom downstream rules](advanced_usage.md/#developer-mode) into the analysis.  

Directory: `../intermediate/`

!!! info
    Here is a list of each subdirectory and a description of what to find in them: 

    **`blastn/`**

    - `Filtered_Annotated_{fragmentsize}_InsertionMatches_{sample}.blastn`: Results from the BLASTn searches after filtering 
    - `Coordinates_{fragmentsize}_InsertionMatches_{sample}.blastn`: Dictionary of the identified FASTA coordinates based on insertions in the reads
    - `ref/`: BLASTN matches of vector fragments with provided ref blastdb (empty files if no `blast_db` provided)


    **`fasta/`**

    - `fragments/`: Constructed BLASTN database based on the query insertion
    - `Modified_{sample}_mod.fa`: Modified FASTA file of input BAM (read modification dependent on `Buffer`, `Split`, or `Join`)
    - `Full_{sample}.fa`: Unmodifed FASTA file of input BAM
    - `Insertion_{sample}.fa`: Detected insertion sequences extacted from the reads

    **`localization/`**

    - `ExactInsertions_{sample}.bed`: File as in final output

    **`log/`**

    - See [Error handling](#log-files)

    **`mapping/`**

    - `Precut_{sample}_sorted.bam`: Unmodified Reads after reference mapping
    - `Postcut_{sample}_unfiltered_sorted.bam`: (Modified) Reads after reference mapping
    - `Postcut_{sample}_sorted.bam`: (Modified) Reads passing the quality filter after reference mapping
    - `Postcut_{sample}_sorted.bed`: Genomic locations of aligned reads

    **`qc/`**

    - `fastqc/`: Fastqc input and raw output
    - `multiqc_data/`: Multiqc raw output
    - `nanoplot/`: Nanoplot raw output
    - `multiqc_report.html`: Report as in final output

---

## Error handling

#### Snakemake
General debugging ressources for everything related to snakemake can be found in the snakemake [FAQ](https://snakemake.readthedocs.io/en/v6.15.5/project_info/faq.html).

#### Log files
The pipeline is designed with rule-specific `log` files, which are stored in the `intermediate` output directory. These logs serve as the primary resource for identifying and addressing any rule-specific issues that arise during execution. If you encounter errors or unexpected behavior, these files should be your first point of reference for debugging. 