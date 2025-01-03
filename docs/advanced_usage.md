## Adcanced usage
### Customization
#### Custom thresholds

There are several possible adjustments that likely need to be made to the thresholds of the pipeline to ensure the best possible tradeoff between sensitivity and specificity.

| Parameter             | Exemplary value | Description                                                                                                                                         | Considerations                                                                                                                                                            |
|-----------------------|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `fragment size`       | 100            | Size of fragments (in bp) for splitting sequences. Fragments of this size will be used to construct the BLASTN database of the insertion.          | Consider the mean read length of the data and the expected insertion size. How likely is it to get a full insertion contained within a full read? Are there many parts of the insertion that can also be expected to be found elsewhere in the reference genome? |
| `bridging_size`       | 300             | Acceptable gap size (in bp) before splitting the longest consecutive interval. It is recommended to customize this parameter according to the underlying read quality and insertion length. | Consider the bridging size in relation to the fragment size. A 100bp fragment should definitely be bridged if the expected insertion is multiple kb long. |
| `MinLength`           | 1             | Minimum read length (in bp) for BLASTN matches processing.                                                                                          | Consider the overlap of the `fragment_size` with the insertion length. A `fragment_size` of 1000bp with an insertion length of 2050 will create two 1000bp fragments and one 50bp fragment. Keep in mind that the BLASTN results in the pipeline are **always** subjected to quality filtering with `evalue < 1e-5` and `bitscore > 50`. |
| `MAPQ`                | 10              | Minimum mapping quality score for reads. Filtering is applied after the modification of the reads with insertions.                                  | Consider the quality of your reference genome.                                                                                                                                 |
| `MinInsertionLength`  | 500            | Minimum length (in bp) of insertions to be detected. This is dependent on the respective insertion and potentially its matches with the reference genome. | Consider the length of the regions with sequence similarity between the reference genome and the insertion reference.                                                              |

#### Custom annotations

Adding custom annotations for the detected insertion coordinates is straightforward. Simply include a sorted `BED6` file in the `config.yml` under one of the annotation keywords (`annotation_1`/`annotation_2`/`annotation_3`/`annotation_4`). For annotation data sourced from repositories like [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) or [GENCODE](https://www.gencodegenes.org/), , ensure the files are sorted before running the pipeline. For instance, the gene annotation file used in the [tutorial](tutorial.md/#genes-in-proximity) was generated following the steps outlined [here](other.md/#annotation-data).

!!! info
    In the `config.yml`, the `annotation_1` field is mandatory. The other annotation file slots are optional and can be used as needed. While it is possible to extend the number of annotation files beyond the predefined slots, it is recommended to implement such customizations by creating custom `rules.smk` files for better flexibility and maintainability..   

<br>

#### Custom rules and functions

One of Snakemake's greatest advantages over standalone software tools is its flexibility, allowing users to easily integrate custom options into any stage of the pipeline. Whether you want to add a new analysis step, build an entirely new branch of functionality, or simply generate additional plots, the pipeline can accommodate your needs.

To achieve this, however, you’ll need to gain a deeper understanding of Snakemake's core [core functionality](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html). This includes learning how to write custom `Snakefile` rules, define input/output relationships, and use configuration files effectively. By doing so, you can benefit from Snakemake's modular design to extend your workflow without compromising its reproducibility and scalability.

As a matter of fact, there is one rule already implemented that generates an additional plot based on the genes in proximity to the detected insertions. The output for this rule is currently just commented out with a `#` in the `rule all:` of the `Snakefile`. For a first step towards your own customization, you can remove the `#`, execute the pipeline, and see whether you can find the rule that generates the output file.  

<details>
  <summary>Output file, rule, and plotting script: </summary>
<br>

The output file will not be generated unless it is included in the <code>rule all</code>. To ensure the file is recognized, remove the <code>#</code> from the relevant line in <code>../workflow/Snakefile</code>.

```python
#expand(f"{outdir}/final/functional_genomics/Plot_Distance_to_Genes_{fragmentsize}_{{sample}}.png", sample=SAMPLES),
```

The rule that generates the output is located in <code>../workflow/rules/functional_genomics.smk</code>. 

```python
rule plot_distance_to_elements:
	input:
		distancetable=f"{outdir}/final/functional_genomics/Functional_distances_to_Insertions_{{sample}}.bed"
	params:
		distances=list(range(-10000, 10001, 2000)),
		threshold=10000
	log:
        	log1=f"{outdir}/intermediate/log/functional_genomics/plot_distance_to_elements/scatter_{{sample}}.log",
	output:
		scatter=report(f"{outdir}/final/functional_genomics/Plot_Distance_to_Genes_{fragmentsize}_{{sample}}.png"),
	run:
	    try:
	        vhf.plot_element_distance(input.distancetable, params.distances, params.threshold, output.scatter, log.log1)
	    except Exception as e:
	        with open(log.log1, "a") as log_file:
                    log_file.write(f"Error: {str(e)}\n")
```

You can find the corresponding plotting function in the <code>python</code> helper script <code>../workflow/scripts/VIS_helper_functions.py</code>:

```python
def plot_element_distance(bed, distances, distance_threshold, output_path, logfile):
    """
    Uses the bed file from the distance calculations and provides a plot to visualize the respective elements with their distance. Entries that are further away than the defined threshold value are excluded.
    """
    # Read the table
    df = pd.read_csv(
        bed,
        sep='\t',
    )

    # Apply threshold filtering if provided
    if distance_threshold is not None:
        df = df[df['Distance'].abs() <= int(distance_threshold)]
  
    # Ensure absolute distance and sort by absolute distance within groups
    df['abs_distance'] = df['Distance'].abs()
    df = df.sort_values(by=['InsertionRead', 'abs_distance']).drop_duplicates(subset=['InsertionRead', 'AnnotationID', "AnnotationSource"], keep='first').reset_index()
    
    # Create scatter plot
    sns.scatterplot(
        data=df,
        x='Distance',
        y='AnnotationID',
        hue='InsertionRead',
        palette='tab10',
        s=100,
        style='AnnotationSource'
    )
    
    # Binned rugplot for distances
    bin_size = 100  # Bin size grouping distances
    df['distance_bin'] = (df['Distance'] // bin_size) * bin_size
    sns.rugplot(x=df['distance_bin'], color='black', height=0.05, linewidth=1)
    
    # Configure plot aesthetics
    plt.xticks(sorted({x for n in distances for x in (n, -n)}), rotation=45)
    plt.xlabel("Distance (bp)")
    plt.ylabel("Element Name")
    plt.title("Distance Distribution to Elements")
    sns.despine()
    plt.legend(title="", fontsize=8)
    
    # Save plot
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Plot saved to {output_path}")
```
</details>

<br>

### Limitations

#### Sequencing depth

The pipeline detects insertions based on the available data: While the number of detected insertions may occasionally appear low, this often reflects the nature of the underlying sequencing data, which typically contains a limited number of insertions. Long-read sequencing is well-suited for detecting and precisely locating insertions at the base level. However, if the base-level probability of the insertion sequences is low, the number of detectable insertions may also be restricted.

#### Input file inconsistencies

The workflow strives to follow best practices for all data types whenever possible. However, this may not always work optimally if the input reference FASTAs, BED files, or BAM files do not strictly adhere to these practices, or if they contain unexpected symbols, which can lead to downstream issues following some kind of malfunctioning string split. Not all such cases can be anticipated in advance. If you encounter problems and suspect that input file names might be the cause, please let us know, ideally with a minimal reproducible example. 

#### Base-level accuracy for Split/Join splitmode

Currently, only the `Buffer` splitmode is fully implemented to allow insertions to be detected with base-level accuracy. In contrast, the `Split` and `Join` modes only report the coordinates of the aligned read(s) that contain the insertion(s). This functionality is expected to be improved in future versions.
