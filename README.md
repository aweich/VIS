# Logo with the name
#github action buttons via snakemake unit tests
# Vector Insertion Site Detection Pipeline
> This is a snakemake workflow for the guided-search for long insertion sequences (e.g. CAR vectors) in long-read DNA Sequencing data.

## Table of Contents
* [General Info](#general-information)
* [Workflow](#workflow)
* [Setup](#setup)
    * [Pipeline](#Pipeline)
    * [Config](#Config)
* [Usage](#usage)
* [Tutorial](#tutorial)
* [Contact](#contact)
* [License](#license)

## General Information
- Biological vector systems are widely used for the targeted expression of genes in biomedical research and therapies (e.g. CAR T cell therapy)
- While insertion sites of different gene delivery systems are classified, there is no concise way to look for the inserted vector sequence and its exact sequence after the transduction
- This pipeline fragments the target sequence of interest into n-sized bins and systematically locates consecutive and isolated vector fragments and their orientation in the long DNA reads
- Furthermore, the insertion sites are located with respect to the reference genome of the underlying organism and, upon providing of annotations of the genome, closeness to transcription factors, genes, or other data sources can be used to annotate each insertion site
- Especially CAR T cell therapy can benefit from such a detailed approach, since the T cells are introduced back into the patients after the transduction and potentially dangerous clones are thereby introduced to the patient. This might be one explanantion of the recent case reports regarding secondary T cell Lymphomas following CAR T cell therapy.

## Workflow
- Overview (create it in a more universal manner though!) + DAG
- Illustration of the workflow

---

## Setup

## Pipeline

This snakemake workflow combines multiple standalone tools in a directed fashion. 

Clone the github repo to your local machine via 

    git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY

Then, proceed by one of the two ways: 

### 1. Classic (snakemake best practices)

The pipeline contains the option for an [integrated package management](https://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management) functionality. 
In other words, the workflow will be installing all required packages for each rule as needed.
This means that we only need to install a lightweight environment ourselves to run the pipeline locally.

**Note:** While the classical way of setting up the pipeline has the advantage of being lightweight in its installation of packagages, the actuall running of the pipeline will be a little more complicated and might require additional modifications depending on your local setup of `conda`/`mamba`. For this reason, we also provide a convenient alternative which works by installing all packages needed for the pipeline in one environment (see below).

Create the minimal environment for running the pipeline via 
    
    mamba env create --name VIS_minimal -f workflow/envs/VIS_minimal_env.yml

**Note:** It is also possible to use `conda` for the creation of the environment (`mamba` is way faster though). 

Then, activate the environment via 

    conda activate VIS_minimal 
 
### 2. Alternative

Snakemake also allows us to use a single environment throughout the whole pipeline. We can also make use of that to create a convenient alternative for the setup of the pipeline. Depending on your system, this might make the downsteam execution of the workflow overall more convenient. 

Create the full environment for running the pipeline via
    
    mamba env create --name VIS_full -f workflow/envs/VIS_full_env.yml

**Note:** It is also possible to use `conda` for the creation of the environment (`mamba` is way faster though). Since we are installing more dependencies here, the use of `mamba` is highly encouraged.

Then, activate the environment via 

    conda activate VIS_full 

## Config
The configuration file is needed to specify which information and files the pipeline needs and where it can find them. You can find the .yaml config used in the tutorial in `config/config.yaml`. A description of the inputs can be seen below. 

Most of the required inputs are mandatory, except for multiple addiitonal BED fiels for the annotation of the identified insertions and the `blastn_db`. 
The different threshholds are described in detail in the tutorial for this workflow. 

| **Parameter**            | **Value**                                                                                      | **Comments**                                                                                              |
|---------------------------|------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| `experiment`             | `tutorial`                                                                                  | Name of the experiment.                                                                                 |
| `samples`                | `S1: /path/to/S1.bam` <br> `S2: /path/to/S2.bam` | Paths to BAM files for the samples. Each sample should have its file path defined.                      |
| `processing_dir`         | `/path/to/outdir`                                                                           | Directory where output files will be saved.                                                             |
| `vector_fasta`           | `/path/to/vectorseq.fa`                                                      | Path to the insertion sequence file in FASTA format.                                                       |
| `blastn_db`              | `/path/to/blastNdb`                                                             | BLASTN database for reference nucleotides to check matches with the insertion sequence. If the insertion sequence contains reference genome parts, the matches here are needed as a baseline.                          |
| `splitmode`              | `Buffer`, `Split`, or `Join`                                                                                      | Mode used for processing reads based on the insertion: <ul><li><strong>Buffer:</strong> Replaces the insertion with "N" and retains the full length of each read. Recommended for the most accurate insertion location (CIGAR-based).</li><li><strong>Split:</strong> Cuts the insertion from the read and creates individual reads from the remaining sequence. The locations of the insertions are reported as the locations of each of the individual reads. Recommended for exploratory search for insertions that might fuse otherwise non-neighboring parts of the genome together.</li><li><strong>Join:</strong> Cuts the insertion from the read and joins the remaining sequence together. The locations of the insertions are reported as the location of the joined read. Recommended for debugging of otherwise unmappable reads.</li></ul> |
| `fragment_size`          | `100`                                                                                         | Size of fragments for splitting sequences. Fragments of this size will be used to construct the BLASTN database of the insertion sequence.                                                              |
| `MinLength`              | `1`                                                                                           | Minimum read length for BLASTN matches processing.                                                      |
| `MAPQ`                   | `10`                                                                                          | Minimum mapping quality score for reads.                                               |
| `MinInsertionLength`     | `500`                                                                                         | Minimum length of insertions to be detected. This is dependent on the respective insertion and potentially its matches with the reference genome.                                                            |
| `ref_genome_ctrl`        | `/path/to/ref.fa`                                                     | Reference genome file in FASTA format.                                                      |
| `ucsc_Genes`             | `/path/to/annotation.bed`                 | Annotation file in [BED6 format](https://samtools.github.io/hts-specs/BEDv1.pdf).                                                          |
| `detection`              | `rules/detection.smk`                                                                         | Snakemake rule file for the detection and localization of insertions.                                                                      |
| `quality_control`        | `rules/qc.smk`                                                                                | Snakemake rule file for the quality control of reads and insertions.                                                                |
| `functional_genomics`    | `rules/functional_genomics.smk`                                                               | Snakemake rule file for the functional annotation of insertions on the genome level.                                                            |

## Usage

After the setup, you are ready to run the pipeline on your specific files and use-cases. Generally, the execution follows the same principles as defined in the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/). However, here is a quick overview about the basic commands for execution, since they slightly differ dependent on the mode of the setup: 

In general, the classic setup needs to be run using `--use-conda` to tell snakemake to build and use the rule-specific envs. 

Exemplarily for both, a dry-run can be executed as follows:

```Classic setup
snakemake --use-conda -n
```

**Note:** If you run into problems with `mamba` or `conda` not being able to create the rule-specific environemnts on to go, check out [this](https://stackoverflow.com/questions/69001097/conda-4-10-3-and-snakemake-5-conda-exe-problem) as a starting point. 

```Alternative setup
snakemake -n
```

To illustrate these rules in a directed acyclic graph, run:

    snakemake --forceall --rulegraph | dot -Tpng > dag.png

To run the full pipeline using 20 cores and automatically generate a report about the run afterwards, run: 

    snakemake --cores 20 && snakemake --report

This will, however, fail since the sample names defined in the default config are dummies.

## Tutorial

A detailed description and an exemplary run with simulated data can be found [here]().

## Recommendations

- Run the pipeline iteratively first on a smaller subset of your data. It is very likely that you need to change the thresholds for the fragment size and MinInsertioNlength multiple times before you figure out the best configuration for your dataset.

## Contact
Created by []() - feel free to contact me!

## License
This project is available under the [... License]().

---

# Tutorial (will be moved to ./docs/tutorial at a later point)

## Setup
Either 

Clone the repository and create a new working environment: 

`git clone https://github.com/yourusername/VIS_pipeline.git`
   
`conda env create -f VIS_env.yml`
   
`conda activate VIS_pipeline

Also make sure that you have [blastn] () and [fastqc] () installed and executable in the previously generated environemnt.

Or

make use of the the modularization character of the pipeline:

`git clone the repository and and conda env create -f VIS_snakemake_base_env.yml´

Activate the environment and dry run it via 

`snakemake --cores 5  --use-conda -n`

In case your conda frontend is conda and not mamba, use it as:

`snakemake --cores 5  --use-conda --conda-frontend conda -n`

For each rule, the pipeline as a dedicated env assigned, if it requires anything that is not already present in VIS_basic_env. 
On the first execution of the pipeline, the requirements of each env will be checked and if needed, packages will be installed. 

For convenience, it is also possible to initialize the full environment once and run the pipeline from within by using  

`snakemake --cores 5 -n`

For clarity, we will use this type of execution for the rest of the tutorial.

## Collect the reference files

To now execute the pipeline for the first time, we need to make sure all dependencies required in the config are defined. Open the `config.yml` file and add the missing paths.

For a minimum reproducible example, you can leave the config as it is for now. 
samples:
    Simulated1: "Simulated.bam"



## Quick start

`snakemake -n`

To illustrate these rules in a directed acyclic graph, run:

`snakemake --forceall --rulegraph | dot -Tpng > dag.png`

To run the full pipeline using 20 cores and automatically generate a report about the run afterwards, run: 
`snakemake --cores 20 && snakemake --report`

## General usage

### Running the pipeline
#### Modify config
#### Run 
 After we have installed everything, it's time to check what the pipeline will be doing after the execution. Let's look at a dry-run of the workflow:
 
 `snakemake -n`
 
Snakemake will run all of these rules on our provided files. Let's check the DAG again to see whether we can understand the route that snakemake will take:

snakemake --forceall --rulegraph | dot -Tpng > dag.png`

Okay, that's enough preparation. Let's run it!

The execution time of the pipeline depends on the amount of physical cores that you can allocate to the workflow. Let's choose a low number for now:

`snakemake --cores 5 && snakemake --report`

Error handling:
You can now follow the different rules in your terminal window. If you encounter errors, make sure to double-check your initial input. If your error is for a specific rule, check the detailed documentation for this step in ./log/rulename

Other general debugging ressources for everything related to snakemake can be found [here] () or [here] ().

If you see the following, the pipeline was executed successfully.
img

So let's take a look at the output. 

### Inspecting the output
#### check report

Since we ran the pipeline with `&& snakemake --report`, we have also automatically generated a general report for the workflow. 
Take a look at the Statistics in report.html. The mapping rules clearly took the longest to finish.

Throughout the pipeline, some simple plots are also generated to give you a glimpse of what your insertions look like. 
Navigate to the results tab and take a look at the detected lengths of your insertions. Some of the lengths are clearly longer than others. It looks like some of the reads only contained parts of the insertion. 

To get an idea about quality control metrics, navigate to the multiqc.html report in the results tab. 

#### check final output files

maybe run `tree´ for the directory structure


Next to this overview, there are also other files generated that can subsequently be used for further investigations. The directory structure provided will pre-sort the outputs into intermediate and final output. Navigate to the final output/final/localization first.

Here, the most interesting one, of course, is the sample-specific BED file in containing the genomic positions of the insertions. This is also used for the annotation of each insertion (see functional genomics).

In output/final/qc, the subfolder fragmentation contains detailed information about each sample's fragmented insertion coverage, i.e. which parts of the insertion sequence was detected in the reads, what are the consecutive intervals for each read's insertion, and if there were also parts of the insertions detected, that match the human reference genome. In summary, this output gives you an idea about the effectiveness of the fragmentation of your inserted sequence for the detection of the read. Depending on the size of the insertion, you might want to run the analysis pipeline again and change the fragment-length in the config to further increase your level of detail.   

In the output/intermediate folder, you can find various subdirectories with intermediate files created during the pipeline run. Most of them are self-explanatory when looking at the workflow of the pipeline. To make things easier, an illustration about which files are created were in the DAG, please see below.  

Congratulations, you have finished the quick start using simulated data! If you want to unleash the full power of the pipeline, feel free to continue with the advanced usage tutorial below.


### functional genomics


## Adcanced usage
#### Custom thresholds
#### Custom insertion annotations
#### Developer mode
- Pre-select output
- Separated or Join functionality
- Add custom rules


to do:

Create minimum reproducible example:
- 2 simulated samples with 10 insertions each; some +, some -, randomly spreach across the reference genome -> make the files small)
- create small lightweight gtf gene reference
- create small lightweight blastndb
- only use chromosome 1 for example for know: randomly smaple reads from chr1 and merge 10 insertions in there (+ and -)
- also add chr1 hg38.fa then


