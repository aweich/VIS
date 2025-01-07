<p align="center">
    <img src="docs/images/logo.png" alt="Logo" width=200>
</p>

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

<p align="center">
    <img src="/docs/introduction/images/Workflow_for_Documentation.svg" alt="Workflow overview" width=300>
</p>

Detailed explanations of the detection pipeline can be found [here](). In brief, a known target insertion (i.e. vector) is fragmented into smaller DNA sequences and used for a sequence similarity search against reads from long-read DNA sequencing. Matching reads are subjected to one of three types of modifications (buffer, split, join). These modified reads are subsequently mapped against their respective reference genome, and the exact genomic location of the insertion is calculated using a CIGAR-based approach. Finally, genomic annotations are utilized to calculate the distances of specific elements (genes, transcription factors, etc.) to the detected insertion sites.

## General Usage

Everything from installation to customization for this pipeline can be found [here](). If you are new to snakemake, check out the general idea behind it [here](). After installation, it is advised to first familiarize yourself with the workflow by following the [tutorial](), where every output file is introduced.   

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

A detailed description and an exemplary run with simulated data can be found [readthedocslink?]().

## Contact & Contribution

## License
This project is available under the [Apache License 2.0](LICENSE).

---


