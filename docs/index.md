# Vector Insertion Site Detection Pipeline

Welcome to the documentation for the **Vector Insertion Site Detection Pipeline**. This Snakemake workflow helps researchers detect and annotate insertion sites of long sequences in DNA sequencing data.

![Pipeline DAG](images/dag.png)

## Features
- Systematic search for vector sequences.
- Annotates insertion sites with reference genome information.
- Configurable thresholds and flexible modes.

---


## Getting Started

### Prerequisites
1. Install [conda](https://conda.io/projects/conda/en/latest/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/).
2. Install the pipeline:
   ```bash
   git clone https://github.com/YOUR-USERNAME/YOUR-REPOSITORY
   cd YOUR-REPOSITORY
´´´
### Setup the Environment

You have two options for setting up the pipeline:

### 1. Classic (snakemake best practices)

The pipeline contains the option for an [integrated package management](https://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management) functionality. 
In other words, the workflow will be installing all required packages for each rule as needed.
This means that we only need to install a lightweight environment ourselves to run the pipeline locally.

**Note:** While the classical way of setting up the pipeline has the advantage of being lightweight in its installation of packagages, the actuall running of the pipeline will be a little more complicated and might require additional modifications depending on your local setup of `conda`/`mamba`. For this reason, we also provide a convenient alternative which works by installing all packages needed for the pipeline in one environment (see below).

Create the minimal environment for running the pipeline via 

```bash
   mamba env create --name VIS_minimal -f workflow/envs/VIS_minimal_env.yml
´´´

**Note:** You can use ´conda´ instead of ´mamba´ if desired, but ´mamba is faster. 

Then, activate the environment via 

```bash
   conda activate VIS_minimal 
´´´
 
### 2. Alternative

Snakemake also allows us to use a single environment throughout the whole pipeline. We can also make use of that to create a convenient alternative for the setup of the pipeline. Depending on your system, this might make the downsteam execution of the workflow overall more convenient. 

Create the full environment for running the pipeline via

```bash    
   mamba env create --name VIS_full -f workflow/envs/VIS_full_env.yml
´´´

**Note:** It is also possible to use `conda` for the creation of the environment (`mamba` is way faster though). Since we are installing more dependencies here, the use of `mamba` is highly encouraged.

Then, activate the environment via 
```bash
   conda activate VIS_full 
´´´

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

