
<p align="center">
    <img src="/images/logo.png" alt="Logo" width=200>
</p>

# Vector Insertion Site – Detection Pipeline
----
* [Introduction](#introduction)
* [Key Features and Applications](#key-features-and-applications)
* [Workflow](#workflow)
* [General Usage](#general-usage)
* [Citation and Contribution](#citation-and-contribution)
* [License](#license)

## Introduction

Welcome to the **Vector Insertion Site (VIS) Detection Pipeline** documentation. This `Snakemake`-based workflow detects and annotates insertion sites in long-read DNA sequencing data. It supports custom functions and follows the `Snakemake` styling guide.  

The diagram below outlines the detection workflow, with key analysis steps on the left and a directed acyclic graph (DAG) of workflow components on the right. *Created in BioRender. [Weich, A. (2025)](https://BioRender.com/z40d414)*  

For a detailed explanation, see our [paper](). In brief, a (partially) known vector sequence is fragmented into kmers and searched for matches in long-read sequencing data. Matching reads are modified and mapped to a reference genome. The **<span style="color:black">detection</span>** of the exact insertion site is implemented using a CIGAR-based reverse calculation and can be **<span style="color:darkblue">functionally annotated</span>** with genomic resources (e.g., genes, transcription factors). Throughout the workflow, multiple **<span style="color:darkred">quality control</span>** steps (e.g., base quality, mapping quality) ensure integrity of input data and results.  


#### Illustrated Core Functionality
<p align="center">
    <a href="/images/Combined_Workflow_for_Documentation.png" target="_blank">
        <img src="/images/Combined_Workflow_for_Documentation.png" alt="Workflow overview">
    </a>
</p>

<br>

## Key Features and Applications

- **Systematic Detection of Inserted Vector Sequences:** The pipeline fragments the target sequence into bins, enabling precise identification of consecutive and isolated vector fragments, including their orientation in long DNA reads.

- **Accurate Localization of Insertion Sites:** By mapping the detected sequences to a reference genome, the pipeline provides exact genomic coordinates of the insertion sites, ensuring high accuracy.

- **Comprehensive Annotation Capabilities:** Leveraging customziable genome annotations, the pipeline determines the proximity of insertion sites to transcription factors, genes, or other biologically relevant elements, allowing for interpretations about genomic alterations caused by the VIS.

- **Customizable and Extendable Framework:** Designed with flexibility in mind, the pipeline allows users to add custom analyses, annotations, or visualizations, catering to specific research needs and enhancing functionality.

- **Applicability to Biomedical Research:** Particularly useful for CAR T cell therapy, the pipeline can, for instance, identify clonal insertions that could lead to therapeutic complications. In addition, it can also be used to analyze lentiviral and retroviral insertion patterns, detect guided insertions, and support other gene delivery studies.

## General Usage

Everything from installation to customization of this pipeline can be found in the [online documentation](). If you are new to `Snakemake`, check out the [introduction to snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) first. 

It is recommended to familiarize yourself with the workflow and its outputs before running it with your own data. A detailed practical example of the workflow and its output files with simulated data can be found in the [tutorial](tutorial/tutorial_intro.md).

## Quick Start

Follow these steps to get the pipeline up and running:

Clone the Repository:

```bash
git clone https://github.com/aweich/VIS
cd VIS
```

Create and Activate the Environment:

```bash
mamba env create --name VIS_minimal -f workflow/envs/VIS_minimal_env.yml
conda activate VIS_minimal
```

Run the Pipeline:

```bash
snakemake --use-conda -n
```

For more detailed instructions, refer to the [online documentation]() and the corresponding [publication]().

## Citation and Contribution

If you are using this pipeline or the search strategy implemented in this pipeline, please cite us using "..." and leave us a star. 

We encourage contributions to enhance and refine this codebase, whether through providing feedback, improving functionality, or sharing domain-specific expertise. If you have suggestions, encounter issues, or require assistance, please feel free to reach out for support or collaboration.

## License
This project is available under the [Apache License 2.0](LICENSE).
