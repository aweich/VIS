
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
Welcome to the documentation for the **Vector Insertion Site (VIS) Detection Pipeline**. This `snakemake`-based workflow helps researchers detect and annotate insertion sites with known sequences in long-read DNA sequencing data. The workflow is fully extensible with customized functions and adheres to the Snakemake styling guide. The basic detection workflow is illustrated in [Figure 1](#figure-1-illustrated-core-functionality) with some of the key analysis steps. A detailed structure for the Snakemake pipeline can be found in [Figure 2](#figure-2-dag-of-the-workflow), where a directed acyclic graph (DAG) of the individual rules provides an overview of the workflow components.

#### Figure 1: Illustrated Core Functionality
<p align="center">
    <img src="/images/Workflow_for_Documentation.svg" alt="Workflow overview" width=500>
</p>

<br>

## Key Features and Applications

- **Systematic Detection of Inserted Vector Sequences:** The pipeline fragments the target sequence into manageable bins, enabling precise identification of consecutive and isolated vector fragments, including their orientation in long DNA reads.

- **Accurate Localization of Insertion Sites:** By mapping the detected sequences to a reference genome, the pipeline provides exact genomic coordinates of the insertion sites, ensuring high accuracy.

- **Comprehensive Annotation Capabilities:** Leveraging customziable genome annotations, the pipeline determines the proximity of insertion sites to transcription factors, genes, or other biologically relevant elements, offering a deeper understanding of potential effects.

- **Customizable and Extendable Framework:** Designed with flexibility in mind, the pipeline allows users to add custom analyses, annotations, or visualizations, catering to specific research needs and enhancing functionality.

- **Applicability to Biomedical Research:** Particularly useful for CAR T cell therapy, the pipeline can, for instance, identify clonal insertions that could lead to therapeutic complications. In addition, it can also be used to analyze lentiviral and retroviral insertion patterns, detect guided insertions and support other gene delivery studies.

## Workflow
Below is the DAG of the Snakemake workflow. Rules are colored according to the different parts of the workflow illustrated in [Figure 1](#figure-1-illustrated-core-functionality). 

Detailed explanations of the detection pipeline can be found [here](). In brief, a known target insertion (i.e. vector) is fragmented into smaller DNA sequences and used for a sequence similarity search against reads from long-read DNA sequencing. Matching reads are subjected to one of three types of modifications (buffer, split, join). These modified reads are subsequently mapped against their respective reference genome, and the exact genomic location of the insertion is calculated using a CIGAR-based approach. Finally, genomic annotations are utilized to calculate the distances of specific elements (genes, transcription factors, etc.) to the detected insertion sites.

#### Figure 2: DAG of the Workflow
<p align="center">
    <img src="/images/DAG_for_Documentation.svg" alt="Pipeline DAG" >
</p>

## General Usage

Everything from installation to customization of this pipeline can be found in the [documentation](index.md). If you are new to `snakemake`, check out the [introduction to snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) first. 

It is recommended to familiarize yourself with the workflow and its outputs before running it with your own data. A detailed description of the workflow and its output files with simulated data can be found in the [tutorial](tutorial/tutorial_intro.md).

## Citation and Contribution

If you are using this pipeline or the search strategy implemented in this pipeline, please cite us using "github.com/aweich/VIS-dp". 

We encourage contributions to enhance and refine this codebase, whether through providing feedback, improving functionality, or sharing domain-specific expertise. If you have suggestions, encounter issues, or require assistance, please feel free to reach out for support and collaboration.

## License
This project is available under the [Apache License 2.0](LICENSE).

---


