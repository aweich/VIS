# Getting Started
#### Prerequisites

!!! info

    Make sure you have a working installation of [conda](https://conda.io/projects/conda/en/latest/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/) (**recommended**)

Clone and navigate into the repository:
   
```bash
   git clone https://github.com/aweich/VIS
   cd VIS
```
---

#### Setting up the Environment

You have two options for setting up the pipeline:

##### 1. Classic (recommended)

The pipeline contains the option for an [integrated package management](https://snakemake.readthedocs.io/en/latest/snakefiles/deployment.html#integrated-package-management) functionality. 
In other words, the workflow will be installing all required packages for each rule as needed.
This means that we only need to install a lightweight environment ourselves to run the pipeline locally.

!!! Info 

    While the traditional method of setting up the pipeline offers the advantage of being lightweight in its package installation, the actual execution of the pipeline can be slightly more complex and may require additional adjustments based on your local setup of `conda` or `mamba`. To address this, we also provide a convenient alternative that installs all required packages for the pipeline within a single environment below.

Create the minimal environment for running the pipeline: 

```bash 
    mamba env create --name VIS_minimal -f workflow/envs/VIS_minimal_env.yml
```

!!! Info 

    You can use `conda` instead of `mamba` if desired, but `mamba` is faster. 

Then, activate the environment: 

```bash
   conda activate VIS_minimal
```
##### 2. Alternative

`Snakemake` also allows us to use a single environment throughout the whole pipeline. We can also make use of that to create a convenient alternative for the setup of the pipeline. Depending on your system, this might make the downsteam execution of the workflow overall more convenient. 

Create the full environment for running the pipeline:

```bash    
   mamba env create --name VIS_full -f workflow/envs/VIS_full_env.yml
```

!!! Info 
    It is also possible to use `conda` for the creation of the environment. Since we are installing more dependencies here, the use of `mamba` is highly encouraged.

Then, activate the environment:

```bash
   conda activate VIS_full 
```