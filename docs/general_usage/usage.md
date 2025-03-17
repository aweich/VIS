# Usage

After the setup, you are ready to run the pipeline on your specific files and use-cases. Generally, the execution follows the same principles as defined in the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/). However, here is a quick overview about the basic commands for execution, since they slightly differ depending on the chosen mode of the setup: 

In general, the classic setup needs to be run using `--use-conda` to tell `Snakemake` to build and use the rule-specific envs. 

Exemplarily for both, a dry-run can be executed as follows:

Classic:
```bash
   snakemake --use-conda -n
```
Alternative:
```bash
   snakemake -n
```

!!! Tip 
    
    If you run into problems with `mamba` or `conda` not being able to create the rule-specific environments on the go, you might first need to allow `conda` or `mamba` to be able to call subprocesses. Check out [this](https://stackoverflow.com/questions/69001097/conda-4-10-3-and-snakemake-5-conda-exe-problem) as a starting point. 

#### Example commands

To illustrate the order of execution for these rules in a directed acyclic graph, run:
```bash
   snakemake --forceall --rulegraph | dot -Tpng > dag.png
```

To run the full pipeline using 20 cores and automatically generate a report about the run afterwards, run: 

```bash
   snakemake --cores 20 && snakemake --report myreport.html
```

This, however, will fail since the sample names defined in the default config are dummies. Let's look at the [tutorial](../tutorial/tutorial_intro.md) for a step-by-step explanation of the the workflow and its outputs. 