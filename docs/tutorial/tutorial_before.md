## Before running the pipeline
---
### Prepare config.yml

To inform the pipeline about the location of our samples and other dependencies, we need to open and edit the configuration file. Open or create a new `config.yml` file and fill in the missing dependencies as outlined below:

```yaml
# tutorial config
experiment: "simulation_tutorial"
samples:
  S1: "tutorial/simulated/S1.bam"
  S2: "tutorial/simulated/S2.bam"
processing_dir: "tutorial/out"
threads: 2
insertion_fasta: "tutorial/references/vectorseq.fa"
splitmode: "Buffer"
fragment_size: 100
bridging_size: 300
MinLength: 1
MAPQ: 10
MinInsertionLength: 500
ref_genome_ctrl: "tutorial/references/chr1_1_50000_ref.fa"
annotate_ucsc_genes: "tutorial/references/UCSC_genes_chr1_0_500000_processed.bed"
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

If you see a list of jobs waiting for [execution](./tutorial_running.md#expected-jobs), you're all set for the next steps.

!!! info 
    
    For this tutorial, we will run the workflow using the recommended setup. If you have chosen the alternative version (i.e., one virtual environment (`VIS_full`)), simply remove `--use-conda` from following `Snakemake` commands. 

<br>
