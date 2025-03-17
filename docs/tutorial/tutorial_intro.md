# Tutorial

Now, let's explore a minimal example of what the pipeline can accomplish. To demonstrate this, we have simulated sequencing data and randomly introduced insertions into some of the samples. For more details on the data, refer to [this section](../other/other_simulation.md).  

#### Background  

Let's assume we have two samples, S1 and S2, which are cells that were previously transduced using a chimeric antigen receptor (CAR) vector construct. We then extracted DNA from these samples and performed long-read DNA sequencing, resulting in BAM files for both samples.  

Our goal is to determine whether the transduction was successful—i.e., whether the transgene was incorporated into the genome—and whether the vector insertion occurred within any known human gene.  

#### Structure  

###### Before Running the Pipeline  
- [Prepare config](tutorial_before.md/#prepare-configyml)  
- [Check setup](tutorial_before.md/#check-setup)  

###### Running the Pipeline  
- [Expected jobs](tutorial_running.md/#expected-jobs)  
- [Execution](tutorial_running.md/#execution)  

###### After Running the Pipeline  
- [Snakemake report](tutorial_after.md/#snakemake-report)  
- [Output directory structure](tutorial_after.md/#output-directory-structure)  
- [Output files](tutorial_after.md/#output-files)  
