import os
import time
from pathlib import Path

configfile: "config.yml"
SRC=config["source_dir"]
SAMPLES = expand(config["samples"]) 
PROCESS = os.path.join(config["processing_dir"],str(config["experiment"]+"/")) #intermediate and results files are stored here

#local functions - path to helper fucntions needs to be added to the sys path, otherwise import won't find the file
rootpath = os.path.join(SRC)
sys.path.append(rootpath)
from helper_functions import parse_FASTQ #functions to make snakemake pipeline leaner

		
rule all:
	input: 
		expand(PROCESS+"FASTA/{sample}.fa", sample=SAMPLES) 

#actual filenames
def get_input_names(wildcards):
    return config["samples"][wildcards.sample]

rule make_FASTA:
	input:
		fq=get_input_names
	output:
		fasta=PROCESS+"FASTA/{sample}.fa"
	run: 
		shell("seqkit fq2fa {input.fq} -o {output.fasta}") #seqkit had to be installed with conda

#prerequisites: build blastdb with vector sequence https://www.ncbi.nlm.nih.gov/books/NBK569841/
#check output file: coordinates of matches? Can I further use them to split the reads?
rule find_vector_BLASTn:
	input:
		fasta=PROCESS+"FASTA/{sample}.fa"
	params:
		vector=config["pro_nucleotidedb"] #vector db
	output:
		#output matches
	run:
		shell("blastn -query {input} -db {params.pro_blastdb} -out {output.pro} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs' -num_threads 10") 

	
