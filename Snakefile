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
		expand(PROCESS+"FASTA/{sample}.fa", sample=SAMPLES),
		expand(PROCESS+"BLASTN/WS5_Annotated_VectorMatches_{sample}.blastn", sample=SAMPLES)

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

''' #to do
rule vector_fragmentation:
	input: ???
	params:
		config["fragment_size"]
	output: #whatever is defined here will be used as a string !!!!
	run:
	fragmentation_fasta(input[0], params[0], output[0]) #output needs to be defined in script! 
	
rule make_BLASTn_DB: #does this work without an output file? Usually, if another rule refers to the output of a rule, the file has to be specified
	input:
		fasta = PROCESS+"FASTA/{sample}.fa"
	output:
		multiext(PROCESS+"FASTA/{sample}.fa",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
	run:
		shell("makeblastdb -in {input.fasta} -dbtype nucl -blastdb_version 5")
'''

#prerequisites: build blastdb with vector sequence https://www.ncbi.nlm.nih.gov/books/NBK569841/
#check output file: coordinates of matches? Can I further use them to split the reads?
rule find_vector_BLASTn:
	input:
		fasta=PROCESS+"FASTA/{sample}.fa"
	params:
		vector=config["blastn_db"] #vector db
	output:
		PROCESS+"BLASTN/WS5_VectorMatches_{sample}.blastn"
	run:
		shell("blastn -query {input} -db {params.vector} -out {output} -evalue 1e-5 -word_size 5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs' -num_threads 10") 

rule hardcode_blast_header:		
	input: 
		PROCESS+"BLASTN/WS5_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/WS5_Annotated_VectorMatches_{sample}.blastn"
	shell:	
		"echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input} > {output}"

