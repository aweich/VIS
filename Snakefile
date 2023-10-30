import os
import time
from pathlib import Path
from time import sleep


configfile: "config.yml"
SRC=config["source_dir"]
SAMPLES = expand(config["samples"]) 
PROCESS = os.path.join(config["processing_dir"],str(config["experiment"]+"/")) #intermediate and results files are stored here
FRAG=config["fragment_size"]

#local functions - path to helper fucntions needs to be added to the sys path, otherwise import won't find the file
rootpath = os.path.join(SRC)
sys.path.append(rootpath)
#print(rootpath)	
import VIS_helper_functions as vhf #functions to make snakemake pipeline leaner
		
rule all:
	input: 
		expand(PROCESS+"FASTA/{sample}.fa", sample=SAMPLES),
		expand(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"),
		expand(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa{ext}", ext=[ ".ndb",".nhr",".nin",".not",".nsq",".ntf",".nto"]), 
		expand(PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		expand(PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES)

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

#fragments the CAR construct reference FASTA to increase detection rates via BLASTn
rule vector_fragmentation:
	input: 
		config["vector_fasta"] #vector fasta sequence
	params:
		FRAG
	output: 
		fasta=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	run:
		vhf.fragmentation_fasta(input[0], params[0], output[0])

	
rule make_BLASTN_DB: #does this work without an output file? Usually, if another rule refers to the output of a rule, the file has to be specified
	input:
		vector_fragmented = PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		multiext(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
	run:
		shell("makeblastdb -in {input.vector_fragmented} -dbtype nucl -blastdb_version 5")


#prerequisites: build blastdb with vector sequence https://www.ncbi.nlm.nih.gov/books/NBK569841/
#check output file: coordinates of matches? Can I further use them to split the reads?
rule find_vector_BLASTn:
	input:
		fasta=PROCESS+"FASTA/{sample}.fa"
	params:
		#vector=config["blastn_db"] #vector db
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:
		shell("blastn -query {input} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs' -num_threads 20") 

rule hardcode_blast_header:		
	input: 
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:	
		"echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input} > {output}"

