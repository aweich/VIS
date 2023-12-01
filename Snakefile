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
		#expand(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"),
		expand(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa{ext}", ext=[ ".ndb",".nhr",".nin",".not",".nsq",".ntf",".nto"]), 
		expand(PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		expand(PROCESS+"MAPPING/BasicMapping_{sample}.qc", sample=SAMPLES),
		#expand(PROCESS+"MAPPING/BasicMapping_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed" , sample=SAMPLES),
		#Methylation
		#expand(PROCESS+"METHYLATION/temp_{sample}/", sample=SAMPLES), #call methylation rule has to be dependent on the index rule. That's why the output is used as a fake input
		#expand(PROCESS+"METHYLATION/{sample}/{sample}_Methylation_pattern_regionBLABLA.tsv", sample=SAMPLES),
		#Visuals
		expand(PROCESS+"LOCALIZATION/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		#PROCESS+"LOCALIZATION/Heatmap/",
		#deeper
		#expand(PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES)

#actual filenames
def get_input_names(wildcards):
    return config["samples"][wildcards.sample]

rule make_FASTA:
	input:
		fq=get_input_names
	output:
		fasta=PROCESS+"FASTA/{sample}.fa"
	run: 
		#shell("seqkit fq2fa {input.fq} -o {output.fasta}") #seqkit had to be installed with conda
		shell("samtools fasta {input} -o {output} > {output}") #takes in bam and outputs fasta, output has to be mentioned twice: 1 for the location, 2 for the option to write "both" reads to one fasta

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

	
rule make_BLASTN_DB:
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

rule find_vector_BLASTn:
	input:
		fasta=PROCESS+"FASTA/{sample}.fa",
		dummy=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		#vector=config["blastn_db"] #vector db
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		temp(PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	run:
		shell("blastn -query {input.fasta} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs'") 

rule hardcode_blast_header:		
	input: 
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:	
		shell("echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input} > {output}")
#the next 3 rules need to be changed: Create bed from bam file, but no mapping necessary in between (since this should be already done)
#mapping without changes to the fasta files
#can be removed after chris acceptance
rule basic_mapping:
	input:
		fasta=get_input_names, # PROCESS+"FASTA/{sample}.fa",
		genome=config["ref_genome"] #cut _index
	output:
		PROCESS+"MAPPING/BasicMapping_{sample}.bam"
	shell:
		"""
		minimap2 -x map-ont -a {input.genome} {input.fasta} | samtools sort -o {output} 
		samtools index {output}
		"""   # alignment map-ont specifies input/task
		
		
		
#add quality ctrl rule
rule mapping_qc:
	input:
		get_input_names
		#PROCESS+"MAPPING/BasicMapping_{sample}.bam"
	output:
		PROCESS+"MAPPING/BasicMapping_{sample}.qc"
	run:
		shell("samtools flagstats {input} > {output}")  
		
rule BAM_to_BED:
	input:
		get_input_names
		#PROCESS+"MAPPING/BasicMapping_{sample}.bam"
	output:
		PROCESS+"MAPPING/BasicMapping_{sample}.bed"
	run:
		shell("bedtools bamtobed -i {input} > {output}")  

rule get_read_identifiers:
	input:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/ID/ID_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell: 
		"cut -f 1  {input} > {output}"

#this can stay				
rule reads_with_BLASTn_matches:
	input:
		refbed=PROCESS+"MAPPING/BasicMapping_{sample}.bed",
		matchreads=PROCESS+"BLASTN/ID/ID_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed" 
	shell:
		"grep -F -f {input.matchreads} {input.refbed} > {output}"

#Methylation
#this has to change or be removed completely: Nanopolish is deprecated for our flow cells
"""
rule index_for_methylation:
	input:
		fast5=config["fast5path"],
		fastq=get_input_names,
		#summary=config["sequencingsummary"]
	output:
		outpath=directory(PROCESS+"METHYLATION/temp_{sample}/")
	run:
		shell("mkdir {output.outpath}")
		shell("nanopolish index -d {input.fast5} {input.fastq}")

rule call_methylation:
	input:
		fastq=get_input_names,
		fake=PROCESS+"METHYLATION/temp_{sample}/",
		bam=PROCESS+"MAPPING/BasicMapping_{sample}.bam",
		ref=config["ref_genome"]
	output:
		PROCESS+"METHYLATION/{sample}/{sample}_Methylation_pattern_regionBLABLA.tsv"
	shell:
		"nanopolish call-methylation -r {input.fastq} -b {input.bam} -g {input.ref} -w 'chr6:1,000,000-20,000,000'> {output}"
"""
#Visuals #they can stay		
rule chromosome_read_plots:
	input:
		bam=get_input_names,
		#bam=PROCESS+"MAPPING/BasicMapping_{sample}.bam",
		bed=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed"
	output:
		outpath=directory(PROCESS+"LOCALIZATION/PLOTS/" + str(FRAG)+"_{sample}")
	params:
		buffer=50000
	shell: 
		r"""
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -buffer {params.buffer} -o {output.outpath}   
		"""	
rule fragmentation_distribution_plots:
	input:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		FRAG
	output:
		outpath=directory(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}")
	run:
		shell("mkdir {output.outpath}")
		vhf.fragmentation_match_distribution(input[0], params[0], output[0])
		vhf.fragmentation_read_match_distribution(input[0], params[0], output[0])

#deeper: BLASTN vector against human genome to see which parts might be matching in the UTD
rule find_vector_BLASTn_in_humanRef:
	input:
		fasta=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	params:
		vector=config["blastn_db"] #full human reference
	output:
		temp(PROCESS+"BLASTN/HUMANREF/"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	run:
		shell("blastn -query {input} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs'")
rule hardcode_blast_header_humanRef:		
	input: 
		PROCESS+"BLASTN/HUMANREF/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:	
		shell("echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input} > {output}")
