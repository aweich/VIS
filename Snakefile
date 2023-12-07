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
		expand(PROCESS+"MAPPING/Normalisation_IPHM_{sample}.txt", sample=SAMPLES),
		#expand(PROCESS+"MAPPING/BasicMapping_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed" , sample=SAMPLES),
		#Methylation
		expand(PROCESS+"METHYLATION/Methyl_{sample}.bed", sample=SAMPLES),
		expand(PROCESS+"METHYLATION/Insertion_fasta_proximity_{sample}.bed", sample=SAMPLES),
		expand(PROCESS+"METHYLATION/MeanMods_Proximity_{sample}.bed", sample=SAMPLES),
		#Visuals
		expand(PROCESS+"LOCALIZATION/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(PROCESS+"BLASTN/HUMANREF/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(PROCESS+"METHYLATION/Heatmap_MeanMods_Proximity_{sample}.png", sample=SAMPLES),
		#PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png",
		#PROCESS+"LOCALIZATION/Heatmap/",
		#deeper
		expand(PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES)

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

rule nBases_for_insertion_count:
	input:
		get_input_names
	output:
		temp(PROCESS+"MAPPING/Number_of_Bases_{sample}.normalisation")
	shell:
		"gatk CountBases -I {input} > {output}"				

rule normalisation_for_insertion_count:
	input:
		insertions=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed",
		number_of_bases=PROCESS+"MAPPING/Number_of_Bases_{sample}.normalisation"
	params:
		scale=100000000 #IPHM
	output:
		PROCESS+"MAPPING/Normalisation_IPHM_{sample}.txt"
	run:
		vhf.insertion_normalisation(input.insertions, input.number_of_bases, params[0], output[0])

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
	
rule reads_with_BLASTn_matches:
	input:
		refbed=PROCESS+"MAPPING/BasicMapping_{sample}.bed",
		matchreads=PROCESS+"BLASTN/ID/ID_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed" 
	shell:
		"grep -F -f {input.matchreads} {input.refbed} > {output}"


#Methylation
rule prepare_BAM:
	input:
		get_input_names
	output:
		PROCESS+"MAPPING/{sample}_sorted.bam"
	shell:
		"""
		samtools sort {input} -o {output}
		samtools index {output}
		"""
			
rule methylation_bedMethyl:
	input:
		bam=PROCESS+"MAPPING/{sample}_sorted.bam",
		#ref=config["ref_genome"]
		
	output:
		PROCESS+"METHYLATION/Methyl_{sample}.bed"
	shell:
		"modkit pileup {input.bam} --filter-threshold C:0.8 {output}" #high threshhold for C modifications
		#"modkit pileup {input.bam} {output} --cpg --ref {input.ref}" 

rule insertion_proximity:
	input:
		PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed"
	params: 
		10000 #this value has to be same as in mean methylation!
	output:
		PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed"
	run:
		vhf.insertion_proximity(input[0], params[0], output[0])


#creates fasta only with insertions (+proximity of 10k up and downstream
rule insertion_methylation_proximity:
    input:
        prox = PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed",
        fasta = config["ref_genome"]
    output:
        PROCESS+"METHYLATION/Insertion_methylation_proximity_{sample}.fa"
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.prox} > {output}"

#adds fasta to genomic gaps
rule fasta_to_insertion_proximity:
    input:
        prox = PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed",
        fasta = PROCESS+"METHYLATION/Insertion_methylation_proximity_{sample}.fa"
    output:
        PROCESS+"METHYLATION/Insertion_fasta_proximity_{sample}.bed"
    run:
        vhf.add_sequence_column(input.prox, input.fasta, output[0])

#this rule is added for performance, so that the following pythons cript only reads in relevant parts of the methylation file and not the full genome

rule methyl_specific:
    input:
        insertions = PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed",
        methyl = PROCESS+"METHYLATION/Methyl_{sample}.bed"
    output:
        PROCESS+"METHYLATION/Specific_Methyl_{sample}.bed"
    shell:
        """
        bedtools intersect -wb -a {input.methyl} -b {input.insertions} > {output}
        """

rule mean_methylation:
	input:
		methbed=PROCESS+"METHYLATION/Specific_Methyl_{sample}.bed",
		insertionfastabed=PROCESS+"METHYLATION/Insertion_fasta_proximity_{sample}.bed"
	params:
		window_size=500,
		max_distance=10000
	output:
		PROCESS+"METHYLATION/MeanMods_Proximity_{sample}.bed"
	run: 
		vhf.methylation_in_insertion_proximity(input.methbed, input.insertionfastabed, params.window_size, params.max_distance,  output[0])  
#Visuals	
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
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		FRAG
	output:
		outpath=directory(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}"),
		outpath2=directory(PROCESS+"BLASTN/HUMANREF/PLOTS/" + str(FRAG)+"_{sample}")
	run:
		shell("mkdir {output.outpath}")
		vhf.fragmentation_match_distribution(input[0], params[0], output[0])
		vhf.fragmentation_read_match_distribution(input[0], params[0], output[0])
		shell("mkdir {output.outpath2}")
		vhf.fragmentation_match_distribution(input[1], params[0], output[1])
		vhf.fragmentation_read_match_distribution(input[1], params[0], output[1])

rule insertion_heatmap:
	input:
		expand(PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed", sample=SAMPLES)
	output:
		PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png"
	run:
		print(input)
		vhf.plot_bed_files_as_heatmap(input, output[0])

rule insertion_modification_heatmap:
	input:
		PROCESS+"METHYLATION/MeanMods_Proximity_{sample}.bed"
	output:
		PROCESS+"METHYLATION/Heatmap_MeanMods_Proximity_{sample}.png"
	run:
		vhf.plot_modification_proximity(input[0], output[0])
		
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
		
'''
#deeper: Sniffles for variant calling: May be BLAST alternative?
rule variant_sniffles:
	input:
		PROCESS+"MAPPING/{sample}_sorted.bam"
	output:
		PROCESS+"MAPPING/Variant_{sample}.vcf"
	shell:
		"sniffles -i {input} -v {output}"
'''
