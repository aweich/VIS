import os
import time
from pathlib import Path
from time import sleep


configfile: "config.yml"
SRC=config["source_dir"]
SAMPLES = expand(config["samples"]) 
OUTDIR = os.path.join(config["processing_dir"],str(config["experiment"]+"/"))
PROCESS = OUTDIR+"intermediate/" #intermediate files are stored here
FINAL = OUTDIR+"final/" # final files are stored here
FRAG=config["fragment_size"]

#local functions - path to helper fucntions needs to be added to the sys path, otherwise import won't find the file
rootpath = os.path.join(SRC)
sys.path.append(rootpath)
#print(rootpath)
print(os.getcwd())
cwd=os.getcwd()    
import VIS_helper_functions as vhf #functions to make snakemake pipeline leaner

#inmport rules
include: config["functional_genomics"]
include: config["quality_control"]
#development
#include: config["epigenetics"]
#include: config["variants"]
#include: config["development"]

#target rule        
rule all:
	input: 
		#Localization
		expand(FINAL+"localization/ExactInsertions_{sample}.bed", sample=SAMPLES),
		FINAL+"localization/Heatmap_Insertion_Chr.png",
		FINAL+"localization/Insertion_length.png",
		#Functional
		#expand(FINAL+"functional_genomics/localization/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(FINAL+"functional_genomics/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png", sample=SAMPLES),
		expand(FINAL+"functional_genomics/Insertion_Scoring_{sample}.png", sample=SAMPLES),
		#Quality
		expand(PROCESS+"qc/mapq/{sample}_mapq_heatmap_image.png", sample=SAMPLES),
		expand(FINAL+"qc/Fragmentation/Insertions/insertions_" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(FINAL+"qc/Fragmentation/Longest_Interval/{sample}/", sample=SAMPLES),
		FINAL+"qc/multiqc_report.html",
		#Process
		OUTDIR+"config_settings.yml",
		#Other
		expand(PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.gff", sample=SAMPLES)
		
#actual filenames
def get_input_names(wildcards):
	return config["samples"][wildcards.sample]


#print("Results in ...")
#print(PROCESS)

######
######
###### Insertion BAM: Convert to fasta 
######
######

rule copy_config_version:
	input:
		config_file="config.yml"
	log:
		log=PROCESS+"log/main/copy_config_version/out.log"
	output:
		OUTDIR+"config_settings.yml" 
	shell:
		"""
		(
		cp {input.config_file} {output} 
		) > {log.log} 2>&1
		"""

rule build_insertion_reference:
	input:
		ref=config["ref_genome_ctrl"],
		vector=config["vector_fasta"]
	log:
		log=PROCESS+"log/main/build_insertion_reference/out.log"
	output:
		PROCESS+"mapping/vector_ref_genome.fa"
	shell:
		"""
		(
		cat {input.ref} {input.vector} > {output}
		) > {log.log} 2>&1
		"""
		

rule minimap_index:
	input:
		ref=PROCESS+"mapping/vector_ref_genome.fa"
	log:
		log=PROCESS+"log/main/minimap_index/out.log"
	output:
		index=temp(PROCESS+"mapping/ref_genome_index.mmi")
	resources:
		mem_mb=5000
	conda:
		"envs/VIS_minimap_env.yml"
	shell:
		"""
		(
		minimap2 -d {output.index} {input.ref} 
		) > {log.log} 2>&1
		"""

rule make_fasta_without_tags: #fasta of raw data no trimming whatsoever
	input:
		fq=get_input_names
	log:
		log=PROCESS+"log/main/make_fasta_without_tags/{sample}.log"
	output:
		fasta=PROCESS+"fasta/Full_{sample}.fa"
	conda:
		"envs/VIS_samtools_env.yml"
	shell: 
		"""
		(
		samtools fasta {input} -o {output} > {output}
		) > {log.log} 2>&1
		"""
######
######
###### "Clean" BAM: Cut-out fasta to BAM via Mapping to reference 
######
######
		
#rule to create the BAM files with the non-insertion reads and the splitted read fragments
rule Non_insertion_mapping: #mapping against the unaltered referenc egenome
	input:
		fasta=PROCESS+"fasta/Cleaved_{sample}_noVector.fa",
		genome=PROCESS+"mapping/vector_ref_genome.fa"
	output:
		PROCESS+"mapping/Postcut_{sample}_unfiltered_sorted.bam"
	log:
		log=PROCESS+"log/main/Non_insertion_mapping/{sample}.log"
	resources:
		mem_mb=5000
	conda:
		"envs/VIS_minimap_env.yml"
	shell: #N=0 instead of default N=1
		"""
		(
		minimap2 -y -ax map-ont --score-N 0 {input.genome} {input.fasta} | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		) > {log.log} 2>&1
		"""

rule insertion_mapping: #conserves tags!
	input:
		bam=get_input_names,
		minimapref=PROCESS+"mapping/ref_genome_index.mmi",
		ref=PROCESS+"mapping/vector_ref_genome.fa"
	output:
		PROCESS+"mapping/Precut_{sample}_sorted.bam"
	log:
		log=PROCESS+"log/main/insertion_mapping/{sample}.log"
	resources:
		mem_mb=5000
	conda:
		"envs/VIS_minimap_env.yml"
	shell:
		"""
		(
		samtools bam2fq -T '*' {input.bam}| minimap2 -y -ax map-ont {input.minimapref} - | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		) > {log.log} 2>&1
		"""

rule clean_postcut_by_maping_quality:
	input:
		PROCESS+"mapping/Postcut_{sample}_unfiltered_sorted.bam"
	params:
		mapq=config["MAPQ"]
	output:
		PROCESS+"mapping/Postcut_{sample}_sorted.bam"
	log:
		log=PROCESS+"log/main/clean_postcut_by_maping_quality/{sample}.log"
	conda:
		"envs/VIS_samtools_env.yml"
	shell:
		"""
		(
		samtools view -h -q {params.mapq} {input} -o {output}
		samtools index {output}
		) > {log.log} 2>&1
		"""
######
######
###### Genomic Coordinates of Reads with and without matches
######
######

rule BAM_to_BED:
	input:
		#get_input_names
		precut=PROCESS+"mapping/Precut_{sample}_sorted.bam",
		postcut=PROCESS+"mapping/Postcut_{sample}_sorted.bam" #Reads that contained an insertion before, are now marked with "_Buffer/Insertion/0,1,2..."
	output:
		postcut=PROCESS+"mapping/Postcut_{sample}.bed",
		precut=PROCESS+"mapping/Precut_{sample}.bed"
	log:
		log1=PROCESS+"log/main/BAM_to_BED/Precut_{sample}.log",
		log2=PROCESS+"log/main/BAM_to_BED/Postcut_{sample}.log"
	conda:
		"envs/VIS_bedtools_env.yml"
	shell:
		"""
		(
		bedtools bamtobed -cigar -i {input.precut} > {output.precut}
		) > {log.log1} 2>&1
		(
		bedtools bamtobed -cigar -i {input.postcut} > {output.postcut}
		) > {log.log2} 2>&1
		"""

######
######
###### fasta preparation: Cut out of blast-detected vector fragments and create new "cut-out" fasta 
######
######

rule get_cleavage_sites_for_fasta: #filters and combines matches
	input:
		PROCESS+"blastn/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		filteroption=True,
		filtervalue=config["MinInsertionLength"], 
		overlap=3*FRAG # 2*FRAG #this is the distance of the start-stop that is allowed to exist to still be combined; This should not be lower than FRAG!
	log:
		log=PROCESS+"log/main/get_cleavage_sites_for_fasta/{sample}.log"
	output:
		cleavage=PROCESS+"blastn/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		reads=PROCESS+"blastn/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt"
	run:
		vhf.splitting_borders(input[0],params.filteroption, params.filtervalue, params.overlap, output.cleavage, output.reads, log.log)

rule split_fasta:
	input:
		breakpoints=PROCESS+"blastn/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		fasta=PROCESS+"fasta/Full_{sample}.fa"
	params:
		mode=config["splitmode"] #if each split fasta substring should be used individually, use "Separated" Join, New mode: Buffer
	log:
		log=PROCESS+"log/main/split_fasta_by_borders/{sample}.log"
	output:
		fasta=PROCESS+"fasta/Cleaved_{sample}_noVector.fa",
		vector=PROCESS+"fasta/Insertion_{sample}_Vector.fa"
	run:
		vhf.split_fasta_by_borders(input.breakpoints, input.fasta, params.mode, output.fasta, output.vector, log.log)

######
######
###### Vector preparation: Fragmentation 
######
######
rule prepare_vector:
	input:
		config["vector_fasta"] #vector fasta sequence
	log:
		log=PROCESS+"log/main/prepare_vector/out.log"
	output: 
		fasta=PROCESS+"fasta/fragments/Forward_Backward_Vector.fa"
	run:
		vhf.reversevector(input[0], output[0], log.log)
		
rule vector_fragmentation:
	input:
		PROCESS+"fasta/fragments/Forward_Backward_Vector.fa" #does not change anything so it can be removed imo
	params:
		FRAG
	log:
		log=PROCESS+"log/main/vector_fragmentation/out.log"
	output: 
		fasta=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	run:
		vhf.fragmentation_fasta(input[0], params[0], output[0], log.log)

######
######
###### BLAST Searches - Against Vector and healthy human reference
######
######

rule make_blastn_DB:
	input:
		vector_fragmented = PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	log:
		log=PROCESS+"log/main/make_blastn_DB/out.log"
	output:
		multiext(PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa",
			".ndb",
			".nhr",
			".nin",
			".not",
			".nsq",
			".ntf",
			".nto"
		)
	conda:
		"envs/VIS_blastn_env.yml"
	shell:
		"""
		(
		makeblastdb -in {input.vector_fragmented} -dbtype nucl -blastdb_version 5 
		) > {log.log} 2>&1
		"""

rule find_vector_BLASTn:
	input:
		fasta=PROCESS+"fasta/Full_{sample}.fa",
		dummy=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		tempdir=PROCESS+"temp_{sample}",
		vector=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	log:
		log=PROCESS+"log/main/find_vector_BLASTn/{sample}.log"
	output:
		PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	conda:
		"envs/VIS_blastn_env.yml"
	shell:
		"""
		(
		mkdir {params.tempdir}
		
        blastn \
        -query {input.fasta} \
        -db {params.vector} \
        -out {params.tempdir}/temp_output.blastn \
        -evalue 1e-5 \
        -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'

        # Filter results based on bitscore > 50
        awk '$16 > 50' {params.tempdir}/temp_output.blastn > {output}

        # Clean up temporary files
        rm -r {params.tempdir}
        ) > {log.log} 2>&1
        """
		
#blastn vector against human genome: Which vector parts are close to human sequences so that they might raise a false positivite BLAST match
rule find_vector_BLASTn_in_humanRef:
	input:
		fasta=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	params:
		vector=config["blastn_db"] #full human reference
	log:
		log=PROCESS+"log/main/find_vector_BLASTn_in_humanRef/{sample}.log"	
	output:
		temp(PROCESS+"blastn/humanref/"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	conda:
		"envs/VIS_blastn_env.yml"
	shell:
		"""
		(
		blastn \
		-query {input} \
		-db {params.vector} \
		-out {output} \
		-evalue 1e-5 \
		-outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'
		) > {log.log} 2>&1
		"""

rule hardcode_blast_header:
	input: 
		vector=PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	log:
		log=PROCESS+"log/main/hardcode_blast_header/{sample}.log"	
	output:
		vector=PROCESS+"blastn/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:
		"""
		(
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.vector} > {output.vector}
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.humanref} > {output.humanref}
		) > {log.log} 2>&1
		"""
		
rule blast_to_gff:
	input:
		ref=PROCESS+"blastn/humanref/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		vector=PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		vector=PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.gff",
		ref=PROCESS+"blastn/humanref/"+str(FRAG)+"_VectorMatches_{sample}.gff"  
	run:
		vhf.blast2gff(input.ref, output.ref)
		vhf.blast2gff(input.vector, output.vector)  

#filter BLAST matches by min length
rule extract_by_length:
	input:
		blast=PROCESS+"blastn/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		threshold=config["MinLength"]
	log:
		log=PROCESS+"log/main/extract_by_length/{sample}.log"	
	output:
		blast=PROCESS+"blastn/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:
		"""
		(
		awk -F'\t' '$11>={params.threshold}' {input.blast} > {output.blast}
		awk -F'\t' '$11>={params.threshold}' {input.humanref} > {output.humanref}
		) > {log.log} 2>&1
		""" 

######
######
###### Visualisation of insertions
######



rule basic_insertion_plots:
	input:
		expand(PROCESS+"localization/ExactInsertions_{sample}.bed", sample=SAMPLES)
	output:
		report(FINAL+"localization/Heatmap_Insertion_Chr.png"),
		report(FINAL+"localization/Insertion_length.png")
	log:
		log1=PROCESS+"log/main/basic_insertion_plots/heat.log",
		log2=PROCESS+"log/main/basic_insertion_plots/length.log"
	run:
		vhf.plot_bed_files_as_heatmap(input, output[0], log.log1)
		vhf.plot_insertion_length(input, output[1], log.log2)

######
######
###### Exact localization
######
######

#exact coordinates of the matching fragments

rule calculate_exact_insertion_coordinates:
	input:
		bed=PROCESS+"mapping/Postcut_{sample}.bed",
		borders=PROCESS+"blastn/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		mode=config["splitmode"]
	log:
		log=PROCESS+"log/main/calculate_exact_insertion_coordinates/{sample}.log"
	output:
		out=PROCESS+"localization/ExactInsertions_{sample}.bed",
	run:
		vhf.reconstruct_coordinates(input.bed, input.borders, params.mode, output.out, log.log)

rule collect_outputs:
	input:
		coordinates=PROCESS+"localization/ExactInsertions_{sample}.bed",
	log:
		log=PROCESS+"log/main/collect_outputs/{sample}.log"
	output:
		coordinates=FINAL+"localization/ExactInsertions_{sample}.bed"
	shell:
		"""
		(
		cp {input.coordinates} {output.coordinates}
		) > {log.log} 2>&1
		"""
