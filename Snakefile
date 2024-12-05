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
import VIS_helper_functions as vhf #functions to make snakemake pipeline leaner

#inmport rules
include: config["functional_genomics"]
include: config["quality_control"]
#development
include: config["epigenetics"]
include: config["variants"]
include: config["development"]

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
	output:
		OUTDIR+"config_settings.yml" 
	shell:
		"cp {input.config_file} {output}"

rule build_insertion_reference:
	input:
		ref=config["ref_genome_ctrl"],
		vector=config["vector_fasta"]
	output:
		PROCESS+"mapping/vector_ref_genome.fa"
	shell:
		"cat {input.ref} {input.vector} > {output}"
		

rule minimap_index:
	input:
		ref=PROCESS+"mapping/vector_ref_genome.fa"
	output:
		index=temp(PROCESS+"mapping/ref_genome_index.mmi")
	resources:
		mem_mb=5000
	shell:
		"minimap2 -d {output.index} {input.ref}"

rule make_fasta_without_tags: #fasta of raw data no trimming whatsoever
	input:
		fq=get_input_names
	output:
		fasta=PROCESS+"fasta/Full_{sample}.fa"
	run: 
		shell("samtools fasta {input} -o {output} > {output}") #takes in bam and outputs fasta, output has to be mentioned twice: 1 for the location, 2 for the option to write "both" reads to one fasta
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
	resources:
		mem_mb=5000
	shell: #added N=0 instead of default N=1
		"""
		minimap2 -y -ax map-ont --score-N 0 {input.genome} {input.fasta} | samtools sort |  samtools view -F 2304 -o {output} #added the removal of sec and suppl alignments 
		samtools index {output}
		"""

rule insertion_mapping: #conserves tags!
	input:
		bam=get_input_names,
		minimapref=PROCESS+"mapping/ref_genome_index.mmi",
		ref=PROCESS+"mapping/vector_ref_genome.fa"
	output:
		PROCESS+"mapping/Precut_{sample}_sorted.bam"
	resources:
		mem_mb=5000
	shell:
		"""
		samtools bam2fq -T '*' {input.bam}| minimap2 -y -ax map-ont {input.minimapref} - | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		"""

rule clean_postcut_by_maping_quality:
	input:
		PROCESS+"mapping/Postcut_{sample}_unfiltered_sorted.bam"
	params:
		mapq=config["MAPQ"]
	output:
		PROCESS+"mapping/Postcut_{sample}_sorted.bam"
	shell:
		"""
		samtools view -h -q {params.mapq} {input} -o {output}
		samtools index {output}
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
	run:
		shell("bedtools bamtobed -i {input.precut} > {output.precut}")
		shell("bedtools bamtobed -i {input.postcut} > {output.postcut}")  

######
######
###### fasta preparation: Cut out of blast-detected vector fragments and create new "cut-out" fasta 
######
######

rule get_cleavage_sites_for_fasta: #filters and combines matches
	input:
		#PROCESS+"blastn/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
		PROCESS+"blastn/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		filteroption=True,
		filtervalue=config["MinInsertionLength"], 
		overlap=3*FRAG # 2*FRAG #this is the distance of the start-stop that is allowed to exist to still be combined; This should not be lower than FRAG!
	output:
		cleavage=PROCESS+"blastn/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		reads=PROCESS+"blastn/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt"
	run:
		vhf.splitting_borders(input[0],params.filteroption, params.filtervalue, params.overlap, output.cleavage, output.reads)

rule split_fasta:
	input:
		breakpoints=PROCESS+"blastn/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		fasta=PROCESS+"fasta/Full_{sample}.fa"
	params:
		mode=config["splitmode"] #if each split fasta substring should be used individually, use "Separated" Join, New mode: Buffer
	output:
		fasta=PROCESS+"fasta/Cleaved_{sample}_noVector.fa",
		vector=PROCESS+"fasta/Insertion_{sample}_Vector.fa"
	run:
		vhf.split_fasta_by_borders(input.breakpoints, input.fasta, params.mode, output.fasta, output.vector)

######
######
###### Vector preparation: Fragmentation 
######
######
rule prepare_vector:
	input:
		config["vector_fasta"] #vector fasta sequence
	output: 
		fasta=PROCESS+"fasta/fragments/Forward_Backward_Vector.fa"
	run:
		vhf.reversevector(input[0], output[0])
		
rule vector_fragmentation:
	input: 
		#config["vector_fasta"] #vector fasta sequence
		PROCESS+"fasta/fragments/Forward_Backward_Vector.fa" #does not change anything so it can be removed imo
	params:
		FRAG
	output: 
		fasta=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	run:
		vhf.fragmentation_fasta(input[0], params[0], output[0])

######
######
###### BLAST Searches - Against Vector and healthy human reference
######
######

rule make_blastn_DB:
	input:
		vector_fragmented = PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
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
	run:
		shell("makeblastdb -in {input.vector_fragmented} -dbtype nucl -blastdb_version 5")

rule find_vector_BLASTn:
	input:
		fasta=PROCESS+"fasta/Full_{sample}.fa",
		dummy=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		tempdir=PROCESS+"temp_{sample}",
		vector=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:
		"""
		mkdir {params.tempdir}
		
        	blastn -query {input.fasta} -db {params.vector} -out {params.tempdir}/temp_output.blastn -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'

        	# Filter results based on bitscore > 50
        	awk '$16 > 50' {params.tempdir}/temp_output.blastn > {output}

        	# Clean up temporary files
        	rm -r {params.tempdir}
        	"""
		
#blastn vector against human genome: Which vector parts are close to human sequences so that they might raise a false positivite BLAST match
rule find_vector_BLASTn_in_humanRef:
	input:
		fasta=PROCESS+"fasta/fragments/" + str(FRAG) + "_Vector_fragments.fa"
	params:
		vector=config["blastn_db"] #full human reference
	output:
		temp(PROCESS+"blastn/humanref/"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	run:
		shell("blastn -query {input} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'")

rule hardcode_blast_header:
	input: 
		vector=PROCESS+"blastn/"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		vector=PROCESS+"blastn/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:
		"""
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.vector} > {output.vector}
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.humanref} > {output.humanref}
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
	output:
		blast=PROCESS+"blastn/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"blastn/humanref/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:
		shell("awk -F'\t' '$11>={params.threshold}' {input.blast} > {output.blast}")
		shell("awk -F'\t' '$11>={params.threshold}' {input.humanref} > {output.humanref}") 

######
######
###### Visualisation of insertions
######



rule insertion_heatmap:
	input:
		expand(PROCESS+"localization/ExactInsertions_{sample}.bed", sample=SAMPLES)
	output:
		report(FINAL+"localization/Heatmap_Insertion_Chr.png")
	run:
		vhf.plot_bed_files_as_heatmap(input, output[0])


rule insertion_length_plot:
	input:
		expand(PROCESS+"localization/ExactInsertions_{sample}.bed", sample=SAMPLES)
	output:
		report(FINAL+"localization/Insertion_length.png")
	run:

		vhf.plot_insertion_length(input, output[0])

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
	output:
		out=PROCESS+"localization/ExactInsertions_{sample}.bed",
	run:
		vhf.exact_insertion_coordinates(input.borders, input.bed, output.out)

rule collect_outputs:
	input:
		coordinates=PROCESS+"localization/ExactInsertions_{sample}.bed",
	output:
		coordinates=FINAL+"localization/ExactInsertions_{sample}.bed"
	shell:
		"""
		cp {input.coordinates} {output.coordinates}
		"""

