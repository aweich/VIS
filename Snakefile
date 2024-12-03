import os
import time
from pathlib import Path
from time import sleep


configfile: "config.yml"
SRC=config["source_dir"]
SAMPLES = expand(config["samples"]) 
OUTDIR = os.path.join(config["processing_dir"],str(config["experiment"]+"/"))
PROCESS = OUTDIR+"Intermediate/" #intermediate files are stored here
FINAL = OUTDIR+"Final/" # final files are stored here
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
		#main output
		OUTDIR+"config_settings.yml",
		FINAL+"LOCALIZATION/Heatmap_Insertion_Chr.png",
		FINAL+"LOCALIZATION/Insertion_length.png",
		FINAL+"QC/multiqc_report.html",
		expand(FINAL+"QC/Fragmentation/Insertions/insertions_" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(FINAL+"QC/Fragmentation/Longest_Interval/{sample}/", sample=SAMPLES),
		expand(PROCESS+"FUNCTIONALGENOMICS/LOCALIZATION/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(FINAL+"LOCALIZATION/ExactInsertions_{sample}.bed", sample=SAMPLES),
		expand(PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png", sample=SAMPLES),
		#expand(PROCESS+"FASTA/InsertionReads/{sample}_Clustalo/", sample=SAMPLES), #multiple sequence alignment
		#expand(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		##expand(PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.gff", sample=SAMPLES),
		#pooling
		#expand(PROCESS+"MAPPING/POOLED/{sample}_sorted.bam", sample=SAMPLES),
		#PROCESS+"MAPPING/POOLED/Pooled_S3.bam",
		#MODULES
		###rules to generate functional genomics output
		#expand(PROCESS+"FUNCTIONALGENOMICS/Functional_distances_to_Insertions_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png", sample=SAMPLES),
		###rules to generate qc output
		#expand(PROCESS+"QC/Nanoplot/{sample}/NanoStats.txt", sample=SAMPLES),
		##expand(PROCESS+"QC/Normalisation_IPG_{sample}.txt", sample=SAMPLES),
		#expand(PROCESS+"QC/Coverage/Genomecoverage_{sample}.bed", sample=SAMPLES),
		#cigar output (as part of qc)
		#expand(PROCESS+"QC/CIGAR/Reads_with_longInsertions_and_vector_{sample}.fastq", sample=SAMPLES),
		#PROCESS+"QC/fastqc/multiqc_report.html",
		#expand(PROCESS + "QC/readlevel_{sample}/", sample=SAMPLES),
		#QC
		expand(PROCESS+"QC/MAPQ/{sample}_mapq_heatmap_image.png", sample=SAMPLES),
		###rules to generate variant output
		#expand(PROCESS+"VARIANTS/BCFTOOLS/Variant_{sample}.vcf", sample=SAMPLES),
		#expand(PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_Variant_{sample}.bed", sample=SAMPLES), #all three callers share this rule
		###rules to generate epigenetics output
		#expand(PROCESS+"METHYLATION/Proximity_ExactInsertions_"+str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		##expand(PROCESS+"METHYLATION/Precut_Methyl_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/PLOTS/InsertionRead_{sample}/",sample=SAMPLES),
		#malignancy score


		
#actual filenames
def get_input_names(wildcards):
	return config["samples"][wildcards.sample]


#print("Results in ...")
#print(PROCESS)

######
######
###### Insertion BAM: Convert to FASTA 
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
		PROCESS+"MAPPING/vector_ref_genome.fa"
	shell:
		"cat {input.ref} {input.vector} > {output}"
		

rule minimap_index:
	input:
		ref=PROCESS+"MAPPING/vector_ref_genome.fa"
	output:
		index=temp(PROCESS+"MAPPING/ref_genome_index.mmi")
	resources:
		mem_mb=5000
	shell:
		"minimap2 -d {output.index} {input.ref}"

rule make_FASTA_without_tags: #fasta of raw data no trimming whatsoever
	input:
		fq=get_input_names
	output:
		fasta=PROCESS+"FASTA/Full_{sample}.fa"
	run: 
		shell("samtools fasta {input} -o {output} > {output}") #takes in bam and outputs fasta, output has to be mentioned twice: 1 for the location, 2 for the option to write "both" reads to one fasta
######
######
###### "Clean" BAM: Cut-out FASTA to BAM via Mapping to reference 
######
######
		
#rule to create the BAM files with the non-insertion reads and the splitted read fragments
rule Non_insertion_mapping: #mapping against the unaltered referenc egenome
	input:
		fasta=PROCESS+"FASTA/Cleaved_{sample}_noVector.fa",
		genome=PROCESS+"MAPPING/vector_ref_genome.fa"
	output:
		PROCESS+"MAPPING/Postcut_{sample}_unfiltered_sorted.bam"
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
		minimapref=PROCESS+"MAPPING/ref_genome_index.mmi",
		ref=PROCESS+"MAPPING/vector_ref_genome.fa"
	output:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
	resources:
		mem_mb=5000
	shell:
		"""
		samtools bam2fq -T '*' {input.bam}| minimap2 -y -ax map-ont {input.minimapref} - | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		"""

rule clean_postcut_by_maping_quality:
	input:
		PROCESS+"MAPPING/Postcut_{sample}_unfiltered_sorted.bam"
	params:
		mapq=config["MAPQ"]
	output:
		PROCESS+"MAPPING/Postcut_{sample}_sorted.bam"
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
		precut=PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
		postcut=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam" #Reads that contained an insertion before, are now marked with "_Buffer/Insertion/0,1,2..."
	output:
		postcut=PROCESS+"MAPPING/Postcut_{sample}.bed",
		precut=PROCESS+"MAPPING/Precut_{sample}.bed"
	run:
		shell("bedtools bamtobed -i {input.precut} > {output.precut}")
		shell("bedtools bamtobed -i {input.postcut} > {output.postcut}")  
 
######
######
###### FASTA preparation: Cut out of blast-detected vector fragments and create new "cut-out" FASTA 
######
######

rule get_cleavage_sites_for_fasta: #filters and combines matches
	input:
		#PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
		PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		filteroption=True,
		filtervalue=config["MinInsertionLength"], 
		overlap=3*FRAG # 2*FRAG #this is the distance of the start-stop that is allowed to exist to still be combined; This should not be lower than FRAG!
	output:
		cleavage=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		reads=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt"
	run:
		vhf.splitting_borders(input[0],params.filteroption, params.filtervalue, params.overlap, output.cleavage, output.reads)

rule split_fasta:
	input:
		breakpoints=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		fasta=PROCESS+"FASTA/Full_{sample}.fa"
	params:
		mode=config["splitmode"] #if each split FASTA substring should be used individually, use "Separated" Join, New mode: Buffer
	output:
		fasta=PROCESS+"FASTA/Cleaved_{sample}_noVector.fa",
		vector=PROCESS+"FASTA/Insertion_{sample}_Vector.fa"
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
		fasta=PROCESS+"FASTA/Fragments/Forward_Backward_Vector.fa"
	run:
		vhf.reversevector(input[0], output[0])
		
rule vector_fragmentation:
	input: 
		#config["vector_fasta"] #vector fasta sequence
		PROCESS+"FASTA/Fragments/Forward_Backward_Vector.fa" #does not change anything so it can be removed imo
	params:
		FRAG
	output: 
		fasta=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	run:
		vhf.fragmentation_fasta(input[0], params[0], output[0])

######
######
###### BLAST Searches - Against Vector and healthy human reference
######
######

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
		fasta=PROCESS+"FASTA/Full_{sample}.fa",
		dummy=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		tempdir=PROCESS+"temp_{sample}",
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:
		"""
		mkdir {params.tempdir}
		
        	blastn -query {input.fasta} -db {params.vector} -out {params.tempdir}/temp_output.blastn -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'

        	# Filter results based on bitscore > 50
        	awk '$16 > 50' {params.tempdir}/temp_output.blastn > {output}

        	# Clean up temporary files
        	rm -r {params.tempdir}
        	"""
		
#BLASTN vector against human genome: Which vector parts are close to human sequences so that they might raise a false positivite BLAST match
rule find_vector_BLASTn_in_humanRef:
	input:
		fasta=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	params:
		vector=config["blastn_db"] #full human reference
	output:
		temp(PROCESS+"BLASTN/HUMANREF/"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	run:
		shell("blastn -query {input} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'")

rule hardcode_blast_header:        
	input: 
		vector=PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"BLASTN/HUMANREF/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		vector=PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	shell:    
		"""
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.vector} > {output.vector}
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.humanref} > {output.humanref}
		"""
		
rule blast_to_gff:
	input: 
		ref=PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		vector=PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		vector=PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.gff",
		ref=PROCESS+"BLASTN/HUMANREF/"+str(FRAG)+"_VectorMatches_{sample}.gff"  
	run:
		vhf.blast2gff(input.ref, output.ref)
		vhf.blast2gff(input.vector, output.vector)  

#filter BLAST matches by min length
rule extract_by_length:
	input:
		blast=PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		threshold=config["MinLength"]
	output:
		blast=PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		humanref=PROCESS+"BLASTN/HUMANREF/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:
		shell("awk -F'\t' '$11>={params.threshold}' {input.blast} > {output.blast}")
		shell("awk -F'\t' '$11>={params.threshold}' {input.humanref} > {output.humanref}") 

######
######
###### Visualisation of insertions
######



rule insertion_heatmap:
	input:
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed", sample=SAMPLES)
	output:
		FINAL+"LOCALIZATION/Heatmap_Insertion_Chr.png"
	run:
		vhf.plot_bed_files_as_heatmap(input, output[0])


rule insertion_length_plot:
	input:
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed", sample=SAMPLES)
	output:
		FINAL+"LOCALIZATION/Insertion_length.png"
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
		bed=PROCESS+"MAPPING/Postcut_{sample}.bed",
		borders=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		out=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
	run:
		vhf.exact_insertion_coordinates(input.borders, input.bed, output.out)

rule collect_outputs:
	input:
		coordinates=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
	output:
		coordinates=FINAL+"LOCALIZATION/ExactInsertions_{sample}.bed"
	shell:
		"""
		cp {input.coordinates} {output.coordinates}
		"""

