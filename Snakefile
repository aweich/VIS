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

#inmport rules
include: config["functional_genomics"]
include: config["quality_control"]
include: config["epigenetics"]
include: config["variants"]

#target rule        
rule all:
	input: 
		#main output
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}_estimated_full_coordinates.bed", sample=SAMPLES),
		PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png",
		##PROCESS+"LOCALIZATION/Insertion_length.png",
		##expand(PROCESS+"FUNCTIONALGENOMICS/LOCALIZATION/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		expand(PROCESS+"BLASTN/PLOTS/Longest_Interval_{sample}", sample=SAMPLES),
		#expand(PROCESS+"FASTA/InsertionReads/{sample}_Clustalo/", sample=SAMPLES), #multiple sequence alignment
		##expand(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		##expand(PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.gff", sample=SAMPLES),
		#expand(PROCESS+"BLASTN/HUMANREF/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		#pooling
		#expand(PROCESS+"MAPPING/POOLED/{sample}_sorted.bam", sample=SAMPLES),
		#PROCESS+"MAPPING/POOLED/Pooled_S3.bam",
		#MODULES
		###rules to generate functional genomics output
		##expand(PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png", sample=SAMPLES),
		#expand(PROCESS+"FUNCTIONALGENOMICS/ORF/PROTEINBLAST/ORFs_{sample}.proteinblast", sample=SAMPLES),
		###rules to generate qc output
		#expand(PROCESS+"QC/Nanoplot/{sample}/Non_weightedHistogramReadlength.png", sample=SAMPLES),
		##expand(PROCESS+"QC/Normalisation_IPG_{sample}.txt", sample=SAMPLES),
		#expand(PROCESS+"QC/Coverage/Genomecoverage_{sample}.bed", sample=SAMPLES),
		#cigar output (as part of qc)
		#expand(PROCESS+"QC/CIGAR/Reads_with_longInsertions_and_vector_{sample}.fastq", sample=SAMPLES),
		###rules to generate variant output
		#expand(PROCESS+"VARIANTS/BCFTOOLS/Variant_{sample}.vcf", sample=SAMPLES),
		#expand(PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_Variant_{sample}.bed", sample=SAMPLES), #all three callers share this rule
		###rules to generate epigenetics output
		#expand(PROCESS+"METHYLATION/Proximity_ExactInsertions_"+str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		##expand(PROCESS+"METHYLATION/Precut_Methyl_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/PLOTS/InsertionRead_{sample}/",sample=SAMPLES),

		
#actual filenames
def get_input_names(wildcards):
	return config["samples"][wildcards.sample]

#BAM Operations:
#Add rule to remove supplementary and secondary alignments
#samtools view -F 2304 -bo filtered.bam original.bam (via picard: supplementary alignment, not primary alignment): Added to prepare BAM ruleF

print("Results in ...")
print(PROCESS)

######
######
###### Insertion BAM: Convert to FASTA 
######
######
#bam pooling
'''
rule prepare_BAM: #sorts, index, and removes supllementary and secondary alignments
	input:
		get_input_names
	output:
		PROCESS+"MAPPING/POOLED/{sample}_sorted.bam"
	shell:
		"""
		samtools sort {input} > {output} 
		samtools index {output}
		"""

rule pool_BAM:
	input:
		expand(PROCESS+"MAPPING/POOLED/{sample}_sorted.bam", sample=SAMPLES)
	output:
		PROCESS+"MAPPING/POOLED/Pooled_S3.bam"
	shell:
		"""
		samtools merge -o {output} {input}
		"""
'''

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
		genome=PROCESS+"MAPPING/vector_ref_genome.fa" #_ctrl, if we still use the modified reference genome, we will detect the reads that have no matches with blast but are still assigned to reference
	output:
		PROCESS+"MAPPING/Postcut_{sample}_sorted.bam"
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

rule reads_with_insertion: #just checks if read id is in the cleavage sites file
	input:
		fasta1=PROCESS+"FASTA/Full_{sample}.fa",
		fasta2=PROCESS+"FASTA/Cleaved_{sample}_noVector.fa",
		breakpoints=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"FASTA/Reads_with_Insertion_{sample}_Vector.fasta",
		PROCESS+"FASTA/Postcut_Reads_with_Insertion_{sample}_Vector.fasta"
	run:
		vhf.reads_with_insertions(input.fasta1, input.breakpoints, output[0])
		vhf.reads_with_insertions(input.fasta2, input.breakpoints, output[1])

rule summary_reads_with_insertion: #just checks if read id is in the cleavage sites file
	input:
		fullvector=config["vector_fasta"],
		insertedvector=PROCESS+"FASTA/Insertion_{sample}_Vector.fa",
		precut_reads=PROCESS+"FASTA/Reads_with_Insertion_{sample}_Vector.fasta",
		postcut_reads=PROCESS+"FASTA/Postcut_Reads_with_Insertion_{sample}_Vector.fasta",
		breakpoints=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		directory(PROCESS+"FASTA/InsertionReads/{sample}/")
	run:
		shell("mkdir {output}")
		vhf.summary_reads_with_insertions(input.fullvector, input.insertedvector, input.precut_reads, input.postcut_reads, input.breakpoints, output[0])

rule clustal_omega:
	input:
		PROCESS+"FASTA/InsertionReads/{sample}/"
	params:
		clustalo='clustalo' #depends on how it is stored in the bashrc
	output:
		directory(PROCESS+"FASTA/InsertionReads/{sample}_Clustalo/")
	shell:
		"""
		mkdir {output}
		for file in {input}/*; do
		echo $file
		clustalo -i $file -o "{output}/$(basename "${{file}}")_clustalo.fasta"
		done
		"""     
######
######
###### Vector preparation: Fragmentation 
######
######
rule reverse_vector:
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
		#vector=config["blastn_db"] #vector db
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:
		shell("blastn -query {input.fasta} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'") 

rule hardcode_blast_header:     
	input: 
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:    
		shell("echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input} > {output}")
		
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
rule hardcode_blast_header_humanRef:        
	input: 
		PROCESS+"BLASTN/HUMANREF/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:    
		shell("echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevaöue\tbitscore' | cat - {input} > {output}")

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
######
######
###### Visualisations of intermediate results
######
###### the follwoing 2 rules are needed to remove the pseudo-alignment of the reads to the modified ref genome s this would otherwise break the plotting function
'''
rule get_reads_of_supplementary_matches:
	input:
		post=PROCESS+"MAPPING/Postcut_{sample}.bed",
		pre=PROCESS+"MAPPING/Precut_{sample}.bed"
	output:
		post=PROCESS+"MAPPING/Postcut_{sample}_lostvector.reads",
		pre=PROCESS+"MAPPING/Precut_{sample}_lostvector.reads"
	run:
		shell("awk '$1 ~ /CAR/ {{print $4}}' {input.post} > {output.post}") #CD19 with CD/CAR, CD123 with V0. Uff that sucks
		shell("awk '$1 ~ /CAR/ {{print $4}}' {input.pre} > {output.pre}")

rule remove_vector_alignments:
	input:
		prebam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
		prelost=PROCESS+"MAPPING/Precut_{sample}_lostvector.reads",
		postbam=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam",
		postlost=PROCESS+"MAPPING/Postcut_{sample}_lostvector.reads"
	output:
		pre=PROCESS+"MAPPING/NoVectorAlignments_Precut_{sample}_sorted.bam",
		post=PROCESS+"MAPPING/NoVectorAlignments_Postcut_{sample}_sorted.bam"
	shell: #V0 for CD123, CAR for CD19
		"""
		samtools view -h {input.postbam} | grep -v 'CAR' | samtools view -bS -o {output.post} -
		samtools index {output.post}
		samtools view -h {input.prebam} | grep -v 'CAR' | samtools view -bS -o {output.pre} -
		samtools index {output.pre}
		"""
'''
rule fragmentation_distribution_plots:
	input:
		PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		PROCESS+"BLASTN/HUMANREF/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
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
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed", sample=SAMPLES)
	output:
		PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png"
	run:
		vhf.plot_bed_files_as_heatmap(input, output[0])


rule insertion_length_plot:
	input:
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}_estimated_full_coordinates.bed", sample=SAMPLES)
	output:
		PROCESS+"LOCALIZATION/Insertion_length.png"
	run:

		vhf.plot_insertion_length(input, output[0])

rule detailed_insertion_length_plot:
	input:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
	params: 
		buffer=3*FRAG
	output:
		outpath=directory(PROCESS+"BLASTN/PLOTS/Longest_Interval_{sample}/")
	run:
		shell("mkdir {output.outpath}")
		vhf.find_and_plot_longest_blast_interval(input[0], params[0], output[0])

######
######
###### WIP rules
######
######
#exact coordinates of the matching fragments

rule exact_insertion_coordinates:
	input:
		bed=PROCESS+"MAPPING/Postcut_{sample}.bed", #full bed, maybe a inbetween step can be replaced!
		borders=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn" #some entries with cleavage site won't be in the output, if they were not mapped to the genome in the postcut sample; but a insertion withput genomic coordinates does not help us anyway
	output:
		out=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		out2=PROCESS+"LOCALIZATION/ExactInsertions_{sample}_estimated_full_coordinates.bed" #only to get FASTA sequence
	run:
		vhf.exact_insertion_coordinates2(input.borders, input.bed, output.out, output.out2)

#filter for BLAST matches
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

'''
#Control of the workflow: How do the genomic coordinates of my two BAMs differ for the insertion vectors! # CUrrently it looks like all the insertion reads do not really differ pre and post cut -> either the cut out is not sensitive ennough (try with split reads) or the vector doesn't really alter the alignment!
rule check_mapping_pre_and_postcut:
	input:
		postcutbed=PROCESS+"MAPPING/Postcut_{sample}.bed",
		precutbed=PROCESS+"MAPPING/Precut_{sample}.bed",
	output:
		notinpostcut=PROCESS+"MAPPING/Not_in_postcut_{sample}.bed",
		notinprecut=PROCESS+"MAPPING/Not_in_precut_{sample}.bed"
	shell:
		"""
			bedtools intersect -v -a {input.precutbed} -b {input.postcutbed} > {output.notinpostcut}
			bedtools intersect -v -a {input.postcutbed} -b {input.precutbed} > {output.notinprecut}
			"""
rule reads_with_matches:
	input:
		readswithmatches=PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Reads_with_VectorMatches_{sample}.blastn"
	shell:
		"cat {input.readswithmatches} | cut -f 1 > {output}"    

rule pre_post_cut_and_reads_with_matches:
	input:
		#notinpostcut=PROCESS+"MAPPING/Not_in_postcut_{sample}.bed",
		#notinprecut=PROCESS+"MAPPING/Not_in_precut_{sample}.bed",
		notinpostcut=PROCESS+"MAPPING/Postcut_{sample}.bed",
		notinprecut=PROCESS+"MAPPING/Precut_{sample}.bed",
		matches=PROCESS+"BLASTN/Reads_with_VectorMatches_{sample}.blastn"
	output:
		readsnotinpostcut=PROCESS+"EVAL/Matches_different_in_postcut_{sample}.bed",
		readsnotinprecut=PROCESS+"EVAL/Matches_different_in_precut_{sample}.bed"
	run: 
		shell("if ! grep -F -f {input.matches} {input.notinpostcut}; then echo not found; fi > {output.readsnotinpostcut}")
		shell("if ! grep -F -f {input.matches} {input.notinprecut}; then echo not found; fi > {output.readsnotinprecut}")

rule pre_post_cut_and_check_FASTA_changes:
	input:
		postcut=PROCESS+"FASTA/Cleaved_{sample}_noVector.fa",
		precut=PROCESS+"FASTA/Full_{sample}.fa",
		matches=PROCESS+"BLASTN/Reads_with_VectorMatches_{sample}.blastn"
	output:
		readsnotinpostcut=PROCESS+"EVAL/Matches_different_in_postcut_FASTA_{sample}.fa",
		readsnotinprecut=PROCESS+"EVAL/Matches_different_in_precut_FASTA_{sample}.fa"
	run: 
		shell("seqkit grep -r -f {input.matches} {input.postcut} -o {output.readsnotinpostcut}")
		shell("seqkit grep -r -f {input.matches} {input.precut} -o {output.readsnotinprecut}")
'''
