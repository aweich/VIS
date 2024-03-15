import os
import time
from pathlib import Path
from time import sleep


configfile: "config.yml"
SRC="./Src"
SAMPLES = expand(config["samples"]) 
PROCESS = os.path.join(config["processing_dir"],str(config["experiment"]+"/")) #intermediate and results files are stored here
FRAG=config["fragment_size"]

#local functions - path to helper fucntions needs to be added to the sys path, otherwise import won't find the file
rootpath = os.path.join(SRC)
sys.path.append(rootpath)
#print(rootpath)	
import VIS_helper_functions as vhf #functions to make snakemake pipeline leaner

#inmport rules
include: SRC+"/epigenetics.smk"

#target rule		
rule all:
	input: 
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed", sample=SAMPLES),
		expand(PROCESS+"FASTA/Postcut_Reads_with_Insertion_{sample}_Vector.fasta", sample=SAMPLES),
		expand(PROCESS+"FASTA/InsertionReads/{sample}_Clustalo/", sample=SAMPLES), #multiple sequence alignment
		expand(PROCESS+"FASTA/Conservedtags_WithTags_Full_{sample}.fa", sample=SAMPLES), #for phred
		#Visuals
		##expand(PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.gff", sample=SAMPLES),
		##expand(PROCESS+"LOCALIZATION/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		##expand(PROCESS+"FUNCTIONALGENOMICS/TF_" + str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		##expand(PROCESS+"FUNCTIONALGENOMICS/Genes_" + str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		##expand(PROCESS+"FUNCTIONALGENOMICS/Formatted_Genes_" + str(FRAG)+"_{sample}.bed", sample=SAMPLES), 
		##expand(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		##expand(PROCESS+"BLASTN/HUMANREF/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		#PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png",
		#PROCESS+"LOCALIZATION/Insertion_length.png",
		#deeper
		##expand(PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		##expand(PROCESS+"VARIANTS/BLASTN/Annotated_SNIFFLES_INS_Variant_{sample}.blastn", sample=SAMPLES),
		expand(PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_Variant_{sample}.bed",sample=SAMPLES),
		#expand(PROCESS+"FUNCTIONALGENOMICS/Genes_" + str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		##expand(PROCESS+"FASTA/Protein_Insertion_{sample}_Vector.orf", sample=SAMPLES),
		#pooling
		#expand(PROCESS+"MAPPING/POOLED/{sample}_sorted.bam", sample=SAMPLES),
		#PROCESS+"MAPPING/POOLED/Pooled_S3.bam",
		#orfs
		##expand(PROCESS+"FASTA/PROTEINBLAST/ORFs_{sample}.proteinblast", sample=SAMPLES),
		##expand(PROCESS+"FASTA/BED_{sample}_Vector.orf", sample=SAMPLES),
		#CIGAR
		##expand(PROCESS+"QC/ReadLevel_MappingQuality_{sample}.txt", sample=SAMPLES),
		##expand(PROCESS+"QC/cigar_{sample}_postcut.txt", sample=SAMPLES),
		##expand(PROCESS+"QC/1000I_cigar_FASTA_{sample}.fa", sample=SAMPLES),
		expand(PROCESS+"BLASTN/Cigar/Reads_cigar1000_VectorMatches_{sample}.fasta", sample=SAMPLES),
		expand(PROCESS+"CIGAR/Cigar_selected_reads_with_Insertions_{sample}.fastq", sample=SAMPLES),
		#expand(PROCESS+"BLASTN/Cigar/cigar_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		#MODULES
		#maybe with dmr but nor sure if its worth the hassle
		#rules from epigenetics file
		expand(PROCESS+"METHYLATION/Proximity_ExactInsertions_"+str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		expand(PROCESS+"METHYLATION/Precut_Methyl_{sample}.bed.gz", sample=SAMPLES),
		expand(PROCESS+"METHYLATION/Precut_Methyl_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"LOCALIZATION/Insertion_fasta_{sample}.bed", sample=SAMPLES),
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
rule minimap_index:
	input:
		ref=config["ref_genome"] #cut _index
	output:
		index=PROCESS+"MAPPING/ref_genome_index.mmi"
	shell:
		"minimap2 -d {output.index} {input.ref}"

rule make_FASTA_with_tags: #needed for precut path
	input:
		fq=get_input_names
	output:
		fasta=PROCESS+"FASTA/Conservedtags_WithTags_Full_{sample}.fa"
	run: 
		shell("samtools bam2fq -T '*' {input} > {output}") 

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
		genome=config["ref_genome"] #_ctrl, if we still use the modified reference genome, we will detect the reads that have no matches with blast but are still assigned to reference! -> should be zero, but lets see
	output:
		PROCESS+"MAPPING/Postcut_{sample}_sorted.bam"
	shell: #added N=0 instead of default N=1
		"""
		minimap2 -y -ax map-ont --score-N 0 {input.genome} {input.fasta} | samtools sort |  samtools view -F 2304 -o {output} #added the removal of sec and suppl alignments 
		samtools index {output}
		"""   

rule insertion_mapping: #conserves tags!
	input:
		bam=get_input_names,
		minimapref=PROCESS+"MAPPING/ref_genome_index.mmi",
		ref=config["ref_genome"]
	output:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
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
		postcut=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam" #Reads that contained an insertion before, are now marked with "_Insertion"
	output:
		postcut=PROCESS+"MAPPING/Postcut_{sample}.bed",
		precut=PROCESS+"MAPPING/Precut_{sample}.bed"
	run:
		shell("bedtools bamtobed -i {input.precut} > {output.precut}")
		shell("bedtools bamtobed -i {input.postcut} > {output.postcut}")  

rule reads_with_BLASTn_matches: #not exclusive but definitely inclusive
	input:
		refbed=PROCESS+"MAPPING/Postcut_{sample}.bed",
	output:
		loc=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed" 
	shell:
		"grep '_' {input.refbed} > {output.loc}"  
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
		mode="Buffer" #if each split FASTA substring should be used individually, use "Separated" Join, New mode: Buffer
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
		
rule orf_prediction:
	input:
		PROCESS+"FASTA/Reads_with_Insertion_{sample}_Vector.fasta"
	output:
		proteinorf=PROCESS+"FASTA/Protein_Insertion_{sample}_Vector.orf",
		fastaorf=PROCESS+"FASTA/FASTA_Insertion_{sample}_Vector.orf"
	shell: #-n = ignore nested ORFs: Should give us only the biggest right=
		"""
		./Src/ORFfinder -s 0 -in {input} -outfmt 0 -out {output.proteinorf}
		./Src/ORFfinder -s 0 -n true -in {input} -outfmt 1 -out {output.fastaorf}
		"""

rule orf_reshape: #start and stop for - strand is eínterchanged to make them plottable
	input:
		orfs=PROCESS+"FASTA/FASTA_Insertion_{sample}_Vector.orf"
	output:
		bed=PROCESS+"FASTA/BED_{sample}_Vector.orf"
	run:
		vhf.orf_reshape(input.orfs, output.bed)

###### ORF to protein BLAST
rule protein_blast:
	input:
		PROCESS+"FASTA/Protein_Insertion_{sample}_Vector.orf"
	params: 
		db=config["proteindb"]
	output:
		temp(PROCESS+"FASTA/PROTEINBLAST/Protein_Insertion_NoHeaderBlast_{sample}.temp")
	shell:
		"blastp -query {input} -db {params.db} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs evalue bitscore' -taxids 9606"
		
				
rule protein_blast_header:		
	input: 
		PROCESS+"FASTA/PROTEINBLAST/Protein_Insertion_NoHeaderBlast_{sample}.temp"
	output:
		PROCESS+"FASTA/PROTEINBLAST/ORFs_{sample}.proteinblast"
	shell:	
		"echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tEvalue\tBitscore' | cat - {input} > {output}"



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
###### Quality control and normalization
######
######
rule nanoplot:
	input:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
	output:
		PROCESS+"QC/{sample}/Non_weightedHistogramReadlength.png"
	params:
		outdir=PROCESS+"QC/{sample}/"
	shell: 
		"""
		NanoPlot --bam {input} -o {params.outdir}
		touch {output}
		"""

rule mapping_qc:
	input:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam" #before cut out, but after removal of secondary and supplementary alignments
	output:
		short=PROCESS+"QC/short_{sample}.qc",
		long=PROCESS+"QC/long_{sample}.qc"
	shell:
		"""
		samtools flagstats {input} > {output.short}
		samtools stats {input} > {output.long}    
		"""
		
rule nReads_for_insertion_count:
	input:
		get_input_names
	output:
		temp(PROCESS+"QC/Number_of_Bases_{sample}.normalisation")
	shell:
		"gatk CountBases -I {input} > {output}"				

rule normalisation_for_insertion_count:
	input:
		insertions=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt",
		number_of_bases=PROCESS+"QC/Number_of_Bases_{sample}.normalisation",
		fasta=PROCESS+"FASTA/Full_{sample}.fa" #for N50
	params:
		scale=3000000000 #IPM #IPG
	output:
		PROCESS+"QC/Normalisation_IPG_{sample}.txt"
	run:
		vhf.insertion_normalisation(input.insertions, input.number_of_bases, params[0], input.fasta, output[0])

rule mapping_for_cigar:
	input:
		bam=get_input_names,
		ref=config["ref_genome_ctrl"]
	output:
		PROCESS+"QC/HealthyRef_cigar_{sample}_sorted.bam"
	shell:
		"""
		samtools bam2fq -T '*' {input.bam}| minimap2 -y -ax map-ont {input.ref} - | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		"""

rule cigar_strings:
	input:
		pre=PROCESS+"QC/HealthyRef_cigar_{sample}_sorted.bam", #PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
		post=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam",
		matches=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt",
	output:
		pre=PROCESS+"QC/cigar_{sample}_precut.txt",
		names=PROCESS+"QC/cigar_{sample}_precut_names.txt",
		post=PROCESS+"QC/cigar_{sample}_postcut.txt",
	shell:
		"""
		samtools view {input.pre} | cut -f 1,3,4,6 | grep '[0-9]\{{4\}}I' > {output.pre}
		samtools view {input.pre} | cut -f 1,3,4,6 | grep '[0-9]\{{4\}}I' | cut -f 1 > {output.names}
		samtools view {input.post} | grep -f {input.matches} | cut -f 1,3,4,6 > {output.post}
		"""
		
rule cigar_fasta:
	input:
		fasta=PROCESS+"FASTA/Full_{sample}.fa",
		matches=PROCESS+"QC/cigar_{sample}_precut_names.txt",
	output:
		PROCESS+"QC/1000I_cigar_FASTA_{sample}.fa",
	run: 
		shell("seqkit grep -r -f {input.matches} {input.fasta} -o {output}")

rule cigar_fasta_blast:
	input:
		fasta=PROCESS+"QC/1000I_cigar_FASTA_{sample}.fa",
		dummy=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		#vector=config["blastn_db"] #vector db
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		PROCESS+"BLASTN/Cigar/cigar"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:
		shell("blastn -query {input.fasta} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'") 

rule cigar_hardcode_blast_header:		
	input: 
		PROCESS+"BLASTN/Cigar/cigar"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Cigar/cigar_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:	
		shell("echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input} > {output}")	
		
rule cigar_reads_and_blast_insertion: #just checks if read id is in the cleavage sites file
	input:
		fasta1=PROCESS+"QC/1000I_cigar_FASTA_{sample}.fa",
		breakpoints=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Cigar/Reads_cigar1000_VectorMatches_{sample}.fasta",
	run:
		vhf.reads_with_insertions(input.fasta1, input.breakpoints, output[0])		

rule fasta_to_fastq: #first line extracts the fastq entries, second line separates them to make the input for fastaqc easier
	input:
		fastq=PROCESS+"FASTA/Conservedtags_WithTags_Full_{sample}.fa",
		fasta=PROCESS+"BLASTN/Cigar/Reads_cigar1000_VectorMatches_{sample}.fasta"
	output:
		full=PROCESS+"CIGAR/Cigar_selected_reads_with_Insertions_{sample}.fastq"
	shell:
		'''
		grep '>' {input.fasta} | cut -c 2- | seqkit grep -f - {input.fastq} -o {output.full}
		csplit -z {output.full} $(awk '/^@/ {{print NR }}' {output.full}) '/^@/' {{*}} --prefix={output.full}
		'''
			
#### read-level stats
rule read_level_stats:
	input:
		bam=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam",
		matches=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt", 
	output:
		PROCESS+"QC/ReadLevel_MappingQuality_{sample}.txt"
	shell:
		"""
		samtools view -h {input.bam} | grep -e '^@' -f {input.matches} | samtools stats | grep '^SN' | cut -f 2- > {output}
		"""


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

rule chromosome_read_plots: 
	input:
		bed=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		bam=PROCESS+"MAPPING/NoVectorAlignments_Postcut_{sample}_sorted.bam",
		H3K4Me1=config["ucsc_H3K4Me1"],
		H3K4Me3=config["ucsc_H3K4Me3"],
		H3K27Ac=config["ucsc_H3K27Ac"],
		gtf=config["ucsc_Genes_gtf"],
		TF=config["ucsc_TF"]#PROCESS+"FUNCTIONALGENOMICS/TF_" + str(FRAG)+"_{sample}.bed",		
	output:
		outpath=directory(PROCESS+"LOCALIZATION/PLOTS/" + str(FRAG)+"_{sample}")
	params:
		buffer=50000
	shell: 
		r"""
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -iH3K4Me1 {input.H3K4Me1} -iH3K4Me3 {input.H3K4Me3} -iH3K27Ac {input.H3K27Ac} -igtf {input.gtf} -iTF {input.TF} -buffer {params.buffer} -o {output.outpath}   
		"""	
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
		expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed", sample=SAMPLES)
	output:
		PROCESS+"LOCALIZATION/Insertion_length.png"
	run:
		vhf.plot_insertion_length(input, output[0])
######
######
###### Variant callers and overlap check of Insertions with BLAST matches
######
######

rule variant_sniffles:
	input:
		bam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
		genome=config["ref_genome_ctrl"]
	output:
		PROCESS+"VARIANTS/SNIFFLES/Variant_{sample}.vcf"
	shell:
		"sniffles -i {input.bam} --reference {input.genome} --output-rnames -v {output}"

rule svim_variants:
	input:
		bam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
		genome=config["ref_genome_ctrl"]
	output:
		PROCESS+"VARIANTS/SVIM_{sample}/variants.vcf", 
	params:
		outdir=PROCESS+"VARIANTS/SVIM_{sample}/"
	shell:
		"""
		svim alignment {params.outdir} {input.bam} {input.genome} --types INS
		touch {output}
		"""

rule nanoVar:
	input:
		ref=config["ref_genome"],
		bam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
	output:
		PROCESS+"VARIANTS/NanoVar_{sample}/Precut_{sample}_sorted.nanovar.pass.vcf"
	params:
		outdir=PROCESS+"VARIANTS/NanoVar_{sample}/"
	shell:
		"""
		nanovar -f hg38 {input.bam} {input.ref} {params.outdir}
		touch {output}
		"""

#reshapes the variants in a more usable format for downstream analysis
rule reshape_variants:
	input:
		nanovar=PROCESS+"VARIANTS/NanoVar_{sample}/Precut_{sample}_sorted.nanovar.pass.vcf",
		#svim=PROCESS+"VARIANTS/SVIM_{sample}/variants.vcf", 
		#sniffles=PROCESS+"VARIANTS/SNIFFLES/Variant_{sample}.vcf"
	output:
		nanovar=PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_Variant_{sample}.bed",
		#svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.bed", 
		#sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.bed"
	run:
		#shell("convert2bed --input=vcf  < {input.nanovar} > {output.nanovar}")
		shell("vcf2bed --do-not-sort < {input.nanovar} > {output.nanovar}") #do not sort is important to generate the file but results in error downstream :/
		#shell("convert2bed --input=vcf --insertions < {input.sniffles} > {output.sniffles}")
		#shell("convert2bed --input=vcf --insertions < {input.svim} > {output.svim}")
		

#extraction of reads with BLAS matches from the variant callers

rule reads_with_BLAST_from_callers:
	input:
		#nanovar=PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_INS_Variant_{sample}.bed",
		svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.bed",
		sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.bed",
		#matches=PROCESS+"BLASTN/Reads_with_VectorMatches_{sample}.blastn"
		matches=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
	output:
		#nanovar=PROCESS+"VARIANTS/Reads_with_BLAST_Nanovar_INS_Variant_{sample}.bed",
		svim=PROCESS+"VARIANTS/Reads_with_BLAST_SVIM_INS_Variant_{sample}.bed",
		sniffles=PROCESS+"VARIANTS/Reads_with_BLAST_SNIFFLES_INS_Variant_{sample}.bed"
	run:
		#shell("bedtools intersect -a {input.matches} -b {input.nanovar} > {output.nanovar}")
		shell("bedtools intersect -a {input.matches} -b {input.svim} > {output.svim}")
		shell("bedtools intersect -a {input.matches} -b {input.sniffles} > {output.sniffles}")

rule variants_with_sequences_fasta:
	input:
		svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.bed", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.bed"
	output:
		svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.fasta", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.fasta"
	run: 
		vhf.variant_bed_to_fasta(input.svim, output.svim)
		vhf.variant_bed_to_fasta(input.sniffles, output.sniffles)

rule blast_vector_against_insertions:
	input:
		svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.fasta", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.fasta",
		dummy=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		svim=PROCESS+"VARIANTS/BLASTN/SVIM_INS_Variant_{sample}.blastn", 
		sniffles=PROCESS+"VARIANTS/BLASTN/SNIFFLES_INS_Variant_{sample}.blastn"
	run:
		shell("blastn -query {input.svim} -db {params.vector} -out {output.svim} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs'")
		shell("blastn -query {input.sniffles} -db {params.vector} -out {output.sniffles} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs'")

rule variants_hardcode_blast_header:		
	input: 
		svim=PROCESS+"VARIANTS/BLASTN/SVIM_INS_Variant_{sample}.blastn", 
		sniffles=PROCESS+"VARIANTS/BLASTN/SNIFFLES_INS_Variant_{sample}.blastn"
	output:
		svim=PROCESS+"VARIANTS/BLASTN/Annotated_SVIM_INS_Variant_{sample}.blastn", 
		sniffles=PROCESS+"VARIANTS/BLASTN/Annotated_SNIFFLES_INS_Variant_{sample}.blastn"
	run:	
		shell("echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input.svim} > {output.svim}")
		shell("echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input.sniffles} > {output.sniffles}")
		

######
######
###### WIP rules
######
######
#exact coordinates of the matching fragments

rule exact_insertion_coordinates:
	input:
		#bed=PROCESS+"MAPPING/Postcut_{sample}.bed", #full bed, maybe a inbetween step can be replaced!
		bed=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed", 
		borders=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn" #some entries with cleavage site won't be in the output, if they were not mapped to the genome in the postcut sample; but a insertion withput genomic coordinates does not help us anyway
	output:
		out=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		out2=PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed" #only to get FASTA sequence
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

rule proximity_generator:
	input:
		bed=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
	params:
		offsets=[0,1000,5000,10000,500000]
	output:
		out=PROCESS+"LOCALIZATION/Proximity_to_ExactInsertions_{sample}.bed"
	run:
		vhf.proximity_generator_for_bed_file(input.bed, output.out, params.offsets)

rule TF_binding: 
	input:
		bed=PROCESS+"LOCALIZATION/Proximity_to_ExactInsertions_{sample}.bed",
		TF=config["ucsc_TF"]
		
	output:
		PROCESS+"FUNCTIONALGENOMICS/TF_" + str(FRAG)+"_{sample}.bed"
	shell: 
		"bedtools intersect -wa -a {input.bed} -wb -b {input.TF} | cut -f1,2,3,4,5,9 | sort -k6 > {output}"
		
rule Genes_in_prox: 
	input:
		bed=PROCESS+"LOCALIZATION/Proximity_to_ExactInsertions_{sample}.bed",
		TF=config["ucsc_Genes"]
		
	output:
		PROCESS+"FUNCTIONALGENOMICS/Genes_" + str(FRAG)+"_{sample}.bed"
	shell: 
		"bedtools intersect -wa -a {input.bed} -wb -b {input.TF} | cut -f1,2,3,4,5,9 | sort -k6 > {output}"	

rule reshape_functional_tables:
	input:
		TF=PROCESS+"FUNCTIONALGENOMICS/TF_" + str(FRAG)+"_{sample}.bed",
		Genes=PROCESS+"FUNCTIONALGENOMICS/Genes_" + str(FRAG)+"_{sample}.bed"
	output:
		TF=PROCESS+"FUNCTIONALGENOMICS/Formatted_TF_" + str(FRAG)+"_{sample}.bed",
		Genes=PROCESS+"FUNCTIONALGENOMICS/Formatted_Genes_" + str(FRAG)+"_{sample}.bed"
	run:
		vhf.reshape_functional_tables(input.TF, output.TF)
		vhf.reshape_functional_tables(input.Genes, output.Genes) 	
'''
#Modules
module epigenetics:
    snakefile:
        SRC+"Snakefile_epigenetics"
    config: SRC+config["config_epigenetics"]
    prefix: "epigenetics"

use rule * from epigenetics as epigenetics_*
'''
