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
		#expand(PROCESS+"FASTA/Full_{sample}.fa", sample=SAMPLES),
		#expand(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"),
		#expand(PROCESS+"FASTA/Cleaved_{sample}_noVector.fa", sample=SAMPLES),
		#expand(PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa{ext}", ext=[ ".ndb",".nhr",".nin",".not",".nsq",".ntf",".nto"]), 
		expand(PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		expand(PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		#expand(PROCESS+"MAPPING/CutOut_{sample}_sorted.bam", sample=SAMPLES),
		#expand(PROCESS+"MAPPING/BasicMapping_{sample}.qc", sample=SAMPLES),
		#expand(PROCESS+"MAPPING/Normalisation_IPHM_{sample}.txt", sample=SAMPLES),
		#expand(PROCESS+"MAPPING/BasicMapping_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"FASTA/Conservedtags_WithTags_Full_{sample}.fa" , sample=SAMPLES),
		#exact coordinates
		#expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"LOCALIZATION/Reshaped_GenomicLocation_"+str(FRAG)+"_{sample}.bed", sample=SAMPLES),
		#Methylation
		#expand(PROCESS+"METHYLATION/Methyl_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/Insertion_fasta_proximity_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/MeanMods_Proximity_{sample}.bed", sample=SAMPLES),
		#Differential Methylation
		#PROCESS+"LOCALIZATION/ExactInsertions_combined.bed",
		#PROCESS+"METHYLATION/Insertion_fasta_proximity_combined.bed",
		#Visuals
		#expand(PROCESS+"LOCALIZATION/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		#expand(PROCESS+"BLASTN/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		#expand(PROCESS+"BLASTN/HUMANREF/PLOTS/" + str(FRAG)+"_{sample}", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/Heatmap_MeanMods_Proximity_{sample}.png", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/All_Insertions/Heatmap_MeanMods_combined_in_{sample}_with_ID.png", sample=SAMPLES),
		#PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png",
		#deeper
		#expand(PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES),
		#sniffles
		#expand(PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.vcf", sample=SAMPLES),
		#svim
		#expand(PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.vcf", sample=SAMPLES),
		#nanovar
		#expand(PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_INS_Variant_{sample}.vcf", sample=SAMPLES),
		#reads in insertion variants
		#expand(PROCESS+"VARIANTS/Reads_with_BLAST_SVIM_INS_Variant_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.fasta", sample=SAMPLES),
		#expand(PROCESS+"VARIANTS/BLASTN/Annotated_SNIFFLES_INS_Variant_{sample}.blastn", sample=SAMPLES),
		#new approach for insertion identification
		#expand(PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES)
		#expand(PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn", sample=SAMPLES)
		#expand(PROCESS+"FASTA/Insertion_{sample}_Vector.fa", sample=SAMPLES),
		#expand(PROCESS+"LOCALIZATION/Insertion_fasta_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/Insertion_MeanMods_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed", sample=SAMPLES),
		#expand(PROCESS+"METHYLATION/MeanModificiation_Insertion_{sample}.png", sample=SAMPLES),
		#BAM evaulation
		#expand(PROCESS+"MAPPING/Matches_different_in_precut_{sample}.bed", sample=SAMPLES),
		#expand(PROCESS+"EVAL/Matches_different_in_precut_FASTA_{sample}.fa", sample=SAMPLES),
		#expand(PROCESS+"EVAL/Matches_different_in_precut_{sample}.bed", sample=SAMPLES)

#actual filenames
def get_input_names(wildcards):
    return config["samples"][wildcards.sample]

#BAM Operations:
#Add rule to remove supplementary and secondary alignments
#samtools view -F 2304 -bo filtered.bam original.bam (via picard: supplementary alignment, not primary alignment): Added to prepare BAM rule

######
######
###### Insertion BAM: Convert to FASTA 
######
######
#bam with insertions
rule prepare_BAM: #sorts, index, and removes supllementary and secondary alignments
	input:
		get_input_names
	output:
		PROCESS+"MAPPING/{sample}_sorted.bam"
	shell:
		"""
		samtools sort {input} > {output} #| samtools view -F 2304 -o {output} #removal for a test run
		samtools index {output}
		"""
rule minimap_index:
	input:
		ref=config["ref_genome"] #cut _index
	output:
		index=PROCESS+"MAPPING/ref_genome_index.mmi"
	shell:
		"minimap2 -d {output.index} {input.ref}"

rule conserve_tags: #works
	input:
		bam=PROCESS+"MAPPING/{sample}_sorted.bam",
		minimapref=PROCESS+"MAPPING/ref_genome_index.mmi",
		ref=config["ref_genome"]
	output:
		PROCESS+"MAPPING/Conservedtags_{sample}_sorted.bam"
	shell:
		"""
		samtools bam2fq -T 'MM' {input.bam}| minimap2 -y -ax map-ont {input.minimapref} - | samtools sort > {output}
		"""

rule make_FASTA_with_tags:
	input:
		#fq=get_input_names
		fq=PROCESS+"MAPPING/Conservedtags_{sample}_sorted.bam"
	output:
		fasta=PROCESS+"FASTA/Conservedtags_WithTags_Full_{sample}.fa"
	run: 
		#shell("samtools bam2fq -T 'MM' {input} > {output}") #seqkit had to be installed with conda
		shell("samtools fasta {input} -o {output} > {output}") #takes in bam and outputs fasta, output has to be mentioned twice: 1 for the location, 2 for the option to write "both" reads to one fasta

rule make_FASTA_without_tags:
	input:
		fq=get_input_names
		#fq=PROCESS+"MAPPING/{sample}_sorted.bam"
	output:
		fasta=PROCESS+"FASTA/Full_{sample}.fa"
	run: 
		#shell("samtools bam2fq -T 'MM' {input} > {output}") #seqkit had to be installed with conda
		shell("samtools fasta {input} -o {output} > {output}") #takes in bam and outputs fasta, output has to be mentioned twice: 1 for the location, 2 for the option to write "both" reads to one fasta
######
######
###### "Clean" BAM: Cut-out FASTA to BAM via Mapping to reference 
######
######
		
#rule to create the BAM files with the non-insertion reads and the splitted read fragments
rule Non_insertion_mapping:
	input:
		fasta=PROCESS+"FASTA/Cleaved_{sample}_noVector.fa",
		genome=PROCESS+"MAPPING/ref_genome_index.mmi"
	output:
		PROCESS+"MAPPING/CutOut_{sample}_sorted.bam"
	shell:
		"""
		minimap2 -y -ax map-ont {input.genome} {input.fasta} | samtools sort |  samtools view -F 2304 -o {output} #added the removal of sec and suppl alignments 
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
		precut=PROCESS+"MAPPING/{sample}_sorted.bam",
		postcut=PROCESS+"MAPPING/CutOut_{sample}_sorted.bam" #Reads that contained an insertion before, are now marked with "_Insertion"
	output:
		postcut=PROCESS+"MAPPING/BasicMapping_{sample}.bed",
		precut=PROCESS+"MAPPING/Precut_{sample}.bed"
	run:
		shell("bedtools bamtobed -i {input.precut} > {output.precut}")
		shell("bedtools bamtobed -i {input.postcut} > {output.postcut}")  

rule reads_with_BLASTn_matches:
	input:
		refbed=PROCESS+"MAPPING/BasicMapping_{sample}.bed",
		#matchreads=PROCESS+"BLASTN/ID/ID_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed" 
	shell:
		"grep '_Insertion' {input.refbed} > {output}"
######
######
###### FASTA preparation: Cut out of blast-detected vector fragments and create new "cut-out" FASTA 
######
######

rule get_cleavage_sites_for_fasta:
	input:
		#PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
		PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		overlap=2*FRAG #this is the distance of the start-stop that is allowed to exist to still be combined; This should not be lower than FRAG!
	output:
		PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:
		vhf.splitting_borders(input[0],params[0], output[0])

rule split_fasta:
	input:
		breakpoints=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		#fasta=PROCESS+"FASTA/Full_{sample}.fa" #removed for test
		fasta=PROCESS+"FASTA/Conservedtags_WithTags_Full_{sample}.fa"
	params:
		mode="Join" #if each split FASTA substring should be used individually, use "Separated"
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
		PROCESS+"FASTA/Fragments/Forward_Backward_Vector.fa"
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
rule mapping_qc:
	input:
		#get_input_names
		PROCESS+"MAPPING/{sample}_sorted.bam" #before cut out, but after removal of secondary and supplementary alignments
	output:
		PROCESS+"MAPPING/BasicMapping_{sample}.qc"
	run:
		shell("samtools flagstats {input} > {output}")  
		
rule nBases_for_insertion_count:
	input:
		#get_input_names
		PROCESS+"MAPPING/{sample}_sorted.bam" #has to be without the cutout, since otherwise we don't get the total aligned bases
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

######
######
###### Base modification file preparation
######
######

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
		PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
	params: 
		10000 #this value has to be same as in mean methylation!
	output:
		temp(PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed")
	run:
		vhf.insertion_proximity(input[0], params[0], output[0])


#creates fasta only with insertions (+proximity of 10k up and downstream
rule insertion_methylation_proximity:
    input:
        prox = PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed",
        fasta = config["ref_genome"]
    output:
        temp(PROCESS+"METHYLATION/Insertion_methylation_proximity_{sample}.fa")
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.prox} > {output}"

#adds fasta to genomic gaps
rule fasta_to_insertion_proximity:
    input:
        prox = PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed",
        fasta = PROCESS+"METHYLATION/Insertion_methylation_proximity_{sample}.fa"
    output:
        temp(PROCESS+"METHYLATION/Insertion_fasta_proximity_{sample}.bed")
    run:
        vhf.add_sequence_column(input.prox, input.fasta, output[0])

#this rule is added for performance, so that the following pythons cript only reads in relevant parts of the methylation file and not the full genome

rule methyl_specific:
    input:
        insertions = PROCESS+"LOCALIZATION/Proximity_GenomicLocation_"+str(FRAG)+"_{sample}.bed", #bed file is based on different genomic coordinates :(
        methyl = PROCESS+"METHYLATION/Methyl_{sample}.bed"
    output:
        temp(PROCESS+"METHYLATION/Specific_Methyl_{sample}.bed")
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

######## INSERTION C_MODIFICATION

#creates fasta only with insertions
#fasta to insertion
rule fasta_to_insertion:
    input:
        ins = PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed", #exact coordinates
        fasta = PROCESS+"FASTA/Insertion_{sample}_Vector.fa"
    output:
        PROCESS+"LOCALIZATION/Insertion_fasta_{sample}.bed"
    run:
        vhf.add_insertion_sequence(input.ins, input.fasta, output[0])


#this rule is added for performance, so that the following pythons cript only reads in relevant parts of the methylation file and not the full genome

rule insertion_methyl_specific: #can be merged into rule above later
    input:
        insertions = PROCESS+"LOCALIZATION/Insertion_fasta_{sample}.bed",
        methyl = PROCESS+"METHYLATION/Methyl_{sample}.bed"
    output:
        PROCESS+"METHYLATION/Insertion_Specific_Methyl_{sample}.bed"
    shell:
        """
        bedtools intersect -wb -a {input.methyl} -b {input.insertions} > {output}
        """

#mean methylation in the insertion
rule insertion_mean_methylation:
	input:
		methbed=PROCESS+"METHYLATION/Insertion_Specific_Methyl_{sample}.bed",
		insertionfastabed=PROCESS+"LOCALIZATION/Insertion_fasta_{sample}.bed"
	params:
		window_size=50,
		max_distance=0
	output:
		PROCESS+"METHYLATION/Insertion_MeanMods_{sample}.bed"
	run: 
		vhf.methylation_in_insertion_proximity(input.methbed, input.insertionfastabed, params.window_size, params.max_distance,  output[0])

rule plot_insertion_mean_methylation:
	input:
		PROCESS+"METHYLATION/Insertion_MeanMods_{sample}.bed"
	params:
		window_size=50
	output:
		PROCESS+"METHYLATION/MeanModificiation_Insertion_{sample}.png"
	run:
		vhf.plot_modification_per_vectorlength(input[0], params.window_size, output[0])

######
######
###### Differential methylation
######
######

#GOAL: Use bed insertions from all exact matches and check for the methylation pattern in healthy and non-healthy
#rule 1: combine the files with ID column
#rule 2: Mean methylation CD123 vs UTD, once in UTD, once in ADS, once in CD123
#rule 3: Heatmap

rule combine_exact_matches:
	input:
		expand(PROCESS+"METHYLATION/Insertion_fasta_proximity_{sample}.bed", sample=SAMPLES)
	output:
		temp(PROCESS+"METHYLATION/All_Insertions/Insertion_fasta_proximity_combined.bed")
	run:
		vhf.combine_beds_add_ID(input, output[0]) 

rule combined_mean_methylation:
	input:
		methbed=PROCESS+"METHYLATION/Specific_Methyl_{sample}.bed",
		insertionfastabed=PROCESS+"METHYLATION/All_Insertions/Insertion_fasta_proximity_combined.bed"
	params:
		window_size=500,
		max_distance=10000
	output:
		temp(PROCESS+"METHYLATION/All_Insertions/MeanMods_Proximity_combined_in_{sample}.bed")
	run: 
		vhf.methylation_in_insertion_proximity(input.methbed, input.insertionfastabed, params.window_size, params.max_distance,  output[0])

rule add_ID:
	input:
		bed1=PROCESS+"METHYLATION/All_Insertions/MeanMods_Proximity_combined_in_{sample}.bed",
		bed2=PROCESS+"METHYLATION/All_Insertions/Insertion_fasta_proximity_combined.bed"
	output:
		PROCESS+"METHYLATION/All_Insertions/MeanMods_Proximity_combined_in_{sample}_with_Insertion_ID.bed"
	run:
		vhf.add_annotation_column_bed(input.bed1, input.bed2, output[0])
		
rule healty_combined_insertion_modification_heatmap:
	input:
		PROCESS+"METHYLATION/All_Insertions/MeanMods_Proximity_combined_in_{sample}_with_Insertion_ID.bed"
	params:
		window_size=500,
		max_distance=10000
	output:
		PROCESS+"METHYLATION/All_Insertions/Heatmap_MeanMods_combined_in_{sample}_with_ID.png"
	run:
		vhf.plot_modification_proximity(input[0], params.window_size, params.max_distance, output[0])

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
		temp(PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	run:
		shell("blastn -query {input.fasta} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs'") 

rule hardcode_blast_header:		
	input: 
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:	
		shell("echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov' | cat - {input} > {output}")
		
#BLASTN vector against human genome: Which vector parts are close to human sequences so that they might raise a false positivite BLAST match
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

rule blast_to_gff:
	input: 
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"BLASTN/"+str(FRAG)+"_VectorMatches_{sample}.gff"	
	run:
		vhf.blast2gff(input[0], output[0])	
######
######
###### Visualisations of intermediate results
######
######
	
rule chromosome_read_plots:
	input:
		bam=PROCESS+"MAPPING/CutOut_{sample}_sorted.bam",
		bed=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		#bed=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed"
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
		#PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		#PROCESS+"BLASTN/HUMANREF/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
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
		expand(PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed", sample=SAMPLES)
	output:
		PROCESS+"LOCALIZATION/Heatmap_Insertion_Chr.png"
	run:
		print(input)
		vhf.plot_bed_files_as_heatmap(input, output[0])

rule insertion_modification_heatmap:
	input:
		single=PROCESS+"METHYLATION/MeanMods_Proximity_{sample}.bed",
	params:
		window_size=500,
		max_distance=10000	
	output:
		singleout=PROCESS+"METHYLATION/Heatmap_MeanMods_Proximity_{sample}.png",
	run:
		vhf.plot_modification_proximity(input.single, params.window_size, params.max_distance, output.singleout)

######
######
###### Variant callers and overlap check of Insertions with BLAST matches
######
######

rule variant_sniffles:
	input:
		bam=PROCESS+"MAPPING/{sample}_sorted.bam",
		genome=config["ref_genome"]
	output:
		PROCESS+"VARIANTS/SNIFFLES/Variant_{sample}.vcf"
	shell:
		"sniffles -i {input.bam} --reference {input.genome} --output-rnames -v {output}"

rule svim_variants:
	input:
		bam=PROCESS+"MAPPING/{sample}_sorted.bam",
		genome=config["ref_genome"]
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
		bam=PROCESS+"MAPPING/{sample}_sorted.bam"
	output:
		PROCESS+"VARIANTS/NanoVar_{sample}/{sample}_sorted.nanovar.pass.vcf"
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
		nanovar=PROCESS+"VARIANTS/NanoVar_{sample}/{sample}_sorted.nanovar.pass.vcf",
		svim=PROCESS+"VARIANTS/SVIM_{sample}/variants.vcf", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/Variant_{sample}.vcf"
	output:
		nanovar=PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_INS_Variant_{sample}.bed",
		svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.bed", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.bed"
	run:
		#shell("convert2bed --input=vcf  < {input.nanovar} | grep 'INS' > {output.nanovar}") #--insertions does not wqork for some reason
		shell("vcf2bed --do-not-sort < {input.nanovar} | grep 'INS' > {output.nanovar}") #do not sort is important to generate the file but results in error downstream :/
		shell("convert2bed --input=vcf --insertions < {input.sniffles} > {output.sniffles}")
		shell("convert2bed --input=vcf --insertions < {input.svim} > {output.svim}")
		

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
		#shell("if ! grep -F -f {input.matches} {input.nanovar}; then echo not found; fi > {output.nanovar}")
		#shell("if ! grep -F -f {input.matches} {input.svim}; then echo not found; fi > {output.svim}")
		#shell("if ! grep -F -f {input.matches} {input.sniffles}; then echo not found; fi > {output.sniffles}")
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
''' temp block until necessary
rule blastn_bed_merger:
	input:
		#blast=PROCESS+"BLASTN/Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		blast=PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		bed=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed"
	output:
		merged=PROCESS+"LOCALIZATION/Exact_GenomicLocation_"+str(FRAG)+"_{sample}.bed"
	run:
		vhf.blastn_bed_merger(input.blast, input.bed, output.merged)
'''
rule exact_insertion_coordinates:
	input:
		#bed=PROCESS+"LOCALIZATION/GenomicLocation_"+str(FRAG)+"_{sample}.bed",
		bed=PROCESS+"MAPPING/BasicMapping_{sample}.bed", #full bed, maybe a inbetween step can be replaced!
		borders=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		diff = 50 #mean difference in the insertion start positions. If less than this threshhold, the first insertion will be treated as a representative
	output:
		out=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		out2=PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed" #only to get FASTA sequence
	run:
		vhf.exact_insertion_coordinates(input.borders, input.bed, params.diff, output.out, output.out2)

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
		shell("awk -F'\t' '$7>={params.threshold}' {input.blast} > {output.blast}")
		shell("awk -F'\t' '$11>={params.threshold}' {input.humanref} > {output.humanref}") 

#####CHecking for overlaps of Insertions detected by BLAST with Insertions detected by the variant callers:
#Strategy: Get READ IDs from blast matches, trace back to dorado basecalled bam reads, extract coordinates of the reads in the dorado bam; then
# check if there are overlaps. General problem: We must assume that the dorado basecalled bam has a low resolution when it comes to the INsertions
# -> this may alter the variant callers ability for the detection of insertions 


#Control of the workflow: How do the genomic coordinates of my two BAMs differ for the insertion vectors! # CUrrently it looks like all the insertion reads do not really differ pre and post cut -> either the cut out is not sensitive ennough (try with split reads) or the vector doesn't really alter the alignment!
rule check_mapping_pre_and_postcut:
	input:
		postcutbed=PROCESS+"MAPPING/BasicMapping_{sample}.bed",
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
		notinpostcut=PROCESS+"MAPPING/BasicMapping_{sample}.bed",
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
