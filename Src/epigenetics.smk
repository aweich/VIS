rule methylation_bedMethyl:
	input:
		bam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
		#bam=PROCESS+"MAPPING/NoVectorAlignments_Precut_{sample}_sorted.bam"
		#ref=config["ref_genome"]
	output:
		PROCESS+"METHYLATION/Methyl_{sample}.bed"
	shell:
		"modkit pileup {input.bam} --filter-threshold C:0.8 {output}" #high threshhold for C modifications, combines all C mods into one
		#"modkit pileup {input.bam} {output} --cpg --ref {input.ref}" 
'''
rule remove_unmapped_from_insertion_list:
	input:
		PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
	output:
		temp(PROCESS+"LOCALIZATION/Onlymapped_ExactInsertions_{sample}.bed")
	run:
		shell("cat {input} | grep -v 'CD' > {output}")
'''
rule insertion_proximity:
	input:
		PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
		#PROCESS+"LOCALIZATION/Onlymapped_ExactInsertions_{sample}.bed"
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
        fasta = config["ref_genome_ctrl"] #healhy genome here, since the one with the vector throws error
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
rule insertion_read_coordinates:
    input:
    	cleavage=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
    	postcut=PROCESS+"MAPPING/Precut_{sample}.bed"
    output:
    	PROCESS+"METHYLATION/Insertion_Read_Coordinates_{sample}.bed"
    run:
    	 vhf.full_coordinates_bed(input.cleavage, input.postcut, output[0])	
	
rule fasta_to_insertion: #this rule must return a bed file with the cooridnates of the mapped FASTA and the FASTA itself
    input:
        #ins = PROCESS+"LOCALIZATION/ExactInsertions_{sample}_full_coordinates_for_methylation.bed", #exact coordinates
        ins=PROCESS+"METHYLATION/Insertion_Read_Coordinates_{sample}.bed",
        fasta = PROCESS+"FASTA/Reads_with_Insertion_{sample}_Vector.fasta"
    output:
        PROCESS+"METHYLATION/Insertion_fasta_{sample}.bed"
    run:
        vhf.add_insertion_sequence(input.ins, input.fasta, output[0])


#this rule is added for performance, so that the following pythons cript only reads in relevant parts of the methylation file and not the full genome

rule insertion_methyl_specific: #can be merged into rule above later
    input:
        insertions = PROCESS+"METHYLATION/Insertion_fasta_{sample}.bed",
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
		insertionfastabed=PROCESS+"METHYLATION/Insertion_fasta_{sample}.bed"
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
		#outpath=PROCESS+"METHYLATION/MeanModificiation_Insertion_{sample}.png"
		outpath = directory(PROCESS+"METHYLATION/PLOTS/InsertionRead_{sample}/")
	run:
		shell("mkdir {output.outpath}")
		vhf.plot_modification_per_vectorlength(input[0], params.window_size, output.outpath)

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
		
rule healthy_combined_insertion_modification_heatmap:
	input:
		PROCESS+"METHYLATION/All_Insertions/MeanMods_Proximity_combined_in_{sample}_with_Insertion_ID.bed"
	params:
		window_size=500,
		max_distance=10000
	output:
		PROCESS+"METHYLATION/All_Insertions/Heatmap_MeanMods_combined_in_{sample}_with_ID.png"
	run:
		vhf.plot_modification_proximity(input[0], params.window_size, params.max_distance, output[0])

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

