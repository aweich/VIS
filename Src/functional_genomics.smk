#Functional genomics rules for VIS
#ORF rules
rule orf_prediction: #to do: add to bashrc
	input:
		PROCESS+"FASTA/Reads_with_Insertion_{sample}_Vector.fasta"
	output:
		proteinorf=PROCESS+"FUNCTIONALGENOMICS/ORF/Protein_Insertion_{sample}_Vector.orf",
		fastaorf=PROCESS+"FUNCTIONALGENOMICS/ORF/FASTA_Insertion_{sample}_Vector.orf"
	shell: #-n = ignore nested ORFs: Should give us only the bigger ones
		"""
		ORFfinder -s 0 -in {input} -outfmt 0 -out {output.proteinorf}
		ORFfinder -s 0 -n true -in {input} -outfmt 1 -out {output.fastaorf}
		"""

rule orf_reshape: #start and stop for - strand is interchanged to make them plottable
	input:
		orfs=PROCESS+"FUNCTIONALGENOMICS/ORF/FASTA_Insertion_{sample}_Vector.orf"
	output:
		bed=PROCESS+"FUNCTIONALGENOMICS/ORF/BED_{sample}_Vector.orf"
	run:
		vhf.orf_reshape(input.orfs, output.bed)

###### ORF to protein BLAST
rule protein_blast:
	input:
		PROCESS+"FUNCTIONALGENOMICS/ORF/Protein_Insertion_{sample}_Vector.orf"
	params: 
		db=config["proteindb"]
	output:
		temp(PROCESS+"FUNCTIONALGENOMICS/ORF/PROTEINBLAST/Protein_Insertion_NoHeaderBlast_{sample}.temp")
	shell:
		"blastp -query {input} -db {params.db} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qlen slen qstart qend length mismatch pident qcovs evalue bitscore' -taxids 9606"
		
				
rule protein_blast_header:		
	input: 
		PROCESS+"FUNCTIONALGENOMICS/ORF/PROTEINBLAST/Protein_Insertion_NoHeaderBlast_{sample}.temp"
	output:
		PROCESS+"FUNCTIONALGENOMICS/ORF/PROTEINBLAST/ORFs_{sample}.proteinblast"
	shell:	
		"echo -e 'QueryID\tSubjectID\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tEvalue\tBitscore' | cat - {input} > {output}"

#Localization
#this ensures that a writable R lib exists where the needed packages will be installed to!
shell.prefix('export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.0 && ')
        
rule chromosome_read_plots: 
	input:
		bed=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		#bam=PROCESS+"MAPPING/NoVectorAlignments_Postcut_{sample}_sorted.bam", #if normal bam would work here, we could cut several upstream rules!
		bam=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam",
		H3K4Me1=config["ucsc_H3K4Me1"],
		H3K4Me3=config["ucsc_H3K4Me3"],
		H3K27Ac=config["ucsc_H3K27Ac"],
		gtf=config["ucsc_Genes_gtf"],
		TF=config["ucsc_TF"],	
	output:
		outpath=directory(PROCESS+"FUNCTIONALGENOMICS/LOCALIZATION/" + str(FRAG)+"_{sample}")
	params:
		buffer=50000
	shell: 
		r"""
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -iH3K4Me1 {input.H3K4Me1} -iH3K4Me3 {input.H3K4Me3} -iH3K27Ac {input.H3K27Ac} -igtf {input.gtf} -iTF {input.TF} -buffer {params.buffer} -o {output.outpath}   
		"""	
#Genes and TF new
rule distance_to_elements:
	input:
		insertions=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		genes=config["ucsc_Genes"],
		tf=config["ucsc_TF"],
		cd4=config["sedb_cd4"],
		cd8=config["sedb_cd8"]
	params:
		distances=[0,1000,5000,10000,50000]
	output:
		genes=PROCESS+"FUNCTIONALGENOMICS/Distance_to_Genes_" + str(FRAG)+"_{sample}.bed",
		#tf=PROCESS+"FUNCTIONALGENOMICS/Distance_to_TF_" + str(FRAG)+"_{sample}.bed",
		#cd4=PROCESS+"FUNCTIONALGENOMICS/Distance_to_CD4_SE_" + str(FRAG)+"_{sample}.bed",
		#cd8=PROCESS+"FUNCTIONALGENOMICS/Distance_to_CD8_SE_" + str(FRAG)+"_{sample}.bed"
	run:
		vhf.calc_distance_to_element(input.insertions, input.genes, params.distances, output.genes)
		#vhf.calc_distance_to_element(input.insertions, input.tf, params.distances, output.tf)
		#vhf.calc_distance_to_element(input.insertions, input.cd4, params.distances, output.cd4)
		#vhf.calc_distance_to_element(input.insertions, input.cd8, params.distances, output.cd8)

rule plot_distance_to_elements:
	input:
		genes=PROCESS+"FUNCTIONALGENOMICS/Distance_to_Genes_" + str(FRAG)+"_{sample}.bed",
		#tf=PROCESS+"FUNCTIONALGENOMICS/Distance_to_TF_" + str(FRAG)+"_{sample}.bed",
		#cd4=PROCESS+"FUNCTIONALGENOMICS/Distance_to_CD4_SE_" + str(FRAG)+"_{sample}.bed",
		#cd8=PROCESS+"FUNCTIONALGENOMICS/Distance_to_CD8_SE_" + str(FRAG)+"_{sample}.bed"
	params:
		distances=[0,1000,5000,10000,50000]
	output:
		genes=PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png",
		#tf=PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_TF_" + str(FRAG)+"_{sample}.png",
		#cd4=PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_CD4_SE_" + str(FRAG)+"_{sample}.png",
		#cd8=PROCESS+"FUNCTIONALGENOMICS/Plot_Distance_to_CD8_SE_" + str(FRAG)+"_{sample}.png"
	run:
		vhf.plot_element_distance(input.genes, params.distances, output.genes)
		#vhf.plot_element_distance(input.tf, params.distances, output.tf)
		#vhf.plot_element_distance(input.cd4, params.distances, output.cd4)
		#vhf.plot_element_distance(input.cd8, params.distances, output.cd8)
	
'''
will be removed if not needed again by 14.08.2024
#Genes and TF
rule proximity_generator:
	input:
		bed=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
	params:
		offsets=[0,1000,5000,10000,500000]
	output:
		out=temp(PROCESS+"LOCALIZATION/Proximity_to_ExactInsertions_{sample}.bed")
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
		"bedtools intersect -wa -a {input.bed} -wb -b {input.TF} > {output}"	
		#"bedtools intersect -wa -a {input.bed} -wb -b {input.TF} | cut -f1,2,3,4,5,9 | sort -k6 > {output}"	

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

### Distance of VIS to genes, TSS, miRNAs, and TF
rule sort_insertion_file:
	input:
		PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed"
	output:
		PROCESS+"LOCALIZATION/Sorted_ExactInsertions_{sample}.bed"
	run:
		shell("sort -k1,1 -k2,2n {input} > {output}")

rule distance_to_regulation:
	input:
		insertions=PROCESS+"LOCALIZATION/Sorted_ExactInsertions_{sample}.bed",
		ref=config["ref_genome_ctrl"]
	params:
		genes=config["ucsc_Genes"],
		tf=config["ucsc_TF"],
		tss=config["ucsc_tss"], 
		mirna=config["ucsc_mirna"]
	output:
		PROCESS+"FUNCTIONALGENOMICS/Functional_distances_to_Insertions_{sample}.bed"
	run:
		shell("bedtools closest -a {input.insertions} -b {params.genes} {params.tf} {params.tss} {params.mirna} -names genes tf tss mirna -D a > {output}")

rule plot_all_regulation:
	input:
		PROCESS+"FUNCTIONALGENOMICS/Functional_distances_to_Insertions_{sample}.bed"
	output:
		PROCESS+"FUNCTIONALGENOMICS/Functional_distances_to_Insertions_{sample}.png"
	run:
		vhf.plot_all_elements_by_distance(input[0], output[0])
