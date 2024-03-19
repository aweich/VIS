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
rule chromosome_read_plots: 
	input:
		bed=PROCESS+"LOCALIZATION/ExactInsertions_{sample}.bed",
		#bam=PROCESS+"MAPPING/NoVectorAlignments_Postcut_{sample}_sorted.bam", #if normal bam would work here, we could cut several upstream rules!
		bam=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam",
		H3K4Me1=config["ucsc_H3K4Me1"],
		H3K4Me3=config["ucsc_H3K4Me3"],
		H3K27Ac=config["ucsc_H3K27Ac"],
		gtf=config["ucsc_Genes_gtf"],
		TF=config["ucsc_TF"]	
	output:
		outpath=directory(PROCESS+"FUNCTIONALGENOMICS/LOCALIZATION/" + str(FRAG)+"_{sample}")
	params:
		buffer=50000
	shell: 
		r"""
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -iH3K4Me1 {input.H3K4Me1} -iH3K4Me3 {input.H3K4Me3} -iH3K27Ac {input.H3K27Ac} -igtf {input.gtf} -iTF {input.TF} -buffer {params.buffer} -o {output.outpath}   
		"""	
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

