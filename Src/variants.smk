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
		ref=config["ref_genome_ctrl"],
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
		svim=PROCESS+"VARIANTS/SVIM_{sample}/variants.vcf", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/Variant_{sample}.vcf"
	output:
		nanovar=PROCESS+"VARIANTS/NanoVar_{sample}/Nanovar_Variant_{sample}.bed",
		svim=PROCESS+"VARIANTS/SVIM_{sample}/SVIM_INS_Variant_{sample}.bed", 
		sniffles=PROCESS+"VARIANTS/SNIFFLES/SNIFFLES_INS_Variant_{sample}.bed"
	run:
		#shell("convert2bed --input=vcf  < {input.nanovar} > {output.nanovar}")
		shell("vcf2bed --do-not-sort < {input.nanovar} > {output.nanovar}") #do not sort is important to generate the file but results in error downstream :/
		shell("convert2bed --input=vcf --insertions < {input.sniffles} > {output.sniffles}")
		shell("convert2bed --input=vcf --insertions < {input.svim} > {output.svim}")
		
#####
# Outdated rules that check if insertions match with BLASTn of vector
#Will be removed if not used by april 5
#####
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
		

