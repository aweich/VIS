######
######
###### Quality control and normalization
######
######
#Nanoplot for general information about the sequencing yield of each sample
#MultiQC for quality of the sequencing yield
#MappingQuality?
#Mosdepth for genome-wide coverage
#CIGAR tools for the identification of high confidence insertions
 
FRAG=100
#actual filenames
def get_input_names(wildcards):
    return config["samples"][wildcards.sample]


rule make_FASTA_with_tags:
	input:
		fq=get_input_names
	output:
		fasta=PROCESS+"FASTA/Conserved_WithTags_Full_{sample}.fa"
	run: 
		shell("samtools bam2fq -T '*' {input} > {output}") 
		
rule nanoplot:
	input:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
	output:
		PROCESS+"QC/Nanoplot/{sample}/Non_weightedHistogramReadlength.png"
	params:
		outdir=PROCESS+"QC/Nanoplot/{sample}/"
	shell: 
		"""
		NanoPlot --bam {input} -o {params.outdir}
		touch {output}
		"""

rule quality_of_mapping_qc:
	input:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam" #before cut out, but after removal of secondary and supplementary alignments
	output:
		#short=PROCESS+"QC/readlevel_{sample}/short_{sample}.qc",
		long=PROCESS+"QC/mapping_qc_{sample}/{sample}.flagstats"
	shell:
		"""
		samtools stats {input} > {output.long}    
		"""
		#samtools flagstats {input} > {output.short}
rule normalisation_for_insertion_count:
	input:
		insertions=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt",
		fasta=PROCESS+"FASTA/Full_{sample}.fa" #for N50
	params:
		scale=3000000000 #IPM #IPG
	output:
		PROCESS+"QC/Normalisation_IPG_{sample}.txt"
	run:
		vhf.insertion_normalisation(input.insertions, params[0], input.fasta, output[0])
			
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

###cigar operations: Checks which of the reads have an insertion >1000 in the cigar upon primary mapping that matches with some part of the vector.
### These are high confidence matches!

rule mapping_for_cigar:
	input:
		bam=get_input_names,
		ref=config["ref_genome_ctrl"]
	output:
		PROCESS+"QC/CIGAR/HealthyRef_cigar_{sample}_sorted.bam"
	shell:
		"""
		samtools bam2fq -T '*' {input.bam}| minimap2 -y -ax map-ont {input.ref} - | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		"""

rule cigar_strings:
	input:
		pre=PROCESS+"MAPPING/Precut_{sample}_sorted.bam", #PROCESS+"QC/CIGAR/HealthyRef_cigar_{sample}_sorted.bam"
		post=PROCESS+"MAPPING/Postcut_{sample}_sorted.bam",
		matches=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt",
	output:
		pre=temp(PROCESS+"QC/CIGAR/cigar_{sample}_precut.txt"),
		names=temp(PROCESS+"QC/CIGAR/cigar_{sample}_precut_names.txt"),
		post=temp(PROCESS+"QC/CIGAR/cigar_{sample}_postcut.txt"),
	shell:
		"""
		samtools view {input.pre} | cut -f 1,3,4,6 | grep '[0-9]\{{4\}}I' > {output.pre}
		samtools view {input.pre} | cut -f 1,3,4,6 | grep '[0-9]\{{4\}}I' | cut -f 1 > {output.names}
		samtools view {input.post} | grep -f {input.matches} | cut -f 1,3,4,6 > {output.post}
		"""
		
rule cigar_fasta:
	input:
		fasta=PROCESS+"FASTA/Full_{sample}.fa",
		matches=PROCESS+"QC/CIGAR/cigar_{sample}_precut_names.txt",
	output:
		PROCESS+"QC/CIGAR/1000I_cigar_FASTA_{sample}.fa"
	run: 
		shell("seqkit grep -r -f {input.matches} {input.fasta} -o {output}")

rule cigar_fasta_blast:
	input:
		fasta=PROCESS+"QC/CIGAR/1000I_cigar_FASTA_{sample}.fa",
		dummy=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa.ndb" #provokes the building of the database first!
	params:
		vector=PROCESS+"FASTA/Fragments/" + str(FRAG) + "_Vector_fragments.fa"
	output:
		temp(PROCESS+"QC/CIGAR/BLASTN/cigar"+str(FRAG)+"_VectorMatches_{sample}.blastn")
	run:
		shell("blastn -query {input.fasta} -db {params.vector} -out {output} -evalue 1e-5 -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'") 

rule cigar_hardcode_blast_header:		
	input: 
		PROCESS+"QC/CIGAR/BLASTN/cigar"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		PROCESS+"QC/CIGAR/BLASTN/cigar_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	run:	
		shell("echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input} > {output}")	
		
rule cigar_reads_and_blast_insertion: #just checks if read id is in the cleavage sites file
	input:
		fasta1=PROCESS+"QC/CIGAR/1000I_cigar_FASTA_{sample}.fa",
		breakpoints=PROCESS+"BLASTN/CleavageSites_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	output:
		temp(PROCESS+"QC/CIGAR/BLASTN/Reads_cigar1000_VectorMatches_{sample}.fasta"),
	run:
		vhf.reads_with_insertions(input.fasta1, input.breakpoints, output[0])		

rule fasta_to_fastq: #first line extracts the fastq entries, second line separates them to make the input for fastaqc easier #this is basically the perfect setup as a multifastqc input
	input:
		fastq=PROCESS+"FASTA/Conserved_WithTags_Full_{sample}.fa",
		fasta=PROCESS+"QC/CIGAR/BLASTN/Reads_cigar1000_VectorMatches_{sample}.fasta"
	output:
		full=PROCESS+"QC/CIGAR/Reads_with_longInsertions_and_vector_{sample}.fastq"
	shell:
		'''
		# Check if the input fasta file is empty
        	if [ ! -s {input.fasta} ]; then
            		echo "No CIGAR - BLAST matches. Creating an empty output file."
            		touch {output.full}
            	else
			grep '>' {input.fasta} | cut -c 2- | seqkit grep -f - {input.fastq} -o {output.full}
			csplit -z {output.full} $(awk '/^@/ {{print NR }}' {output.full}) '/^@/' {{*}} --prefix={output.full}
		fi
		'''

rule bam_coverage: 
	input:
		bam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
	output:
		covbed=PROCESS+"QC/Coverage/Genomecoverage_{sample}.bed"
	run: 
		shell("bedtools genomecov -ibam {input} -bga > {output}")

### Read level fastqc analysis

rule extract_fastq_insertions:
    input:
    	bam=PROCESS+"MAPPING/Precut_{sample}_sorted.bam",
        readnames=PROCESS+"BLASTN/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt"
    params:
        tempdir=PROCESS+"temp_{sample}"  # Temporary directory for intermediate files
    output:
        fastq=PROCESS + "QC/{sample}_filtered.fastq"
    shell:
        '''
        mkdir -p {params.tempdir}

        # Extract header from BAM or generate it from reference genome
        samtools view -H {input.bam} > {params.tempdir}/header.sam

        samtools view -N {input.readnames} {input.bam} > {params.tempdir}/filtered_reads.sam
        
        # Concatenate all filtered reads and add header
        cat {params.tempdir}/header.sam {params.tempdir}/filtered_reads.sam | samtools view -b > {params.tempdir}/filtered_withheader.bam
        # Convert to FASTQ
        samtools fastq {params.tempdir}/filtered_withheader.bam > {output.fastq}

        # Clean up temporary files
        rm -r {params.tempdir}
        '''

rule read_level_fastqc:
    input:
        PROCESS+"QC/{sample}_filtered.fastq"
    params:
        prefix=PROCESS + "QC/readlevel_{sample}/{sample}_read_"
    output:
        directory(PROCESS + "QC/readlevel_{sample}/")
    shell:
        """
        mkdir -p {output}

        # Split the FASTQ file into individual entries
        csplit -z -f {params.prefix} {input} '/^@/' '{{*}}'

        # Run FastQC on each split FASTQ file
        for seq in {params.prefix}*; do
            fastqc -f fastq --noextract -o {output} "$seq"
        done
        """

rule multiqc:
    input: 
        fastqc=expand(PROCESS + "QC/readlevel_{sample}/", sample=SAMPLES)
    output: 
        PROCESS+"QC/Insertion_Reads_multiqc_report.html"
    wrapper: "0.31.1/bio/multiqc"

### Read level overview of mapping quality before and after the cut out of the insertions
rule extract_mapping_quality:
    input:
        bam=PROCESS + "MAPPING/Precut_{sample}_sorted.bam",
        bam2=PROCESS + "MAPPING/Postcut_{sample}_sorted.bam",
        bam3=PROCESS+"MAPPING/Postcut_{sample}_unfiltered_sorted.bam",
        readnames=PROCESS + "BLASTN/Readnames_" + str(FRAG) + "_VectorMatches_{sample}.txt"
    params:
        tempdir=PROCESS+"temp_{sample}"
    output:
        quality_scores=PROCESS + "QC/{sample}_precut_mapping_quality.txt",
        quality_scores2=PROCESS + "QC/{sample}_postcut_mapping_quality.txt",
        quality_scores3=PROCESS + "QC/{sample}_postcut_unfiltered_mapping_quality.txt"
    shell:
        '''
        mkdir -p {params.tempdir}

        # Extract reads of interest
        samtools view {input.bam} | grep -F -f {input.readnames} > {params.tempdir}/temp_precut_reads.sam
        samtools view {input.bam2} | grep -F -f {input.readnames} > {params.tempdir}/temp_postcut_reads.sam
        samtools view {input.bam3} | grep -F -f {input.readnames} > {params.tempdir}/temp_postcut_unfiltered_reads.sam

        # Extract mapping quality column (5th field) and read name
        awk '{{print $1,$3,$5}}' {params.tempdir}/temp_precut_reads.sam > {output.quality_scores}
        awk '{{print $1,$3,$5}}' {params.tempdir}/temp_postcut_reads.sam > {output.quality_scores2}
        awk '{{print $1,$3,$5}}' {params.tempdir}/temp_postcut_unfiltered_reads.sam > {output.quality_scores3}

        # Clean up
        rm -r {params.tempdir}
        '''
