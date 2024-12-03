######
######
###### Quality control and normalization
######
######
 
#FRAG=config["fragment_size"]

#full stats of input data		
rule nanoplot:
	input:
		PROCESS+"MAPPING/Precut_{sample}_sorted.bam"
	output:
		PROCESS+"QC/Nanoplot/{sample}/NanoStats.txt"
	params:
		outdir=PROCESS+"QC/Nanoplot/{sample}/"
	shell: 
		"""
		NanoPlot --bam {input} -o {params.outdir}
		touch {output}
		"""

'''
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
        tempdir=PROCESS+"temp_fastq_{sample}"  # Temporary directory for intermediate files
    output:
        fastq=PROCESS + "QC/FASTQC/{sample}_filtered.fastq"
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
        PROCESS + "QC/FASTQC/{sample}_filtered.fastq"
    params:
        prefix=PROCESS + "QC/FASTQC/readlevel_{sample}/{sample}_read_"
    output:
        directory(PROCESS + "QC/FASTQC/readlevel_{sample}/")
    shell:
        """
        mkdir -p {output}

        # Split the FASTQ file into individual entries and replace spaces with underscores in the read names
        awk 'BEGIN {{ OFS=""; }}
        NR % 4 == 1 {{
            if (filename) close(filename);
            # Replace spaces in the header (after @) with underscores
            $0 = "@" substr($0, 2);
            gsub(" ", "_", $0);  # Replace spaces with underscores
            filename = "{params.prefix}" substr($0, 2) ".fastq";  # Use modified header as filename
        }}
        {{ print > filename }}' {input}

        # Run FastQC on each split FASTQ file
        for seq in {params.prefix}*; do
            fastqc -f fastq --noextract -o {output} "$seq"
        done
        """

rule multiqc:
    input: 
        fastqc=expand(PROCESS+"QC/FASTQC/readlevel_{sample}/", sample=SAMPLES),
        nanoplot=expand(PROCESS+"QC/Nanoplot/{sample}/NanoStats.txt", sample=SAMPLES)
    params:
        qc_output_dir=PROCESS+"QC/",
        qc_report_location=FINAL+"QC/"
    output: 
        PROCESS + "QC/multiqc_report.html",
        FINAL + "QC/multiqc_report.html"
    shell:
        """
        multiqc {input.fastqc} --dirs {input.nanoplot} --force -o {params.qc_output_dir}
        
        # copy report to final location
        cp {output[0]} {output[1]}
        """

### Read level overview of mapping quality before and after the cut out of the insertions
rule extract_mapping_quality:
    input:
        bam=PROCESS + "MAPPING/Precut_{sample}_sorted.bam",
        bam2=PROCESS + "MAPPING/Postcut_{sample}_sorted.bam",
        bam3=PROCESS+"MAPPING/Postcut_{sample}_unfiltered_sorted.bam",
        readnames=PROCESS + "BLASTN/Readnames_" + str(FRAG) + "_VectorMatches_{sample}.txt"
    params:
        tempdir=PROCESS+"temp_mapping_{sample}"
    output:
        quality_scores=temp(PROCESS + "QC/{sample}_precut_mapping_quality.txt"),
        quality_scores2=temp(PROCESS + "QC/{sample}_postcut_mapping_quality.txt"),
        quality_scores3=temp(PROCESS + "QC/{sample}_postcut_unfiltered_mapping_quality.txt")
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
rule finalize_mapping_quality:
    input:
        quality_scores_pre=PROCESS + "QC/{sample}_precut_mapping_quality.txt",
        quality_scores_filtered=PROCESS + "QC/{sample}_postcut_unfiltered_mapping_quality.txt",
        quality_scores_post=PROCESS + "QC/{sample}_postcut_mapping_quality.txt"
    params:
        prefixes=["Precut", "Postcut", "Postcut_filtered"]
    output:
        outfile=PROCESS + "QC/MAPQ/Insertions_{sample}_mapq.txt"
    run:
        vhf.join_read_mapq(input[0:3], params.prefixes, output.outfile)

rule generate_mapq_heatmap:
    input:
        table=PROCESS+"QC/MAPQ/Insertions_{sample}_mapq.txt"
    output:
        heatmap=PROCESS+"QC/MAPQ/{sample}_mapq_heatmap_image.png"
    run:
        vhf.plot_mapq_changes(input.table, output.heatmap)

#### Visualize fragmentation and longest consecutive interval per read

rule fragmentation_distribution_plots:
	input:
		PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		PROCESS+"BLASTN/HUMANREF/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		FRAG
	output:
		outpath=directory(FINAL+"QC/Fragmentation/Insertions/insertions_" + str(FRAG)+"_{sample}"),
		outpath2=directory(FINAL+"QC/Fragmentation/Reference/reference_" + str(FRAG)+"_{sample}")
	run:
		shell("mkdir {output.outpath}")
		vhf.fragmentation_match_distribution(input[0], params[0], output[0])
		vhf.fragmentation_read_match_distribution(input[0], params[0], output[0])
		shell("mkdir {output.outpath2}")
		vhf.fragmentation_match_distribution(input[1], params[0], output[1])
		vhf.fragmentation_read_match_distribution(input[1], params[0], output[1])

rule detailed_insertion_length_plot:
    input:
        matches=PROCESS+"BLASTN/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
    params: 
        buffer=3*FRAG,
        threshold=config["MinInsertionLength"]
    output:
        outpath=directory(FINAL+"QC/Fragmentation/Longest_Interval/{sample}/")
    run:
        shell("mkdir -p {output.outpath}")
        
        vhf.find_and_plot_longest_blast_interval(
            input.matches,
            params.buffer,
            params.threshold,
            output.outpath
        )
