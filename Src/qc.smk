######
######
###### Quality control and normalization
######
######
 
#FRAG=config["fragment_size"]

#full stats of input data		
rule nanoplot:
	input:
		PROCESS+"mapping/Precut_{sample}_sorted.bam"
	output:
		PROCESS+"qc/nanoplot/{sample}/NanoStats.txt"
	log:
		log=PROCESS+"log/qc/nanoplot/{sample}.log"
	params:
		outdir=directory(PROCESS+"qc/nanoplot/{sample}/"),
	conda:
		"../envs/VIS_nanoplot_env.yml"
	shell: 
		"""
		(
		NanoPlot --bam {input} -o {params.outdir}
		touch {output}
		) > {log.log} 2>&1
		"""

rule bam_coverage: 
	input:
		bam=PROCESS+"mapping/Precut_{sample}_sorted.bam"
	log:
		log=PROCESS+"log/qc/bam_coverage/{sample}.log"
	output:
		covbed=PROCESS+"qc/Coverage/Genomecoverage_{sample}.bed"
	conda:
		"../envs/VIS_bedtools_env.yml"
	shell: 
		"""
		(
		bedtools genomecov -ibam {input} -bga > {output}
		) > {log.log} 2>&1
		"""

### Read level fastqc analysis

rule extract_fastq_insertions:
    input:
    	bam=PROCESS+"mapping/Precut_{sample}_sorted.bam",
        readnames=PROCESS+"blastn/Readnames_"+str(FRAG)+"_VectorMatches_{sample}.txt"
    params:
        tempdir=PROCESS+"temp_fastq_{sample}"  # Temporary directory for intermediate files
    log:
    	log=PROCESS+"log/qc/extract_fastq_insertions/{sample}.log"
    output:
        fastq=PROCESS + "qc/fastqc/{sample}_filtered.fastq"
    conda:
    	"../envs/VIS_samtools_env.yml"
    shell:
        '''
        (
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
        ) > {log.log} 2>&1
        '''

rule read_level_fastqc:
    input:
        PROCESS + "qc/fastqc/{sample}_filtered.fastq"
    params:
        prefix=PROCESS + "qc/fastqc/readlevel_{sample}/{sample}_read_"
    log:
    	log=PROCESS+"log/qc/read_level_fastqc/{sample}.log"
    output:
        directory(PROCESS + "qc/fastqc/readlevel_{sample}/")
    conda:
    	f"{cwd}/envs/VIS_fastqc_env.yml"
    shell:
        """
        (
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

        # Run Fastqc on each split FASTQ file
        for seq in {params.prefix}*; do
            fastqc -f fastq --noextract -o {output} "$seq"
        done
        ) > {log.log} 2>&1
        """

rule multiqc:
    input: 
        fastqc=expand(PROCESS+"qc/fastqc/readlevel_{sample}/", sample=SAMPLES),
        nanoplot=expand(PROCESS+"qc/nanoplot/{sample}/NanoStats.txt", sample=SAMPLES)
    params:
        qc_output_dir=PROCESS+"qc/",
        qc_report_location=FINAL+"qc/"
    log:
    	log=PROCESS+"log/qc/multiqc/out.log"
    output: 
        PROCESS + "qc/multiqc_report.html",
        report(FINAL + "qc/multiqc_report.html")
    conda:
    	f"{cwd}/envs/VIS_multiqc_env.yml"
    shell:
        """
        (
        multiqc {input.fastqc} --dirs {input.nanoplot} --force -o {params.qc_output_dir}
        
        # copy report to final location
        cp {output[0]} {output[1]}
        ) > {log.log} 2>&1
        """

### Read level overview of mapping quality before and after the cut out of the insertions
rule extract_mapping_quality:
    input:
        bam=PROCESS + "mapping/Precut_{sample}_sorted.bam",
        bam2=PROCESS + "mapping/Postcut_{sample}_sorted.bam",
        bam3=PROCESS+"mapping/Postcut_{sample}_unfiltered_sorted.bam",
        readnames=PROCESS + "blastn/Readnames_" + str(FRAG) + "_VectorMatches_{sample}.txt"
    params:
        tempdir=PROCESS+"temp_mapping_{sample}"
    log:
    	log=PROCESS+"log/qc/extract_mapping_quality/{sample}.log"
    output:
        quality_scores=temp(PROCESS + "qc/{sample}_precut_mapping_quality.txt"),
        quality_scores2=temp(PROCESS + "qc/{sample}_postcut_mapping_quality.txt"),
        quality_scores3=temp(PROCESS + "qc/{sample}_postcut_unfiltered_mapping_quality.txt")
    conda:
    	"../envs/VIS_samtools_env.yml"
    shell:
        '''
        (
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
        ) > {log.log} 2>&1
        '''
rule finalize_mapping_quality:
    input:
        quality_scores_pre=PROCESS + "qc/{sample}_precut_mapping_quality.txt",
        quality_scores_filtered=PROCESS + "qc/{sample}_postcut_unfiltered_mapping_quality.txt",
        quality_scores_post=PROCESS + "qc/{sample}_postcut_mapping_quality.txt"
    params:
        prefixes=["Precut", "Postcut", "Postcut_filtered"]
    log:
    	log=PROCESS+"log/qc/finalize_mapping_quality/{sample}.log"
    output:
        outfile=PROCESS + "qc/mapq/Insertions_{sample}_mapq.txt"
    run:
        vhf.join_read_mapq(input[0:3], params.prefixes, output.outfile, log.log)

rule generate_mapq_heatmap:
    input:
        table=PROCESS+"qc/mapq/Insertions_{sample}_mapq.txt"
    log:
    	log=PROCESS+"log/qc/generate_mapq_heatmap/{sample}.log"
    output:
        heatmap=report(PROCESS+"qc/mapq/{sample}_mapq_heatmap_image.png")
    run:
        vhf.plot_mapq_changes(input.table, output.heatmap, log.log)

#### Visualize fragmentation and longest consecutive interval per read

rule fragmentation_distribution_plots:
	input:
		PROCESS+"blastn/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn",
		PROCESS+"blastn/humanref/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
	params:
		FRAG
	log:
		log1=PROCESS+"log/qc/fragmentation_distribution_plots/fragmentation_match_distribution_{sample}.log",
		log2=PROCESS+"log/qc/fragmentation_distribution_plots/fragmentation_read_match_distribution_{sample}.log",
		log3=PROCESS+"log/qc/fragmentation_distribution_plots/fragmentation_match_distribution_{sample}.log",
		log4=PROCESS+"log/qc/fragmentation_distribution_plots/fragmentation_read_match_distribution_{sample}.log"
	output:
		outpath=directory(FINAL+"qc/Fragmentation/Insertions/insertions_" + str(FRAG)+"_{sample}"),
		outpath2=directory(FINAL+"qc/Fragmentation/Reference/reference_" + str(FRAG)+"_{sample}")
	run:
		shell("mkdir {output.outpath}")
		vhf.fragmentation_match_distribution(input[0], params[0], output[0], log.log1)
		vhf.fragmentation_read_match_distribution(input[0], params[0], output[0], log.log2)
		shell("mkdir {output.outpath2}")
		vhf.fragmentation_match_distribution(input[1], params[0], output[1], log.log3)
		vhf.fragmentation_read_match_distribution(input[1], params[0], output[1], log.log4)

rule detailed_fragmentation_length_plot:
    input:
        matches=PROCESS+"blastn/Filtered_Annotated_"+str(FRAG)+"_VectorMatches_{sample}.blastn"
    params: 
        buffer=3*FRAG,
        threshold=config["MinInsertionLength"]
    log:
    	log=PROCESS+"log/qc/detailed_fragmentation_length_plot/{sample}.log"
    output:
        outpath=directory(FINAL+"qc/Fragmentation/Longest_Interval/{sample}/")
    run:
        shell("mkdir -p {output.outpath}")
        vhf.find_and_plot_longest_blast_interval(input.matches, params.buffer, params.threshold, output.outpath,log.log)
