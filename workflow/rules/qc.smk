######
######
###### Quality control and normalization
######
######
 
#FRAG=config["fragment_size"]

#full stats of input data		
rule nanoplot:
	input:
		f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam"
	output:
		f"{outdir}/intermediate/qc/nanoplot/{{sample}}/NanoStats.txt"
	log:
		log=f"{outdir}/intermediate/log/qc/nanoplot/{{sample}}.log"
	params:
		outdir=directory(f"{outdir}/intermediate/qc/nanoplot/{{sample}}/"),
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
		bam=f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam"
	log:
		log=f"{outdir}/intermediate/log/qc/bam_coverage/{{sample}}.log"
	output:
		covbed=f"{outdir}/intermediate/qc/Coverage/Genomecoverage_{{sample}}.bed"
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
    	bam=f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam",
        readnames=f"{outdir}/intermediate/blastn/Readnames_{fragmentsize}_InsertionMatches_{{sample}}.txt"
    params:
        tempdir=f"{outdir}/intermediate/temp_fastq_{{sample}}"  # Temporary directory for intermediate files
    log:
    	log=f"{outdir}/intermediate/log/qc/extract_fastq_insertions/{{sample}}.log"
    output:
        fastq=f"{outdir}/intermediate/qc/fastqc/{{sample}}_filtered.fastq"
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
        f"{outdir}/intermediate/qc/fastqc/{{sample}}_filtered.fastq"
    params:
        prefix=f"{outdir}/intermediate/qc/fastqc/readlevel_{{sample}}/{{sample}}_read_"
    log:
    	log=f"{outdir}/intermediate/log/qc/read_level_fastqc/{{sample}}.log"
    output:
        directory(f"{outdir}/intermediate/qc/fastqc/readlevel_{{sample}}/")
    conda:
    	"../envs/VIS_fastqc_env.yml"
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
        fastqc=expand(f"{outdir}/intermediate/qc/fastqc/readlevel_{{sample}}/", sample=SAMPLES),
        nanoplot=expand(f"{outdir}/intermediate/qc/nanoplot/{{sample}}/NanoStats.txt", sample=SAMPLES)
    params:
        qc_output_dir=lambda wildcards, input: f"{outdir}/intermediate/qc/", #lambda added to please the linter
        qc_report_location=lambda wildcards, input: f"{outdir}/final/qc/", #lambda added to please the linter
    log:
    	log=f"{outdir}/intermediate/log/qc/multiqc/out.log"
    output: 
        intermediate=f"{outdir}/intermediate/qc/multiqc_report.html",
        final=report(f"{outdir}/final/qc/multiqc_report.html")
    conda:
    	"../envs/VIS_multiqc_env.yml"
    shell:
        """
        (
        multiqc {input.fastqc} --dirs {input.nanoplot} --force -o {params.qc_output_dir}
        
        # copy report to final location
        cp {output.intermediate} {output.final}
        ) > {log.log} 2>&1
        """

### Read level overview of mapping quality before and after the cut out of the insertions
rule extract_mapping_quality:
    input:
        bam=f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam",
        bam2=f"{outdir}/intermediate/mapping/Postcut_{{sample}}_sorted.bam",
        bam3=f"{outdir}/intermediate/mapping/Postcut_{{sample}}_unfiltered_sorted.bam",
        readnames=f"{outdir}/intermediate/blastn/Readnames_{fragmentsize}_InsertionMatches_{{sample}}.txt"
    params:
        tempdir=f"{outdir}/intermediate/temp_mapping_{{sample}}"
    log:
    	log=f"{outdir}/intermediate/log/qc/extract_mapping_quality/{{sample}}.log"
    output:
        quality_scores=f"{outdir}/intermediate/qc/{{sample}}_precut_mapping_quality.txt",
        quality_scores2=f"{outdir}/intermediate/qc/{{sample}}_postcut_mapping_quality.txt",
        quality_scores3=f"{outdir}/intermediate/qc/{{sample}}_postcut_unfiltered_mapping_quality.txt"
    conda:
    	"../envs/VIS_samtools_env.yml"
    shell:
        '''
        (
        mkdir -p {params.tempdir}

        # Extract reads of interest
        samtools view {input.bam} | grep -wF -f {input.readnames} > {params.tempdir}/temp_precut_reads.sam
	samtools view {input.bam2} | grep 'Insertion' > {params.tempdir}/temp_postcut_reads.sam
	samtools view {input.bam3} | grep 'Insertion' > {params.tempdir}/temp_postcut_unfiltered_reads.sam

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
        quality_scores_pre=f"{outdir}/intermediate/qc/{{sample}}_precut_mapping_quality.txt",
        quality_scores_filtered=f"{outdir}/intermediate/qc/{{sample}}_postcut_unfiltered_mapping_quality.txt",
        quality_scores_post=f"{outdir}/intermediate/qc/{{sample}}_postcut_mapping_quality.txt"
    params:
        prefixes=["Precut", "Postcut", "Postcut_filtered"]
    log:
    	log=f"{outdir}/intermediate/log/qc/finalize_mapping_quality/{{sample}}.log"
    output:
        outfile=f"{outdir}/final/qc/mapq/Insertions_{{sample}}_mapq.txt"
    run:
    	try:
    	    vhf.join_read_mapq(input[0:3], params.prefixes, output.outfile, log.log)
    	except Exception as e:
            with open(log.log, "a") as log_file:
                log_file.write(f"Error: {str(e)}\n")

rule generate_mapq_heatmap:
    input:
        table=f"{outdir}/final/qc/mapq/Insertions_{{sample}}_mapq.txt"
    log:
    	log=f"{outdir}/intermediate/log/qc/generate_mapq_heatmap/{{sample}}.log"
    output:
        heatmap=report(f"{outdir}/final/qc/mapq/{{sample}}_mapq_heatmap_image.png")
    run:
        try:
            vhf.plot_mapq_changes(input.table, output.heatmap, log.log)
        except Exception as e:
            with open(log.log, "a") as log_file:
                log_file.write(f"Error: {str(e)}\n")

#### Visualize fragmentation and longest consecutive interval per read

rule fragmentation_distribution_plots:
	input:
		f"{outdir}/intermediate/blastn/Filtered_Annotated_{fragmentsize}_InsertionMatches_{{sample}}.blastn",
		f"{outdir}/intermediate/blastn/ref/Filtered_Annotated_{fragmentsize}_InsertionMatches_{{sample}}.blastn"
	params:
		fragmentsize
	log:
		log1=f"{outdir}/intermediate/log/qc/fragmentation_distribution_plots/fragmentation_match_distribution_{{sample}}.log",
		log2=f"{outdir}/intermediate/log/qc/fragmentation_distribution_plots/fragmentation_read_match_distribution_{{sample}}.log",
		log3=f"{outdir}/intermediate/log/qc/fragmentation_distribution_plots/fragmentation_match_distribution_{{sample}}.log",
		log4=f"{outdir}/intermediate/log/qc/fragmentation_distribution_plots/fragmentation_read_match_distribution_{{sample}}.log"
	output:
		outpath=directory(f"{outdir}/final/qc/Fragmentation/Insertions/insertions_{fragmentsize}_{{sample}}"),
		outpath2=directory(f"{outdir}/final/qc/Fragmentation/Reference/reference_{fragmentsize}_{{sample}}")
	run:
		shell("mkdir {output.outpath}")
		try:
		    vhf.fragmentation_match_distribution(input[0], params[0], output[0], log.log1)
		except Exception as e:
		    with open(log.log1, "a") as log_file:
                        log_file.write(f"Error: {str(e)}\n")
              
		try:
		    vhf.fragmentation_read_match_distribution(input[0], params[0], output[0], log.log2)
		except Exception as e:
		    with open(log.log2, "a") as log_file:
                        log_file.write(f"Error: {str(e)}\n")
                
		shell("mkdir {output.outpath2}")
		
		try:
		    vhf.fragmentation_match_distribution(input[1], params[0], output[1], log.log3)
		except Exception as e:
		    with open(log.log3, "a") as log_file:
                        log_file.write(f"Error: {str(e)}\n")
                
		try:
		    vhf.fragmentation_read_match_distribution(input[1], params[0], output[1], log.log4)
		except Exception as e:
		    with open(log.log4, "a") as log_file:
                        log_file.write(f"Error: {str(e)}\n")
                
rule detailed_fragmentation_length_plot:
    input:
        matches=f"{outdir}/intermediate/blastn/Filtered_Annotated_{fragmentsize}_InsertionMatches_{{sample}}.blastn"
    params: 
        bridge=config["bridging_size"],
        threshold=config["MinInsertionLength"]
    log:
    	log=f"{outdir}/intermediate/log/qc/detailed_fragmentation_length_plot/{{sample}}.log"
    output:
        outpath=directory(f"{outdir}/final/qc/Fragmentation/Longest_Interval/{{sample}}/")
    run:
        shell("mkdir -p {output.outpath}")
        try:
            vhf.find_and_plot_longest_blast_interval(input.matches, params.bridge, params.threshold, output.outpath,log.log)
        except Exception as e:
            with open(log.log, "a") as log_file:
                log_file.write(f"Error: {str(e)}\n")
