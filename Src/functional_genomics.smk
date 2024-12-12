#Functional genomics rules for VIS

#sort insertion file
rule sort_insertion_file:
	input:
		PROCESS+"localization/ExactInsertions_{sample}.bed"
	log:
		log=PROCESS+"log/functional_genomics/sort_insertion_file/{sample}.log"
	output:
		PROCESS+"localization/Sorted_ExactInsertions_{sample}.bed"
	shell:
		"""
		(
		sort -k1,1 -k2,2n {input} > {output}
		) > {log.log} 2>&1
		"""

#Localization
#this ensures that a writable R lib exists where the needed packages will be installed to!
shell.prefix('export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.0 && ')      
rule chromosome_read_plots: 
	input:
		bed=PROCESS+"localization/ExactInsertions_{sample}.bed",
		#bam=PROCESS+"mapping/NoVectorAlignments_Postcut_{sample}_sorted.bam", #if normal bam would work here, we could cut several upstream rules!
		bam=PROCESS+"mapping/Postcut_{sample}_sorted.bam",
		H3K4Me1=config["ucsc_H3K4Me1"],
		H3K4Me3=config["ucsc_H3K4Me3"],
		H3K27Ac=config["ucsc_H3K27Ac"],
		gtf=config["ucsc_Genes_gtf"],
		#TF=config["ucsc_TF"],	-iTF {input.TF}
	log:
		log=PROCESS+"log/functional_genomics/chromosome_read_plots/{sample}.log"
	output:
		outpath=directory(FINAL+"functional_genomics/localization/" + str(FRAG)+"_{sample}")
	params:
		buffer=50000
	shell: 
		r"""
		(
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -iH3K4Me1 {input.H3K4Me1} -iH3K4Me3 {input.H3K4Me3} -iH3K27Ac {input.H3K27Ac} -igtf {input.gtf} -buffer {params.buffer} -o {output.outpath}   
		) > {log.log} 2>&1
		"""	

### Distance of VIS to genes, TSS, miRNAs, and TF: Annotation of the insertions
rule calc_distance_to_elements:
    input:
        insertions=PROCESS + "localization/Sorted_ExactInsertions_{sample}.bed",
        ref=config["ref_genome_ctrl"]
    params:
        genes=config["ucsc_Genes"],
        tf=config.get("ucsc_TF",""),
        tss=config.get("encode_hic", ""), 
        mirna=config.get("cosmic_genes", ""),
        exons=config.get("ucsc_exons", "")
    log:
        log=PROCESS+"log/functional_genomics/calc_distance_to_elements/{sample}.log"
    output:
        FINAL+"functional_genomics/Functional_distances_to_Insertions_{sample}.bed"
    run:
        vhf.calculate_element_distance(input.insertions, output[0], log.log, params[0:5])

rule plot_distance_to_elements:
	input:
		distancetable=FINAL+"functional_genomics/Functional_distances_to_Insertions_{sample}.bed"
	params:
		distances=list(range(-50000, 50001, 5000)),
		threshold=50000
	log:
        	log1=PROCESS+"log/functional_genomics/plot_distance_to_elements/scatter_{sample}.log",
        	log2=PROCESS+"log/functional_genomics/plot_distance_to_elements/violin_{sample}.log"
	output:
		scatter=report(FINAL+"functional_genomics/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png"),
		violin=report(FINAL+"functional_genomics/BarPlot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png"),
	run:
		vhf.plot_element_distance(input.distancetable, params.distances, params.threshold, output.scatter, log.log1)
		vhf.plot_element_distance_violin(input.distancetable, params.distances, params.threshold, output.violin, log.log2)

rule plot_scoring:
    input:
        FINAL+"functional_genomics/Functional_distances_to_Insertions_{sample}.bed"
    log:
        log=PROCESS+"log/functional_genomics/plot_scoring/{sample}.log"
    output:
        plot=report(FINAL+"functional_genomics/Insertion_Scoring_{sample}.png")
    run:
        vhf.scoring_insertions(input[0], output.plot, log.log)

