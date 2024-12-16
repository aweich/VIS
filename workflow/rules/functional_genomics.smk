#Functional genomics rules for VIS

#sort insertion file
rule sort_insertion_file:
	input:
		f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed"
	log:
		log=f"{outdir}/intermediate/log/functional_genomics/sort_insertion_file/{{sample}}.log"
	output:
		f"{outdir}/intermediate/localization/Sorted_ExactInsertions_{{sample}}.bed"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		sort -k1,1 -k2,2n {input} > {output}
		) > {log.log} 2>&1
		"""

### Distance of VIS to genes, TSS, miRNAs, and TF: Annotation of the insertions
rule calc_distance_to_elements:
    input:
        insertions=f"{outdir}/intermediate/localization/Sorted_ExactInsertions_{{sample}}.bed",
        ref=config["ref_genome_ctrl"]
    params:
        genes=config["ucsc_Genes"],
        tf=config.get("ucsc_TF",""),
        tss=config.get("encode_hic", ""), 
        mirna=config.get("cosmic_genes", ""),
        exons=config.get("ucsc_exons", "")
    log:
        log=f"{outdir}/intermediate/log/functional_genomics/calc_distance_to_elements/{{sample}}.log"
    output:
        f"{outdir}/final/functional_genomics/Functional_distances_to_Insertions_{{sample}}.bed"
    run:
        vhf.calculate_element_distance(input.insertions, output[0], log.log, params[0:5])

rule plot_distance_to_elements:
	input:
		distancetable=f"{outdir}/final/functional_genomics/Functional_distances_to_Insertions_{{sample}}.bed"
	params:
		distances=list(range(-50000, 50001, 5000)),
		threshold=50000
	log:
        	log1=f"{outdir}/intermediate/log/functional_genomics/plot_distance_to_elements/scatter_{{sample}}.log",
        	log2=f"{outdir}/intermediate/log/functional_genomics/plot_distance_to_elements/violin_{{sample}}.log"
	output:
		scatter=report(f"{outdir}/final/functional_genomics/Plot_Distance_to_Genes_{fragmentsize}_{{sample}}.png"),
		violin=report(f"{outdir}/final/functional_genomics/BarPlot_Distance_to_Genes_{fragmentsize}_{{sample}}.png"),
	run:
		vhf.plot_element_distance(input.distancetable, params.distances, params.threshold, output.scatter, log.log1)
		vhf.plot_element_distance_violin(input.distancetable, params.distances, params.threshold, output.violin, log.log2)

rule plot_scoring:
    input:
        f"{outdir}/final/functional_genomics/Functional_distances_to_Insertions_{{sample}}.bed"
    log:
        log=f"{outdir}/intermediate/log/functional_genomics/plot_scoring/{{sample}}.log"
    output:
        plot=report(f"{outdir}/final/functional_genomics/Insertion_Scoring_{{sample}}.png")
    run:
        vhf.scoring_insertions(input[0], output.plot, log.log)

