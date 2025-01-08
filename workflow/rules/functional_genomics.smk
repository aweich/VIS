#Functional genomics rules for VIS

#sort insertion file
rule sort_insertion_file:
	input:
		f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed"
	log:
		log=f"{outdir}/intermediate/log/functional_genomics/sort_insertion_file/{{sample}}.log"
	output:
		temp(f"{outdir}/intermediate/localization/Sorted_ExactInsertions_{{sample}}.bed")
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
        genes=config["annotation_1"],
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
		distances=list(range(-10000, 10001, 2000)),
		threshold=10000
	log:
        	log1=f"{outdir}/intermediate/log/functional_genomics/plot_distance_to_elements/scatter_{{sample}}.log",
        	log2=f"{outdir}/intermediate/log/functional_genomics/plot_distance_to_elements/violin_{{sample}}.log"
	output:
		scatter=report(f"{outdir}/final/functional_genomics/Plot_Distance_to_Genes_{fragmentsize}_{{sample}}.png"),
	run:
	    try:
	        vhf.plot_element_distance(input.distancetable, params.distances, params.threshold, output.scatter, log.log1)
	    except Exception as e:
	        with open(log.log1, "a") as log_file:
                    log_file.write(f"Error: {str(e)}\n")
