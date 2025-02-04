#Functional genomics rules for VIS

#sort insertion file
rule sort_insertion_file:
	input:
		exact=f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed",
		point=f"{outdir}/intermediate/localization/annotation/temp_Insertions_{{sample}}.bed"
	log:
		log=f"{outdir}/intermediate/log/functional_genomics/sort_insertion_file/{{sample}}.log",
		log2=f"{outdir}/intermediate/log/functional_genomics/sort_insertion_file/temp_{{sample}}.log"
	output:
		exact=temp(f"{outdir}/intermediate/localization/Sorted_ExactInsertions_{{sample}}.bed"),
		point=f"{outdir}/intermediate/localization/temp_Sorted_ExactInsertions_{{sample}}.bed"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		sort -k1,1 -k2,2n {input.exact} > {output.exact}
		) > {log.log} 2>&1
		(
		sort -k1,1 -k2,2n {input.point} > {output.point}
		) > {log.log2} 2>&1
		"""

### Distance of VIS to genes, TSS, miRNAs, and TF: Annotation of the insertions
rule calc_distance_to_elements:
    input:
        #insertions=f"{outdir}/intermediate/localization/Sorted_ExactInsertions_{{sample}}.bed",
        insertions=f"{outdir}/intermediate/localization/temp_Sorted_ExactInsertions_{{sample}}.bed",
        ref=config["ref_genome_ctrl"]
    params:
        genes=config["ucsc_introns"],
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
		#violin=report(f"{outdir}/final/functional_genomics/BarPlot_Distance_to_Genes_{fragmentsize}_{{sample}}.png"),
	run:
	    try:
	        vhf.plot_element_distance(input.distancetable, params.distances, params.threshold, output.scatter, log.log1)
		#vhf.plot_element_distance_violin(input.distancetable, params.distances, params.threshold, output.violin, log.log2)
	    except Exception as e:
	        with open(log.log1, "a") as log_file:
                    log_file.write(f"Error: {str(e)}\n")
                
rule plot_scoring:
    input:
        f"{outdir}/final/functional_genomics/Functional_distances_to_Insertions_{{sample}}.bed"
    log:
        log=f"{outdir}/intermediate/log/functional_genomics/plot_scoring/{{sample}}.log"
    output:
        plot=report(f"{outdir}/final/functional_genomics/Insertion_Scoring_{{sample}}.svg"),
        data=(f"{outdir}/final/functional_genomics/Insertion_Scoring_Data_{{sample}}.txt")
    run:
        try:
            vhf.scoring_insertions(input[0], output.plot, output.data, log.log)
        except Exception as e:
            with open(log.log, "a") as log_file:
                log_file.write(f"Error: {str(e)}\n")
###
rule modify_insertions:
	input:
		f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed"
	output:
		temp(f"{outdir}/intermediate/localization/annotation/temp_Insertions_{{sample}}.bed")
	log:
		f"{outdir}/intermediate/log/functional_genomics/modify_insertions/{{sample}}.log"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		awk '{{OFS="\t"; $3 = $2 + 1; print $0}}' {input} > {output}
		) > {log} 2>&1
		"""
		        
rule annotation_overlap_insertion:
	input:
		f"{outdir}/intermediate/localization/temp_Sorted_ExactInsertions_{{sample}}.bed"
		#f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed"
	params:
        	exons=config.get("ucsc_exons"),
        	introns=config.get("ucsc_introns"),
        	promoter=config.get("ucsc_promoter"),
        	gene=config.get("annotation_1"),
	log:
		log=f"{outdir}/intermediate/log/functional_genomics/annotation_overlap_insertion/introns_{{sample}}.log",
		log2=f"{outdir}/intermediate/log/functional_genomics/annotation_overlap_insertion/exons_{{sample}}.log",
		log3=f"{outdir}/intermediate/log/functional_genomics/annotation_overlap_insertion/promoter_{{sample}}.log",
		log4=f"{outdir}/intermediate/log/functional_genomics/annotation_overlap_insertion/gene_{{sample}}.log"
	output:
		introns=f"{outdir}/intermediate/localization/annotation/Annotation_introns_Insertions_{{sample}}.bed",
		exons=f"{outdir}/intermediate/localization/annotation/Annotation_exons_Insertions_{{sample}}.bed",
		promoter=f"{outdir}/intermediate/localization/annotation/Annotation_promoter_Insertions_{{sample}}.bed",
		gene=f"{outdir}/intermediate/localization/annotation/Annotation_gene_Insertions_{{sample}}.bed"
	conda:
		"../envs/VIS_bedtools_env.yml"
	shell: #1 base is enough
		"""
		(
		bedtools intersect -a {input} -b {params.introns} -wb > {output.introns}
		) > {log.log} 2>&1
		(
		bedtools intersect -a {input} -b {params.exons} -wb > {output.exons}
		) > {log.log2} 2>&1
		(
		bedtools intersect -a {input} -b {params.promoter} -wb > {output.promoter}
		) > {log.log3} 2>&1
		(
		bedtools intersect -a {input} -b {params.gene} -wb > {output.gene}
		) > {log.log4} 2>&1
		"""
