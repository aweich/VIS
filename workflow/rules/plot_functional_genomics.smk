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
	        vhf_fa.plot_element_distance(input.distancetable, params.distances, params.threshold, output.scatter, log.log1)
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
            vhf_fa.scoring_insertions(input[0], output.plot, output.data, log.log)
        except Exception as e:
            with open(log.log, "a") as log_file:
                log_file.write(f"Error: {str(e)}\n")