#Functional genomics rules for VIS

#sort insertion file
rule sort_insertion_file:
	input:
		PROCESS+"localization/ExactInsertions_{sample}.bed"
	output:
		PROCESS+"localization/Sorted_ExactInsertions_{sample}.bed"
	run:
		shell("sort -k1,1 -k2,2n {input} > {output}")

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
	output:
		outpath=directory(FINAL+"functional_genomics/localization/" + str(FRAG)+"_{sample}")
	params:
		buffer=50000
	shell: 
		r"""
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -iH3K4Me1 {input.H3K4Me1} -iH3K4Me3 {input.H3K4Me3} -iH3K27Ac {input.H3K27Ac} -igtf {input.gtf} -buffer {params.buffer} -o {output.outpath}   
		"""	

### Distance of VIS to genes, TSS, miRNAs, and TF
rule distance_to_regulation: #this rule is strictly boudn to the 6 field BED format: Chr-Start-Stop-Name-Score-Strand; consider using "awk -v OFS="\t" {print $1, $2, $3, $4, ".", "."}'" to fill the columns 
    input:
        insertions=PROCESS+"localization/Sorted_ExactInsertions_{sample}.bed",
        ref=config["ref_genome_ctrl"]
    params:
        genes=config["ucsc_Genes"],
        tf=config.get("ucsc_TF", ""),
        tss=config.get("encode_hic", ""), 
        mirna=config.get("cosmic_genes", ""),
        exons=config.get("ucsc_exons", "")
    output:
        FINAL+"functional_genomics/Functional_distances_to_Insertions_{sample}.bed"
    run:
        # this command helps to balance if any of the other lists are not there. 
        valid_inputs = [params.genes, params.tf, params.tss, params.mirna, params.exons]
        valid_inputs = [file for file in valid_inputs if file]  # Exclude empty parameters

        # Build bedtools command dynamically
        input_files = " ".join(valid_inputs)
        if not input_files:
            raise ValueError("No valid input files provided for bedtools closest. Keep in mind that at least 'UCSC_genes' needs to be provided for this kind of analysis.")
        shell(f"""
        bedtools closest -a {input.insertions} -b {input_files} -D a -filenames | cut -f 1,2,3,4,7,11,14 > {output}
        """)

rule plot_distance_to_elements:
	input:
		distancetable=FINAL+"functional_genomics/Functional_distances_to_Insertions_{sample}.bed"
	params:
		distances=list(range(-50000, 50001, 5000)),
		threshold=50000
	output:
		genes=report(FINAL+"functional_genomics/Plot_Distance_to_Genes_" + str(FRAG)+"_{sample}.png"),
	run:
		vhf.plot_element_distance(input.distancetable, params.distances, params.threshold, output.genes)

rule plot_scoring:
    input:
        FINAL+"functional_genomics/Functional_distances_to_Insertions_{sample}.bed"
    output:
        plot=report(FINAL+"functional_genomics/Insertion_Scoring_{sample}.png")
    run:
        # Call your scoring function
        vhf.scoring_insertions(input[0], output.plot)
