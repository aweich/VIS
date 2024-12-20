# Collection of rules for analysis-specific plotting or pre-processing

# for the insertion specific gviz plot
#this ensures that a writable R lib exists where the needed packages will be installed to!
shell.prefix('export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.0 && ')      
rule chromosome_read_plots: 
	input:
		bed=f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed",
		#bam=f"{outdir}/intermediate/mapping/NoVectorAlignments_Postcut_{sample}_sorted.bam", #if normal bam would work here, we could cut several upstream rules!
		bam=f"{outdir}/intermediate/mapping/Postcut_{{sample}}_sorted.bam",
		H3K4Me1=config["ucsc_H3K4Me1"],
		H3K4Me3=config["ucsc_H3K4Me3"],
		H3K27Ac=config["ucsc_H3K27Ac"],
		gtf=config["ucsc_Genes_gtf"],
		#TF=config["ucsc_TF"],	-iTF {input.TF}
	log:
		log=f"{outdir}/intermediate/log/functional_genomics/chromosome_read_plots/{{sample}}.log"
	output:
		outpath=directory(f"{outdir}/final/functional_genomics/localization/{fragmentsize}_{{sample}}")
	params:
		buffer=50000
	shell: 
		r"""
		(
		mkdir {output.outpath}	#required, otherwise snakemake doesn't find the output folder and reports missing output
		Src/BAM_Inspection.R -ibam {input.bam} -ibed {input.bed} -iH3K4Me1 {input.H3K4Me1} -iH3K4Me3 {input.H3K4Me3} -iH3K27Ac {input.H3K27Ac} -igtf {input.gtf} -buffer {params.buffer} -o {output.outpath}   
		) > {log.log} 2>&1
		"""	


#for the circular insertion maps
rule blast_to_gff:
	input:
		ref=f"{outdir}/intermediate/blastn/humanref/Annotated_{fragmentsize}_InsertionMatches_{{sample}}.blastn",
		insertion=f"{outdir}/intermediate/blastn/{fragmentsize}_InsertionMatches_{{sample}}.blastn"
	output:
		insertionr=f"{outdir}/intermediate/blastn/{fragmentsize}_InsertionMatches_{{sample}}.gff",
		ref=f"{outdir}/intermediate/blastn/humanref/{fragmentsize}_InsertionMatches_{{sample}}.gff"  
	run:
		vhf.blast2gff(input.ref, output.ref)
		vhf.blast2gff(input.insertion, output.insertion)  
