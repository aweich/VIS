#VIS detection rules

######
######
###### Insertion BAM: Convert to fasta 
######
######

rule copy_config_version:
	input:
		"config/config.yml"
	log:
		log=f"{outdir}/intermediate/log/detection/copy_config_version/out.log"
	output:
		f"{outdir}/config_settings.yml"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		cp {input} {output} 
		) > {log.log} 2>&1
		"""

rule build_insertion_reference:
	input:
		ref=config["ref_genome_ctrl"],
		vector=config["vector_fasta"]
	log:
		log=f"{outdir}/intermediate/log/detection/build_insertion_reference/out.log"
	output:
		f"{outdir}/intermediate/mapping/vector_ref_genome.fa"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		cat {input.ref} {input.vector} > {output}
		) > {log.log} 2>&1
		"""
		

rule minimap_index:
	input:
		ref=f"{outdir}/intermediate/mapping/vector_ref_genome.fa"
	log:
		log=f"{outdir}/intermediate/log/detection/minimap_index/out.log"
	output:
		index=temp(f"{outdir}/intermediate/mapping/ref_genome_index.mmi")
	resources:
		mem_mb=5000
	conda:
		"../envs/VIS_minimap_env.yml"
	shell:
		"""
		(
		minimap2 -d {output.index} {input.ref} 
		) > {log.log} 2>&1
		"""

rule make_fasta_without_tags: #fasta of raw data no trimming whatsoever
	input:
		fq=lambda wildcards: config["samples"][wildcards.sample]
	log:
		log=f"{outdir}/intermediate/log/detection/make_fasta_without_tags/{{sample}}.log"
	output:
		fasta=f"{outdir}/intermediate/fasta/Full_{{sample}}.fa"
	conda:
		"../envs/VIS_samtools_env.yml"
	shell: 
		"""
		(
		samtools fasta {input} -o {output} > {output}
		) > {log.log} 2>&1
		"""
######
######
###### "Clean" BAM: Cut-out fasta to BAM via Mapping to reference 
######
######
		
#rule to create the BAM files with the non-insertion reads and the splitted read fragments
rule Non_insertion_mapping: #mapping against the unaltered referenc egenome
	input:
		fasta=f"{outdir}/intermediate/fasta/Cleaved_{{sample}}_noVector.fa",
		genome=f"{outdir}/intermediate/mapping/vector_ref_genome.fa"
	output:
		f"{outdir}/intermediate/mapping/Postcut_{{sample}}_unfiltered_sorted.bam"
	log:
		log=f"{outdir}/intermediate/log/detection/Non_insertion_mapping/{{sample}}.log"
	resources:
		mem_mb=5000
	conda:
		"../envs/VIS_minimap_env.yml"
	shell: #N=0 instead of default N=1
		"""
		(
		minimap2 -y -ax map-ont --score-N 0 {input.genome} {input.fasta} | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		) > {log.log} 2>&1
		"""

rule insertion_mapping: #conserves tags!
	input:
		bam=lambda wildcards: config["samples"][wildcards.sample],
		minimapref=f"{outdir}/intermediate/mapping/ref_genome_index.mmi",
		ref=f"{outdir}/intermediate/mapping/vector_ref_genome.fa"
	output:
		f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam"
	log:
		log=f"{outdir}/intermediate/log/detection/insertion_mapping/{{sample}}.log"
	resources:
		mem_mb=5000
	conda:
		"../envs/VIS_minimap_env.yml"
	shell:
		"""
		(
		samtools bam2fq -T '*' {input.bam}| minimap2 -y -ax map-ont {input.minimapref} - | samtools sort |  samtools view -F 2304 -o {output}
		samtools index {output}
		) > {log.log} 2>&1
		"""

rule clean_postcut_by_maping_quality:
	input:
		f"{outdir}/intermediate/mapping/Postcut_{{sample}}_unfiltered_sorted.bam"
	params:
		mapq=config["MAPQ"]
	output:
		f"{outdir}/intermediate/mapping/Postcut_{{sample}}_sorted.bam"
	log:
		log=f"{outdir}/intermediate/log/detection/clean_postcut_by_maping_quality/{{sample}}.log"
	conda:
		"../envs/VIS_samtools_env.yml"
	shell:
		"""
		(
		samtools view -h -q {params.mapq} {input} -o {output}
		samtools index {output}
		) > {log.log} 2>&1
		"""
######
######
###### Genomic Coordinates of Reads with and without matches
######
######

rule BAM_to_BED:
	input:
		#get_input_names
		precut=f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam",
		postcut=f"{outdir}/intermediate/mapping/Postcut_{{sample}}_sorted.bam" #Reads that contained an insertion before, are now marked with "_Buffer/Insertion/0,1,2..."
	output:
		postcut=f"{outdir}/intermediate/mapping/Postcut_{{sample}}.bed",
		precut=f"{outdir}/intermediate/mapping/Precut_{{sample}}.bed"
	log:
		log1=f"{outdir}/intermediate/log/detection/BAM_to_BED/Precut_{{sample}}.log",
		log2=f"{outdir}/intermediate/log/detection/BAM_to_BED/Postcut_{{sample}}.log"
	conda:
		"../envs/VIS_bedtools_env.yml"
	shell:
		"""
		(
		bedtools bamtobed -cigar -i {input.precut} > {output.precut}
		) > {log.log1} 2>&1
		(
		bedtools bamtobed -cigar -i {input.postcut} > {output.postcut}
		) > {log.log2} 2>&1
		"""

######
######
###### fasta preparation: Cut out of blast-detected vector fragments and create new "cut-out" fasta 
######
######

rule get_cleavage_sites_for_fasta: #filters and combines matches
	input:
		f"{outdir}/intermediate/blastn/Filtered_Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn"
	params:
		filteroption=True,
		filtervalue=config["MinInsertionLength"], 
		overlap=3*fragmentsize #this is the distance of the start-stop that is allowed to exist to still be combined; This should not be lower than FRAG!
	log:
		log=f"{outdir}/intermediate/log/detection/get_cleavage_sites_for_fasta/{{sample}}.log"
	output:
		cleavage=f"{outdir}/intermediate/blastn/CleavageSites_{fragmentsize}_VectorMatches_{{sample}}.blastn",
		reads=f"{outdir}/intermediate/blastn/Readnames_{fragmentsize}_VectorMatches_{{sample}}.txt"
	run:
		vhf.splitting_borders(input[0],params.filteroption, params.filtervalue, params.overlap, output.cleavage, output.reads, log.log)

rule split_fasta:
	input:
		breakpoints=f"{outdir}/intermediate/blastn/CleavageSites_{fragmentsize}_VectorMatches_{{sample}}.blastn",
		fasta=f"{outdir}/intermediate/fasta/Full_{{sample}}.fa"
	params:
		mode=config["splitmode"] #if each split fasta substring should be used individually, use "Separated" Join, New mode: Buffer
	log:
		log=f"{outdir}/intermediate/log/detection/split_fasta_by_borders/{{sample}}.log"
	output:
		fasta=f"{outdir}/intermediate/fasta/Cleaved_{{sample}}_noVector.fa",
		vector=f"{outdir}/intermediate/fasta/Insertion_{{sample}}_Vector.fa"
	run:
		vhf.split_fasta_by_borders(input.breakpoints, input.fasta, params.mode, output.fasta, output.vector, log.log)

######
######
###### Vector preparation: Fragmentation 
######
######
rule prepare_vector:
	input:
		config["vector_fasta"] #vector fasta sequence
	log:
		log=f"{outdir}/intermediate/log/detection/prepare_vector/out.log"
	output: 
		fasta=f"{outdir}/intermediate/fasta/fragments/Forward_Backward_Vector.fa"
	run:
		vhf.reversevector(input[0], output[0], log.log)
		
rule vector_fragmentation:
	input:
		f"{outdir}/intermediate/fasta/fragments/Forward_Backward_Vector.fa" #does not change anything so it can be removed imo
	params:
		fragmentsize
	log:
		log=f"{outdir}/intermediate/log/detection/vector_fragmentation/out.log"
	output: 
		fasta=f"{outdir}/intermediate/fasta/fragments/{fragmentsize}_Vector_fragments.fa"
	run:
		vhf.fragmentation_fasta(input[0], params[0], output[0], log.log)

######
######
###### BLAST Searches - Against Vector and healthy human reference
######
######

rule make_blastn_DB:
	input:
		vector_fragmented = f"{outdir}/intermediate/fasta/fragments/{fragmentsize}_Vector_fragments.fa"
	log:
		log=f"{outdir}/intermediate/log/detection/make_blastn_DB/out.log"
	output:
		multiext(f"{outdir}/intermediate/fasta/fragments/{fragmentsize}_Vector_fragments.fa",
			".ndb",
			".nhr",
			".nin",
			".not",
			".nsq",
			".ntf",
			".nto"
		)
	conda:
		"../envs/VIS_blastn_env.yml"
	shell:
		"""
		(
		makeblastdb -in {input.vector_fragmented} -dbtype nucl -blastdb_version 5 
		) > {log.log} 2>&1
		"""

rule find_vector_BLASTn:
	input:
		fasta=f"{outdir}/intermediate/fasta/Full_{{sample}}.fa",
		dummy=f"{outdir}/intermediate/fasta/fragments/{fragmentsize}_Vector_fragments.fa.ndb", #provokes the building of the database first!
		vector=f"{outdir}/intermediate/fasta/fragments/{fragmentsize}_Vector_fragments.fa"
	params:
		tempdir=f"{outdir}/intermediate/temp_{{sample}}",
	log:
		log=f"{outdir}/intermediate/log/detection/find_vector_BLASTn/{{sample}}.log"
	output:
		f"{outdir}/intermediate/blastn/{fragmentsize}_VectorMatches_{{sample}}.blastn"
	conda:
		"../envs/VIS_blastn_env.yml"
	shell:
		"""
		(
	mkdir {params.tempdir}	
        
        blastn \
        -query {input.fasta} \
        -db {input.vector} \
        -out {params.tempdir}/temp_output.blastn \
        -evalue 1e-5 \
        -outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'

        # Filter results based on bitscore > 50
        awk '$16 > 50' {params.tempdir}/temp_output.blastn > {output}

        # Clean up temporary files
        rm -r {params.tempdir}
        	) > {log.log} 2>&1
        	"""
		
#blastn vector against human genome: Which vector parts are close to human sequences so that they might raise a false positivite BLAST match
rule find_vector_BLASTn_in_humanRef:
	input:
		fasta=f"{outdir}/intermediate/fasta/fragments/{fragmentsize}_Vector_fragments.fa"
	params:
		vector=config["blastn_db"] #full human reference
	log:
		log=f"{outdir}/intermediate/log/detection/find_vector_BLASTn_in_humanRef/{{sample}}.log"	
	output:
		temp(f"{outdir}/intermediate/blastn/humanref/{fragmentsize}_VectorMatches_{{sample}}.blastn")
	conda:
		"../envs/VIS_blastn_env.yml"
	shell:
		"""
		(
		blastn \
		-query {input} \
		-db {params.vector} \
		-out {output} \
		-evalue 1e-5 \
		-outfmt '6 qseqid sseqid qseq sseq qlen slen qstart qend sstart send length mismatch pident qcovs evalue bitscore'
		) > {log.log} 2>&1
		"""

rule hardcode_blast_header:
	input: 
		vector=f"{outdir}/intermediate/blastn/{fragmentsize}_VectorMatches_{{sample}}.blastn",
		humanref=f"{outdir}/intermediate/blastn/humanref/{fragmentsize}_VectorMatches_{{sample}}.blastn"
	log:
		log=f"{outdir}/intermediate/log/detection/hardcode_blast_header/{{sample}}.log"	
	output:
		vector=f"{outdir}/intermediate/blastn/Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn",
		humanref=f"{outdir}/intermediate/blastn/humanref/Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.vector} > {output.vector}
		echo -e 'QueryID\tSubjectID\tQueryAligned\tSubjectAligned\tQueryLength\tSubjectLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tLength\tMismatch\tPercentageIdentity\tQueryCov\tevalue\tbitscore' | cat - {input.humanref} > {output.humanref}
		) > {log.log} 2>&1
		"""
		
#filter BLAST matches by min length
rule extract_by_length:
	input:
		blast=f"{outdir}/intermediate/blastn/Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn",
		humanref=f"{outdir}/intermediate/blastn/humanref/Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn"
	params:
		threshold=config["MinLength"]
	log:
		log=f"{outdir}/intermediate/log/detection/extract_by_length/{{sample}}.log"	
	output:
		blast=f"{outdir}/intermediate/blastn/Filtered_Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn",
		humanref=f"{outdir}/intermediate/blastn/humanref/Filtered_Annotated_{fragmentsize}_VectorMatches_{{sample}}.blastn"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		awk -F'\t' '$11>={params.threshold}' {input.blast} > {output.blast}
		awk -F'\t' '$11>={params.threshold}' {input.humanref} > {output.humanref}
		) > {log.log} 2>&1
		""" 

######
######
###### Visualisation of insertions
######



rule basic_insertion_plots:
	input:
		expand(f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed", sample=SAMPLES)
	output:
		report(f"{outdir}/final/localization/Heatmap_Insertion_Chr.png"),
		report(f"{outdir}/final/localization/Insertion_length.png")
	log:
		log1=f"{outdir}/intermediate/log/detection/basic_insertion_plots/heat.log",
		log2=f"{outdir}/intermediate/log/detection/basic_insertion_plots/length.log"
	run:
		vhf.plot_bed_files_as_heatmap(input, output[0], log.log1)
		vhf.plot_insertion_length(input, output[1], log.log2)

######
######
###### Exact localization
######
######

#exact coordinates of the matching fragments

rule calculate_exact_insertion_coordinates:
	input:
		bed=f"{outdir}/intermediate/mapping/Postcut_{{sample}}.bed",
		borders=f"{outdir}/intermediate/blastn/CleavageSites_{fragmentsize}_VectorMatches_{{sample}}.blastn"
	params:
		mode=config["splitmode"]
	log:
		log=f"{outdir}/intermediate/log/detection/calculate_exact_insertion_coordinates/{{sample}}.log"
	output:
		out=f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed",
	run:
		vhf.reconstruct_coordinates(input.bed, input.borders, params.mode, output.out, log.log)

rule collect_outputs:
	input:
		coordinates=f"{outdir}/intermediate/localization/ExactInsertions_{{sample}}.bed",
	log:
		log=f"{outdir}/intermediate/log/detection/collect_outputs/{{sample}}.log"
	output:
		coordinates=f"{outdir}/final/localization/ExactInsertions_{{sample}}.bed"
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		cp {input.coordinates} {output.coordinates}
		) > {log.log} 2>&1
		"""
