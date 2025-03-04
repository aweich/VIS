######
######
###### Base modifications of reads with isnertions
######
######

# Notes: Methylartist works in principle, but there is one very interesting problem: SInce the Insertion are usually cut out from the optimal alignment, their modifications cannot be seen anymore.
# 		This means that it is not possible to have an idea about what is going on at the exact insertiopn if we use the common tools! -> this makes manual labor necessary

rule readnames_final:
	input:
		final_insertions=f"{outdir}/final/localization/ExactInsertions_{{sample}}.bed"
	log:
		log=f"{outdir}/intermediate/log/base_modifications/readnames_final/{{sample}}.log"
	output:
		readnames=temp(f"{outdir}/intermediate/base_modifications/final_insertion_readnames_{{sample}}.txt")
	conda:
		"../envs/VIS_dummy_env.yml"
	shell:
		"""
		(
		cut {input} -f 4 > {output}
		) > {log.log} 2>&1
		"""

rule only_keep_filtered_insertions:
	input:
		isobam=f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam",
		readnames=f"{outdir}/intermediate/cmod/final_insertion_readnames_{{sample}}.txt"
	log:
		log=f"{outdir}/intermediate/log/base_modifications/only_keep_filtered_insertions/{{sample}}.log"
	output:
		bam=f"{outdir}/intermediate/base_modifications/Final_Isolated_Reads_{{sample}}.bam"
	conda:
		"../envs/VIS_samtools_env.yml"
	shell:
		"""
		(
		samtools view -h -b -N {input.readnames} {input.isobam} | samtools sort > {output.bam}
		samtools index {output.bam}
		) > {log.log} 2>&1
		"""

rule modkit:
	input:
		isobam=f"{outdir}/intermediate/base_modifications/Final_Isolated_Reads_{{sample}}.bam"
	log:
		log=f"{outdir}/intermediate/log/base_modifications/modkit/{{sample}}.log"
	output:
		tsv=f"{outdir}/final/base_modifications/Isolated_Reads_{{sample}}.tsv"
	conda:
		"../envs/VIS_modkit_env.yml"
	threads: config["threads"]
	shell:
		"""
		(
		modkit extract full -t {threads} {input.isobam} {output.tsv}
		) > {log.log} 2>&1
		"""

rule call_modkit:
	input:
		isobam=f"{outdir}/intermediate/base_modifications/Final_Isolated_Reads_{{sample}}.bam"
	log:
		log=f"{outdir}/intermediate/log/base_modifications/call_modkit/{{sample}}.log"
	output:
		tsv=f"{outdir}/final/base_modifications/Calls_Isolated_Reads_{{sample}}.tsv"
	conda:
		"../envs/VIS_modkit_env.yml"
	threads: config["threads"]
	shell:
		"""
		(
		modkit extract calls -t {threads} {input.isobam} {output.tsv}
		) > {log.log} 2>&1
		"""

'''
rule specific_methylartist:
	input:
		bam=f"{outdir}/intermediate/mapping/Precut_{{sample}}_sorted.bam",
		ref=f"{outdir}/intermediate/mapping/insertion_ref_genome.fa" #must be indexed!
	log:
		log=f"{outdir}/intermediate/log/cmod/specific_methylartist/{{sample}}.log"
	output:
		plot=f"{outdir}/intermediate/cmod/Region3_Reads_{{sample}}.png"
	conda:
		"../envs/VIS_methylartist_env.yml"
	shell:
		"""
		(
		methylartist locus -b {input.bam} -i chr17:31124037-31180287 -n C -r {input.ref} -l 31154037-31159287 -o {output.plot}
		) > {log.log} 2>&1
		"""


rule inserted_seq:
	input:
		insertion=f"{outdir}/intermediate/fasta/Isolated_Reads_{{sample}}.fa",
		coordinates=f"{outdir}/intermediate/blastn/Coordinates_{fragmentsize}_InsertionMatches_{{sample}}.blastn"
	log:
		log=f"{outdir}/intermediate/log/cmod/inserted_seq/{{sample}}.log"
	output:
		fasta=f"{outdir}/intermediate/fasta/Inserted_sequence_{{sample}}.fa"
	run:
		try:
			vhf.get_inserted_fasta_seq(input.insertion, input.coordinates, output[0], log.log)
		except Exception as e:
			with open(log.log, "a") as log_file:
					log_file.write(f"Error: {str(e)}\n")
'''