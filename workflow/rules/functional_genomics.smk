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

rule calc_distance_to_elements:
    input:
        #insertions=f"{outdir}/intermediate/localization/Sorted_ExactInsertions_{{sample}}.bed",
        insertions=f"{outdir}/intermediate/localization/temp_Sorted_ExactInsertions_{{sample}}.bed",
        ref=config["ref_genome_ctrl"]
    params:
        annotation_files={k: v for k, v in config.items() if k.startswith("annotate_")}
    log:
        log=f"{outdir}/intermediate/log/functional_genomics/calc_distance_to_elements/{{sample}}.log"
    output:
        f"{outdir}/final/functional_genomics/Functional_distances_to_Insertions_{{sample}}.bed"
    run:
        vhf.calculate_element_distance(input.insertions, output[0], log.log, params.annotation_files)
              
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

#kinda 'hacky' solution for flexibility: Now it does not matter if one annotation or 20 are defined. 
rule annotation_overlap_insertion:
    input:
        insertions_bed=f"{outdir}/intermediate/localization/temp_Sorted_ExactInsertions_{{sample}}.bed"
    params:
        annotation_files={k.replace("annotate_", ""): v for k, v in config.items() if k.startswith("annotate_")}
    log:
        log=f"{outdir}/intermediate/log/functional_genomics/annotation_overlap_insertion/{{sample}}.log"
    output:
        **{k.replace("annotate_", ""): f"{outdir}/intermediate/localization/annotation/Annotation_{k.replace('annotate_', '')}_Insertions_{{sample}}.bed"
           for k in config if k.startswith("annotate_")}
    run:
        vhf.run_bedtools_intersect(input.insertions_bed, output, log.log, params.annotation_files)

