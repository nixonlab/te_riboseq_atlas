#! /usr/bin/env python
# -*- coding utf-8 -*-

################################## HERV QUANT ##################################

rule alignment_filtering:
    conda:
        "../envs/samtools.yaml"
    input: "results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.Aligned.out.sam"
    output: "results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.Aligned.out.filtered.bam"
    params:
        sam_file_filtered = "results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.Aligned.out.filtered.sam"
    shell:
        """
        sed '/uc.*/d' {input} > {params.sam_file_filtered}
        samtools view -bS {params.sam_file_filtered} > {output}
        rm {params.sam_file_filtered}
        """

rule salmon:
    conda: "../envs/salmon.yaml"
    input: "results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.Aligned.out.filtered.bam"
    output: "results/salmon/{samid}/quant.sf"
    params:
        hervquant_ref = "databases/remotefiles/hervquant_final_reference.fa",
        out_dir =  "results/salmon/{samid}",
        out_file = "results/salmon/{samid}/{samid}_quant.sf"
    threads: 20
    shell:
        """
        salmon quant \
          -t {params.hervquant_ref} \
          -l A \
          -a {input} \
          -o {params.out_dir} \
          -p {threads}

        cp {output} {params.out_file}
        """

rule salmon_complete:
    input:
        rules.salmon.output
    output:
        touch("results/salmon/{samid}_salmon_completed.txt")
