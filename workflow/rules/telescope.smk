#! /usr/bin/env python
# -*- coding utf-8 -*-

################################### TELESCOPE ###################################

rule telescope:
    conda: "../envs/telescope.yaml"
    output:
        "results/telescope/{samid}/{samid}-telescope_report.tsv"
    input:
        bam = "results/star_alignment/{samid}/{samid}_GDC38.Aligned.out.bam",
        annotation = rules.repeat_annotation_retro38.output[0]
    benchmark: "benchmarks/telescope/{samid}_telescope.tsv"
    log:
        "results/telescope/{samid}/telescope.log"
    threads: config["telescope_threads"]
    params:
        outdir = "results/telescope/{samid}",
        exp_tag = "{samid}"
    shell:
        """
        telescope assign\
         --exp_tag {params.exp_tag}\
         --theta_prior 200000\
         --max_iter 200\
         --updated_sam\
         --outdir {params.outdir}\
         {input[0]}\
         {input[1]}\
         2>&1 | tee {log[0]}
        chmod 660 {output[0]}
        """

rule complete:
    input:
        expand("results/telescope/{samid_paired}_completed.txt", samid_paired=SAMPLES_PAIRED),
        expand("results/telescope/{samid_single}_completec.txt", samid_single=SAMPLES_SINGLE)
