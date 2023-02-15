#! /usr/bin/env python
# -*- coding utf-8 -*-

################################ STAR ALIGNMENT ################################

rule star_alignment_hervquant:
    conda:
        "../envs/star.yaml"
    input:
        single = "samples/{samid}_noribo_notrna_trimmed.fastq",
        genome = config['indexes']['star_hervquant']
    output:
        aligned_bam = temp("results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.Aligned.out.sam")
    params:
        out_prefix="results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.",
        outFilterMultimapNmax=3,
        outFilterMismatchNmax=1
    threads: config['star_alignment_threads']
    benchmark: "benchmarks/star_alignment/{samid}_star_alignment_hervquant.tsv"
    shell:
        '''
        STAR\
            --runThreadN {threads}\
            --outFileNamePrefix {params.out_prefix}\
            --outFilterMultimapNmax {params.outFilterMultimapNmax}\
            --outFilterMismatchNmax {params.outFilterMismatchNmax}\
            --genomeDir {input.genome}\
            --readFilesIn {input.single}
        '''
