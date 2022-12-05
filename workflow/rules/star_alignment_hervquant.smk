#! /usr/bin/env python
# -*- coding utf-8 -*-

# Add quant mode for gene TE_counts
# Make sure that star runs with stringtie
################################ STAR ALIGNMENT ################################

rule star_alignment:
    conda:
        "../envs/star.yaml"
    input:
        single = "samples/{samid}.fastq",
        genome = config['indexes']['star_hervquant']
    output:
        aligned_bam = "results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.Aligned.out.bam"
    params:
        out_prefix="results/star_alignment/{samid}/{samid}_STAR_hervquant_hg19.",
        outFilterMultimapNmax=10,
        outFilterMismatchNmax=7
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
