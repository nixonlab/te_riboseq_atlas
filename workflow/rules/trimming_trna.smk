#! /usr/bin/env python
# -*- coding utf-8 -*-

################################ TRNA TRIMMING ################################

rule trna_trimming:
    conda:
        "../envs/seqkit.yaml"
    input:
        single = "samples/{samid}_noribo.fastq"
    output:
        ribo_fastq = temp("samples/{samid}_noribo_trimmed.fastq")
    threads: 20
    benchmark: "benchmarks/seqkit/{samid}_noribo_notrna.tsv"
    shell:
        '''
        seqkit seq {input} -m 25 -M 35 > {output}
        '''

############################## BOWTIE2 BUILD TRNA ##############################

rule bowtie_build:
    conda:
        "../envs/bowtie.yaml"
    input:
        fasta="custom_databases/tRNA_rRNA_hg19_ND.fa"
    output:
        directory=directory("databases/indexes/tRNArRNAref")
    shell:
        '''
        bowtie2-build {input.fasta} {output.directory}
        '''

############################### BOWTIE2 TRIMMING ###############################

rule bowtie_align:
    conda:
        "../envs/bowtie.yaml"
    input:
        fasta="samples/{samid}_noribo_trimmed.fastq",
        ref="databases/indexes/tRNArRNAref"
    output:
        fasta= temp("samples/{samid}_noribo_notrna_trimmed.fastq")
    threads: 8
    benchmark: "benchmarks/bowtie/{samid}_noribo_notrna_trimmed.tsv"
    shell:
        '''
        bowtie2 --phred33 --very-sensitive-local -p 8 --quiet -x {input.ref} -U {input.fasta} --un {output.fasta}
        '''
