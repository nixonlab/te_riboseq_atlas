#! /usr/bin/env python
# -*- coding utf-8 -*-

################################ RIBO DETECTOR ################################

rule ribo_detector:
    conda:
        "../envs/ribo_detector.yaml"
    input:
        single = "samples/{samid}.fastq"
    output:
        ribo_fastq = temp("samples/{samid}_noribo.fastq")
    threads: 20
    benchmark: "benchmarks/ribo_detector/{samid}_ribo_detector.tsv"
    shell:
        '''
        ribodetector_cpu -t 20 \
        -l 30 \
	--chunk_size 256 \
        -i {input.single} \
        -o {output.ribo_fastq}
        '''
