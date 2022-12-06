#! /usr/bin/env python
# -*- coding utf-8 -*-

################################## INDEX REFS ##################################

rule star_index_hervquant_hg19:
    input:
        fasta = 'databases/remotefiles/hervquant_hg19_reference.fa'
    output:
        directory("databases/indexes/STAR_hervquant_hg19"),
        "databases/indexes/STAR_hervquant_hg19/SAindex"
    conda: "../envs/star_hervquant.yaml"
    threads: snakemake.utils.available_cpu_count()
    shell:
        '''
STAR\
 --runThreadN {threads}\
 --runMode genomeGenerate\
 --genomeDir {output[0]}\
 --genomeSAindexNbases 7\
 --genomeChrBinNbits 8\
 --outFileNamePrefix {output[0]}\
 --genomeFastaFiles {input.fasta}\
 --limitGenomeGenerateRAM 54000000000
        '''
