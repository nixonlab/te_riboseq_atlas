"""
Download FASTQs
"""

rule fasterq_dump_paired:
    conda: "../envs/sra.yaml"
    output:
        temp("runs/{s_paired}/{s_paired}/{s_paired}.sra"),
        temp("runs/{s_paired}/{s_paired}_1.fastq"),
        temp("runs/{s_paired}/{s_paired}_2.fastq")
    params:
        tmpdir = config["fasterq_dump_tmp"],
        outdir = "runs/{s.paired}",
        sra_outdir = "runs/{s_paired}/{s_paired}"
    threads: 8
    resources:
        mem_mb = 10000, disk_mb = 60000
    log: "runs/{s_paired}/fasterq_sra_to_fastq.log"
    shell:
        """
        prefetch --ngc {params.dbgap_key} -O {params.outdir} {wildcards.s}
        fastq-dump --split-3 --skip-technical --origfmt -O {params.outdir} {params.sra_outdir} &> {log[0]}
        """

rule fasterq_dump_single:
    conda: "../envs/sra.yaml"
    output:
        temp("runs/{s_single}/{s_single}/{s_single}.sra"),
        temp("runs/{s_single}/{s_single}.fastq")
    params:
        tmpdir = config["fasterq_dump_tmp"],
        outdir = "runs/{s_single}",
        sra_outdir = "runs/{s_single}/{s_single}"
    threads: 8
    resources:
        mem_mb = 10000, disk_mb = 60000
    log: "runs/{s.single}/fasterq_sra_to_fastq.log"
    shell:
        """
        prefetch -O {params.outdir} {wildcards.s}
        fastq-dump --origfmt -O {params.outdir} {params.sra_outdir} &> {log[0]}
        """

rule cat_runids_to_samples_paired:
    input:
        R1 = lambda wc: expand("runs/{s_paired}/{s_paired}_1.fastq", s=SAMPLE_RUN_PAIRED[wc.samid_paired]), # input fastqs per run
        R2 = lambda wc: expand("runs/{s_paired}/{s_paired}_2.fastq", s=SAMPLE_RUN_PAIRED[wc.samid_paired])
    output:
        R1 = temp("samples/{samid_paired}_1.fastq"), # output fastqs per sample
        R2 = temp("samples/{samid_paired}_2.fastq")
    shell:
        """
        cat {input.R1} > {output.R1} # fastqs belonging to same sample (multipe runs) are combined
        cat {input.R2} > {output.R2} # fastqs belonging to same sample (multipe runs) are combined
        rm {input.R1} {input.R2}
        """

rule cat_runids_to_samples_single:
    input:
        single = lambda wc: expand("runs/{s_single}/{s_single}.fastq", s=SAMPLE_RUN_SINGLE[wc.samid_single]) # input fastqs per run
    output:
        single = temp("samples/{samid_single}.fastq") # output fastqs per sample
    shell:
        """
        cat {input.R1} > {output.R1} # fastqs belonging to same sample (multipe runs) are combined
        cat {input.R2} > {output.R2} # fastqs belonging to same sample (multipe runs) are combined
        rm {input.R1} {input.R2}
        """

rule conversion_complete:
    input:
        expand("samples/{samid_paired}_1.fastq", samid_paired=SAMPLES_PAIRED),
        expand("samples/{samid_single}.fastq", samid_single=SMAMPLES_SINGLE)
