"""
Download FASTQs
"""

rule fasterq_dump:
    conda: "../envs/sra.yaml"
    output:
        temp("runs/{s}/{s}/{s}.sra"),
        temp("runs/{s}/{s}.fastq")
    params:
        tmpdir = config["fasterq_dump_tmp"],
        outdir = "runs/{s}",
        sra_outdir = "runs/{s}/{s}"
    threads: 8
    resources:
        mem_mb = 10000, disk_mb = 60000
    log: "runs/{s}/fasterq_sra_to_fastq.log"
    shell:
        """
        prefetch -O {params.outdir} {wildcards.s}
        fastq-dump --origfmt -O {params.outdir} {params.sra_outdir} &> {log[0]}
        """

rule cat_runids_to_samples_single:
    input:
        single = lambda wc: expand("runs/{s}/{s}.fastq", s=SAMPLE_RUN[wc.samid]) # input fastqs per run
    output:
        single = temp("samples/{samid}.fastq") # output fastqs per sample
    shell:
        """
        cat {input.single} > {output.single}
        """

rule conversion_complete:
    input:
        expand("samples/{samid}.fastq", samid=SAMPLES)
