# TE RiboSeq Atlas (Dopkins)

A transposable element atlas of ribo-seq data from multiple tissue types. (Project ID: PRJNA756018)

**To get DAG:** 

```
snakemake --profile profiles/aws  --forceall --dag | dot -Tpdf > dag.pdf
```

**To get rule graph:** 

```
snakemake --profile profiles/aws  --forceall --rulegraph | dot -Tpdf > rulegraph.pdf
```

**To get file graph:** 

```
snakemake --profile profiles/aws  --forceall --filegraph | dot -Tpdf > filegraph.pdf
```

**To run pipeline:**

```
snakemake --profile profiles/aws/ all
```

**To run only Telescope in a sample-specific manner:**

``` 
snakemake --profile profiles/aws/ results/telescope/{sampleid1}_telescope_completed.txt results/telescope/{sampleid2}_telescope_completed.txt
```

**To run only hervquant in a sample-specific manner:**

``` 
snakemake --profile profiles/aws/ results/salmon/{sampleid1}_salmon_completed.txt results/salmon/{sampleid2}_salmon_completed.txt
```

**Example to run both hervquant and telescope for 3 samples**:

```
snakemake --profile profiles/aws results/telescope/SAMN20849944_telescope_completed.txt results/salmon/SAMN20849944_salmon_completed.txt results/telescope/SAMN20849945_telescope_completed.txt results/salmon/SAMN20849945_salmon_completed.txt results/telescope/SAMN20849946_telescope_completed.txt results/salmon/SAMN20849946_salmon_completed.txt 
```

**Reference genomes used for telescope and hervquant:**

For a full list of reference genomes used, please refer to the config.yaml file. For STAR, we used [GRCh38.d1.vd1.fa.tar.gz](https://api.gdc.cancer.gov/data/254f697d-310d-4d7d-a27b-27fbf767a834). For Telescope, we used [retro.hg38.v1.gtf](https://github.com/mlbendall/telescope_annotation_db/raw/master/builds/retro.hg38.v1/transcripts.gtf). For hervquant, we used the reference from the original paper, which can be found [here](https://unclineberger.org/vincentlab/wp-content/uploads/sites/1083/2020/10/hervquant_hg19_reference.fa_.zip) and [here](https://unclineberger.org/vincentlab/wp-content/uploads/sites/1083/2020/10/hervquant-reference-file.zip).

**To modify pipeline:**

Change sample download table. 
