# TE RiboSeq Atlas (Dopkins)

A transposable element atlas of ribo-seq data from multiple tissue types. (Project ID: PRJNA756018)

**To get DAG:** 

``` snakemake --profile profiles/aws  --forceall --dag | dot -Tpdf > dag.pdf   ```

**To get rule graph:** 

``` snakemake --profile profiles/aws  --forceall --rulegraph | dot -Tpdf > rulegraph.pdf   ```

**To get file graph:** 

``` snakemake --profile profiles/aws  --forceall --filegraph | dot -Tpdf > filegraph.pdf   ```

**To run pipeline:**

``` snakemake --profile profiles/aws/ all ```

**To run only Telescope in a sample-specific manner:**

``` snakemake --profile profiles/aws/ results/telescope/{sampleid1}_telescope_completed.txt results/telescope/{sampleid2}_telescope_completed.txt```

**To run only hervquant in a sample-specific manner:**

``` snakemake --profile profiles/aws/ results/salmon/{sampleid1}_salmon_completed.txt results/salmon/{sampleid2}_salmon_completed.txt ```

**Example to run both hervquant and telescope for 3 samples***:

```
snakemake --profile profiles/aws results/telescope/SAMN20849944_telescope_completed.txt results/salmon/SAMN20849944_salmon_completed.txt results/telescope/SAMN20849945_telescope_completed.txt results/salmon/SAMN20849945_salmon_completed.txt results/telescope/SAMN20849946_telescope_completed.txt results/salmon/SAMN20849946_salmon_completed.txt ```

**To modify pipeline:**

Change sample download table. 
