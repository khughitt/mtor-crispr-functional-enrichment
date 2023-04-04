# CRISPR-guided multi-omics: RNA-Seq Differential Expression & Functional Enrichment Analysis Pipeline

The pipeline in this directory accepts as input processed RNA-Seq reads from a previous step and
performs differential expression analysis using DESeq2 and functional enrichment analysis using
clusterProfiler.

To run the pipeline, start by creating a conda environment with the necessary requirements:

```
conda create -n nguyen22-rnaseq --file requirements.txt
```

Activate the environment using:

```
conda activate nguyen22-rnaseq
```

Next, copy the config file named "config/config.example.yml" to "config/config.yml", and modify the
filepaths to point to the location of the transcript abundance estimates from the previous step.

```
cp config/config.example.yml config/config.yml
$EDITOR config/config.yml
```

Finally, launch the pipeline with a desired number of jobs:

```
snakemake -j <num_jobs>
```
