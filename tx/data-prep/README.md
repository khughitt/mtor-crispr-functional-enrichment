# CRISPR-guided multi-omics: RNA-Seq transcript abundance estimation

## Overview

This directory contains a Snakemake pipeline used to perform basic RNA-Seq read pre-processing and
QA, and to generate transcript abundance estimates.

The pipeline was designed to be used with the raw reads (.fastq) from Nguyen et al. (2023).

Steps included:

1. QA ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Trimming ([cutadapt](https://cutadapt.readthedocs.io/en/stable/))
3. Abundance estimation ([HISAT2](https://daehwankimlab.github.io/hisat2/manual/))
4. SAM -> BAM ([samtools]())

## Usage

To run the pipeline, start by creating a conda environment with the necessary requirements:

```
conda create -n nguyen23-rnaseq-prep --file requirements.txt
```

Activate the environment using:

```
conda activate nguyen23-rnaseq-prep 
```

Next, copy the example config file located at "config/config.example.yml" to "config.yml", and
modify the directories to indicate where the pipeline should expect to find the required input
reads and HISAT2 index, and where outputs should be saved to:


```
cp config/config.example.yml config/config.yml
$EDITOR config/config.yml
```

Finally, launch the pipeline with a desired number of jobs:

```
snakemake -j <num_jobs>
```
