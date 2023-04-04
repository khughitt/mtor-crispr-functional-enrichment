# CRISPR-guided multiomics unveils direct mTOR inhibition of Hippo signaling through CDCP1: RNA-Seq / Phosphoproteomics Data Analysis

## Overview

This repo contains code for the RNA-Seq and phosphoproteomics analysis portions of the paper
"CRISPR-guided multiomics unveils direct mTOR inhibition of Hippo signaling through CDCP1" (Nguyen
_et al._, 2023)

The code has been organized into several task-specific directories:

1. RNA-Seq data preparation and transcript abundance estimation (`tx/data-prep`)
2. RNA-Seq differential expression & functional enrichment analysis (`tx/downstream`)
3. Phosphoproteomics functional enrichment analysis (`phospho/`)
4. PubMed mTOR/Hippo gene co-occurrence analysis (`pubmed/`)
5. Summary report (`summary/`)

The RNA-Seq analyses were performed using [Snakemake](https://snakemake.readthedocs.io/en/stable/)
pipelines.

The phosphoproteomics analysis and summary report are provided as RMarkdown files, while the PubMed
literature analysis is provided as a jupyter notebook.

In general, efforts have been made to use open-source libraries and to share all relevant analysis
code in this repo, in order to clearly document what has been done and to make it possible to
reproduce as much of the analyses as possible.

One part of the analysis that was performed separately and is not documented here, however, is the
phosphoproteomics data preparation and abundance estimation. This was performed using [Proteome
Discover](https://www.thermofisher.com/us/en/home/industrial/mass-spectrometry/liquid-chromatography-mass-spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/proteome-discoverer-software.html),
and the outputs from that analysis are what is included in the `data/phospho-xx.tsv` input datasets.
