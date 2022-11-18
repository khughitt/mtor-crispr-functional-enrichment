#
# create batch adjusted version of RNA-Seq by creating a model with only batch and then extracting
# the model residuals.
#
# this is useful for visualization purposes, in order to estimate the impact of batch
# adjustment on the raw counts.
#
library(tidyverse)
library(limma)

set.seed(1)

raw_counts <- read_tsv(snakemake@input[[1]], col_types = cols()) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# load count table & sample metadata
sample_mdata <- read_tsv(snakemake@input[[2]], col_types = 'fff')

# batch model
batch <- sample_mdata$batch
batch_design <- model.matrix(~batch)

batch_voom <- voom(raw_counts, batch_design, normalize.method = 'none', plot = FALSE)

# create a linear model for batch and get residuals
batch_fit <- lmFit(batch_voom, design = batch_design)

limma_residuals <- residuals(batch_fit, batch_voom$E) %>%
  as.data.frame() %>%
  rownames_to_column("Geneid")

write_tsv(limma_residuals, snakemake@output[[1]])
write_tsv(sample_mdata, snakemake@output[[2]])
