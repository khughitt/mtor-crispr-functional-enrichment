#
# functional enrichment analysis
#
library(tidyverse)
library(GSEABase)
library(clusterProfiler)
library(openxlsx)

cfg <- snakemake@config$functional_enrichment

# load MSigDB gene sets
gene_sets <- geneIds(getGmt(cfg$gmt))

gene_sets <- stack(gene_sets) %>%
  select(term = ind, gene = values)

# load deseq2 results
deseq2_res <- read_tsv(snakemake@input[[1]], col_types = cols())

# collapse multiple entries mapping to the same gene symbols;
# mostly impacts non-coding RNA genes (SNORDs, Y_RNA, etc.)
deseq2_res <- deseq2_res %>%
  group_by(symbol) %>%
  summarize(log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
            padj = mean(padj, na.rm = TRUE)) %>%
  arrange(-log2FoldChange) %>%
  ungroup()

# drop entry corresponding to genes which could not be mapped to symbols
deseq2_res <- deseq2_res %>%
  filter(!is.na(symbol))

# GSEA
gene_scores <- setNames(deseq2_res$log2FoldChange, deseq2_res$symbol)

set.seed(cfg$random_seed)

gsea_res <- GSEA(gene_scores,
                 TERM2GENE = gene_sets,
                 seed = TRUE,
                 pvalueCutoff = cfg$min_padj,
                 nPermSimple = cfg$nperm)

gsea_res_df <- gsea_res@result %>%
  arrange(p.adjust)

# save result as .tsv
write_tsv(gsea_res_df, snakemake@output[[1]])

# create a combined xlsx file with results
tbls = list("GSEA" = gsea_res_df)
write.xlsx(tbls, file = snakemake@output[[2]])

# generate plot of significant results
dotplot(gsea_res) +
  ggplot2::xlim(0, 0.8)

ggsave(snakemake@output[[3]], width = 1080, height = 800, units = "px", dpi = 192)
