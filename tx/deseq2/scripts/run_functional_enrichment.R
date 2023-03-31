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

# GSEA (all genes)
gene_scores <- setNames(deseq2_res$log2FoldChange, deseq2_res$symbol)

set.seed(cfg$random_seed)

gsea_res <- GSEA(gene_scores,
                 TERM2GENE = gene_sets,
                 seed = TRUE,
                 pvalueCutoff = cfg$min_padj,
                 nPermSimple = cfg$nperm)

gsea_res <- gsea_res@result %>%
  arrange(p.adjust)

write_tsv(gsea_res, snakemake@output[[1]])

# Gene set over-representation analysis (up- and down-regulated gene sets)
up_reg <- deseq2_res %>%
  filter(padj <= cfg$min_padj & log2FoldChange > 0) %>%
  pull(symbol)

down_reg <- deseq2_res %>%
  filter(padj <= cfg$min_padj & log2FoldChange < 0) %>%
  pull(symbol)

up_reg_gene_sets <- enricher(up_reg, TERM2GENE = gene_sets)@result %>%
  arrange(p.adjust)

down_reg_gene_sets <- enricher(down_reg, TERM2GENE = gene_sets)@result %>%
  arrange(p.adjust)

write_tsv(up_reg_gene_sets, snakemake@output[[2]])
write_tsv(down_reg_gene_sets, snakemake@output[[3]])

# create a combined xlsx file with results
tbls = list(
  "GSEA" = gsea_res,
  "over_rep_upreg" = up_reg_gene_sets,
  "over_rep_downreg" = down_reg_gene_sets
)

write.xlsx(tbls, file = snakemake@output[[4]])
