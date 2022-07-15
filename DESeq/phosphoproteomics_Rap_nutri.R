#read in and preview the raw counts file
setwd("~/Desktop/Desktop Files/Phosphoproteomics")
mtor_raw_counts <- read.csv('final_phosphoproteomics_caltech_nutri.csv' , row.names = 1)
head(mtor_raw_counts)
str(mtor_raw_counts)
#mtor_raw_counts$symbol <- str_extract(mtor_raw_counts$Phospho_Site, "[^ ]+")
#mtor_raw_counts <- mtor_raw_counts %>% relocate(symbol, .before = Phospho_Site)
#mtor_raw_counts <- mtor_raw_counts %>% relocate(Phospho_Site, .after = Uniprot_ID)

mtor_raw_counts <- dplyr::select(mtor_raw_counts, Phospho_Site, control_nutri_1, control_nutri_2, control_nutri_3, Rap_nutri_1, Rap_nutri_2, Rap_nutri_3, Rap_nutri_4)
row.names(mtor_raw_counts) <- mtor_raw_counts$Phospho_Site
mtor_raw_counts$Phospho_Site <- NULL #erase column but leave the row names
#write.csv(mtor_raw_counts, row.names = TRUE, "final_phosphoproteomics_caltech_nutri_phospho_names_included.csv")

#mtor_raw_counts <- read.csv('final_phosphoproteomics_caltech_nutri_phospho_names_included.csv', header = T, row.names = 1)
#full_nutri_global_protein <- read_excel('Nutrient group global 2020_0214ProteomeFractionatedFullProteinList.xlsx')
#full_nutri_global_protein <- as.data.frame(full_nutri_global_protein)
#intersection <- intersect(full_nutri_global_protein$`Gene Symbol`, mtor_raw_counts$symbol) #find the intersection of the global proteomics list and the phosphoproteomics list

# TAILOR RAW COUNTS

library(tibble)  # for `rownames_to_column` and `column_to_rownames`
#mtor_raw_counts<- as.matrix(as.data.frame(mtor_raw_counts))
Rap_nutri <- dplyr::select(mtor_raw_counts, control_nutri_1, control_nutri_2, control_nutri_3, Rap_nutri_1, Rap_nutri_2, Rap_nutri_3, Rap_nutri_4)

#Rap_nutri <- arrange(Rap_nutri, Phospho_Site)
head(Rap_nutri)
Rap_nutri<- as.data.frame(Rap_nutri)
#rownames_to_column(Rap_nutri$Phospho_Site)


# SET UP METADATA
#create the Rictors for the metadata
cell_type <- c("Vector","Vector","Vector", "Raptor","Raptor","Raptor","Raptor")
starv_status <- c("Nutri","Nutri","Nutri", "Nutri","Nutri","Nutri","Nutri")
all_comb <- c("control_nutri","control_nutri","control_nutri", "Rap_Nutri","Rap_Nutri","Rap_Nutri","Rap_Nutri")

#combine into a dataframe for the metadata
Rap_v_Control_nutri_metadata <- data.frame(cell_type, starv_status, all_comb)
#specify row names to match up to count matrix column names and preview the metadata
rownames(Rap_v_Control_nutri_metadata) <- c("control_nutri_1","control_nutri_2","control_nutri_3", "Rap_nutri_1","Rap_nutri_2","Rap_nutri_3","Rap_nutri_4")
head(Rap_v_Control_nutri_metadata)


# CHECK IF METADATA CORRESPONDS WITH RAW COUNTS
#call row and column names of the metadata and raw counts, respectively
rownames(Rap_v_Control_nutri_metadata)
colnames(Rap_nutri)

#check if sample names match (between the column names of the data matrix and row names of the metadata)
all(rownames(Rap_v_Control_nutri_metadata)==colnames(Rap_nutri))
#returns TRUE, so we do not have to reorder and can move on to making the DESeq2 object

# CHECK IF METADATA CORRESPONDS WITH RAW COUNTS

#View both the count matrix and metadata
#View(Rap_nutri)
#View(Rap_v_Control_nutri_metadata)

#call row and column names of the metadata and raw counts, respectively
rownames(Rap_v_Control_nutri_metadata)
colnames(Rap_v_Control_nutri_metadata)

#check if sample names match (between the column names of the data matrix and row names of the metadata)
all(rownames(Rap_v_Control_nutri_metadata)==colnames(Rap_nutri))
#returns TRUE, so we do not have to reorder and can move on to making the DESeq2 object

#load libraries
library(DESeq2)
library(tidyverse)

# BEGIN DESEQ2 ANALYSIS
#create the DESeq2 object
dds_Rapnutri <- DESeqDataSetFromMatrix(countData = round(Rap_nutri),
                                       colData = Rap_v_Control_nutri_metadata,
                                       design = ~ cell_type)


# NORMALIZED COUNTS

#calculate normalized counts using size factors
dds_Rapnutri <- estimateSizeFactors(dds_Rapnutri)
sizeFactors(dds_Rapnutri)

#extract normalized counts from DESeq2 object
normalized_Rapnutri_counts <- counts(dds_Rapnutri, normalized = TRUE)
View(normalized_Rapnutri_counts)

#log transform the normalized counts to improve visualization of clustering in visuals
vsd_Rapnutri <- vst(dds_Rapnutri, blind = TRUE)

#extract vst transformed normalized counts as a matrix from the vsd object
vsd_mat_Rapnutri <- assay(vsd_Rapnutri)
#compute pairwise correlation values between each pair of samples
vsd_cor_Rapnutri <- cor(vsd_mat_Rapnutri)

# CORRELATION HEATMAP, PCA PLOT, AND OTHER VISUALS

#load pheatmap library
library(pheatmap)

#plot the correlation heatmap
pheatmap(vsd_cor_Rapnutri, annotation = dplyr::select(Rap_v_Control_nutri_metadata, all_comb))
#pheatmap(vsd_cor_Rapnutri)

#load tools for PCA
#install.packages("ggfortify")
library(ggfortify)

#plot PCA
#autoplot(prcomp(vsd_cor_Rapnutri), data=Rapnutri_metadata, colour="all_comb", label=TRUE, labl.size=3)
pc<-autoplot(prcomp(vsd_cor_Rapnutri), data=Rap_v_Control_nutri_metadata, colour="all_comb", label=FALSE, size=4.00)
pc+scale_color_manual(values=c("deepskyblue1", "green1"))+theme_classic()
ggsave("RapnutriPCA_highres.png", width=7, height=7, dpi=300)

#MODEL FITTING & VISUALIZATION

#perform model fitting with DESeq
dds_Rapnutri <- DESeq(dds_Rapnutri)

#calculate mean and variance for each gene of the samples
Rapnutri_mean_counts <- apply(Rap_nutri, 1, mean)
Rapnutri_variance_counts <- apply(Rap_nutri, 1, var)
#create a data frame to plot the relationship between mean and variance for each gene
df <- data.frame(Rapnutri_mean_counts, Rapnutri_variance_counts)

#plot the mean and variance for each gene
ggplot(df)+
  geom_point(aes(x=Rapnutri_mean_counts, y=Rapnutri_variance_counts))+
  scale_y_log10()+
  scale_x_log10()+
  xlab("Mean counts per gene")+
  ylab("Variance per gene")

#plot dispersion estimates
#plotDispEsts(dds_Rapnutri)

# DE ANALYSIS RESULTS
contrast1<-c("cell_type","Raptor", "Vector")

#extract the results of the differential expression analysis
Rapnutri_res <- results(dds_Rapnutri, contrast=contrast1, alpha = 0.05)


# SHRINK LOG2FC

#shrink log2 fold changes
Rapnutri_res <- lfcShrink(dds_Rapnutri, contrast=contrast1, type="ashr", res=Rapnutri_res)
#plot MA plot again 
#plotMA(Rapnutri_res, ylim = c(-2,2))

#descriptions for columns in results table
#mcols(Rapnutri_res)
#view the first few rows of the results table
#head(Rapnutri_res)
#preview DEGs and genes filtered
summary(Rapnutri_res)

# SIGNIFICANT GENES

#test for significant genes
Rapnutri_res <- results(dds_Rapnutri, alpha = 0.05, contrast=contrast1, lfcThreshold = 0) ####### could just leave the default value. 
#reshrink fold changes with modified results
Rapnutri_res <- lfcShrink(dds_Rapnutri, contrast=contrast1, type="ashr", res=Rapnutri_res)

#run summary again
summary(Rapnutri_res)




# ANNOTABLES FOR GENE SPECIFICATIONS             

#to better understand which genes the results pertain to, use annotables
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
#load libraries
library(dplyr)
library(annotables)
grch38

#extract gene name from the phosphopeptide name
#to annotate genes, first turn results table into data frame
library(stringr)
Rapnutri_res_all <- data.frame(Rapnutri_res)
Rapnutri_res_all$symbol <- str_extract(rownames(mtor_raw_counts), "[^ ]+") #extract the symbol from the rowname

library(dplyr)
require(dplyr)
Rapnutri_res_all <- Rapnutri_res_all %>% relocate(symbol, .before = baseMean)

#SIGNIFICANT GENES
Rapnutri_res_all<-filter(Rapnutri_res_all, !is.na(pvalue))
Rapnutri_res_all<-filter(Rapnutri_res_all, !is.na(padj))

#extract significant DE genes
Rapnutri_res_sig <- subset(Rapnutri_res_all, padj <= 0.05)

#order genes by padj to generate final table of significant results
Rapnutri_res_sig <- arrange(Rapnutri_res_sig, padj)
#View(Rapnutri_res_sig)

#subset normalized counts to only significant DE genes
normalized_Rapnutri_counts <- as.data.frame(normalized_Rapnutri_counts) #have the phospho site as row names
symbol <- str_extract(rownames(mtor_raw_counts), "[^ ]+") #extract the symbol from the rowname
normalized_Rapnutri_counts$symbol <- symbol
normalized_Rapnutri_counts <- normalized_Rapnutri_counts %>% relocate(symbol, .before = control_nutri_1)

sig_norm_counts_Rapnutri <- normalized_Rapnutri_counts[rownames(Rapnutri_res_sig),] #extract just the significant phosphopeptides from the normalized data frame
sig_norm_counts_Rapnutri$Phospho_site <- rownames(sig_norm_counts_Rapnutri)
#sig_norm_counts_Rapnutri<- left_join(x = sig_norm_counts_Rapnutri,y = grch38[,c("ensgene","symbol","description")], by = "symbol")
#sig_norm_counts_Rapnutri <- filter(sig_norm_counts_Rapnutri, !is.na(ensgene)) #filter out if phosphopeptide corresponding gene doesn't have ensemble gene IDs

ACTsignormcounts<- left_join(x = sig_norm_counts_Rapnutri,y = grch38[,c("ensgene","symbol","description")], by = "symbol")
#ACTsignormcounts <- filter(ACTsignormcounts, !is.na(ensgene)) #filter out if phosphopeptide corresponding gene doesn't have ensemble gene IDs

#ACTsignormcounts$symbol <- NULL
sig_norm_counts_Rapnutri$symbol <- NULL
sig_norm_counts_Rapnutri$Phospho_site <- NULL

#HEATMAP
#choose a color palette from RColorBrewer
library(RColorBrewer)
heat_colors <- brewer.pal(6, "Blues")
my_colors <- list(all_comb=c(Rap_Nutri="green1", control_nutri="deepskyblue1"))

#plot the heatmap
#rownames(sig_norm_counts_Rapnutri) = colnames(mat)

expheat<-pheatmap(mat=sig_norm_counts_Rapnutri[1:20, ],
                  color = heat_colors,
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  labels_row = ACTsignormcounts[1:20, ]$Phospho_site,
                  show_rownames = TRUE,
                  annotation = dplyr::select(Rap_v_Control_nutri_metadata, all_comb),
                  annotation_colors = my_colors,
                  scale = "row")
#dev.off()
ggsave("Rapnutri_EXPHEAT50_highres.png", expheat, width=8, height=8, dpi=300)



#TOP 20 DEGs

#gather top 20 significant genes from normalized counts
top_20_Rapnutri <- data.frame(sig_norm_counts_Rapnutri)[1:20, ]
top_20_Rapnutri <- rownames_to_column(top_20_Rapnutri, var = "Phospho_Site")
top_20_Rapnutri <- gather(top_20_Rapnutri, key = "sample_name", value = "normalized_counts", 2:8)

#merge the metadata to color the plot by sample group
Rapnutri_meta_new <- rownames_to_column(Rap_v_Control_nutri_metadata, var = "sample_name")
top_20_Rapnutri <- inner_join(top_20_Rapnutri, Rapnutri_meta_new, by = "sample_name")

#create the expression plot
m=ggplot(top_20_Rapnutri)+
  geom_point(aes(x=reorder(Phospho_Site,-normalized_counts), y=normalized_counts, colour= all_comb))+
  scale_y_log10()+
  xlab("Phosphorylation Site")+
  ylab("Normalized Counts")+
  ggtitle("Top 20 Significant DE Phosphopeptides")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.05))
m+scale_color_manual(values=c("deepskyblue1", "green1"))

ggsave("rapnutriDEgen4.png", width=8.5, height=7, dpi=300)
#write.csv(as.data.frame(Rapnutri_res), "/Users/holzergc/Desktop/Rapnutri_res.csv")


#UP AND DOWNREGULATED GENES 
#downregulated phosphopeptides
rapnutri_downreg <- filter(Rapnutri_res_all, log2FoldChange <= -1.0, padj <= 0.05)
rapnutri_downreg <- arrange(rapnutri_downreg, padj)

#export significantly downreg genes
write.csv(row.names(rapnutri_downreg), "phospho_Rap_nutri_sig_downreg.csv")

unique_downreg_Rap_nutri <- unique(rapnutri_downreg$symbol)
write.csv(unique_downreg_Rap_nutri, "phospho_Rap_nutri_sig_downreg_symbol.csv")


#normalized_Rapnutri_counts <- rownames_to_column(normalized_Rapnutri_counts, var = "Phospho_Site") 
rapnutri_downreg <- rownames_to_column(rapnutri_downreg, var = "Phospho_Site") 

sig_norm_counts_downreg <- normalized_Rapnutri_counts[rapnutri_downreg$Phospho_Site, ]

sig_norm_counts_downreg <- as.data.frame(sig_norm_counts_downreg)
sig_norm_counts_downreg <- rownames_to_column(sig_norm_counts_downreg, var = "Phospho_Site")
sig_norm_counts_downreg <- left_join(x = sig_norm_counts_downreg,y = grch38[,c("ensgene","symbol","description")], by = "symbol")
sig_norm_counts_downreg <- sig_norm_counts_downreg %>% mutate_all(na_if,"")

#sig_norm_counts_downreg <- filter(sig_norm_counts_downreg, !is.na(ensgene)) #filter out if phosphopeptide corresponding gene doesn't have ensemble gene IDs


#sig_norm_counts_downreg <- filter(sig_norm_counts_downreg, !is.na(description))
#sig_norm_counts_downreg <- dplyr::select(sig_norm_counts_downreg, -(symbol:description))

#gather top 20 significant genes from normalized counts
top_20_downreg <- data.frame(sig_norm_counts_downreg)[1:20, ]
top_20_downreg <- data.frame(sig_norm_counts_downreg)[1:20, ]

write.csv(top_20_downreg, row.names = TRUE, "Rap_nutri_top_20_downreg.csv")

#merge the metadata to color the plot by sample group
top_20_downreg <- gather(top_20_downreg, key = "sample_name", value = "normalized_counts", 3:9)
Rapnutri_meta_new <- rownames_to_column(Rap_v_Control_nutri_metadata, var = "sample_name")
top_20_downreg <- inner_join(top_20_downreg, Rapnutri_meta_new, by = "sample_name")

#top_20_downreg <- filter(top_20_downreg, !is.na(Phospho_site))

#top_20_downreg <- data.frame(top_20_downreg)[1:20, ]


#create the expression plot
d=ggplot(top_20_downreg)+
  geom_point(aes(x=reorder(Phospho_Site,normalized_counts), y=normalized_counts, colour= all_comb))+
  scale_y_log10()+
  xlab("Phosphopeptides")+
  ylab("Normalized Counts")+
  ggtitle("Top 20 Significantly Downregulated Phosphopeptides")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.05))
d+scale_color_manual(values=c("deepskyblue1", "green1"))

ggsave("rapnutriDEdownreg4.png", width=8.5, height=7, dpi=300)


#UPregulated phosphopeptides
rapnutri_upreg <- filter(Rapnutri_res_all, log2FoldChange >= 1.00, padj <= 0.05)
rapnutri_upreg <- arrange(rapnutri_upreg, padj)

#export significantly upreg phosphopeptides
write.csv(row.names(rapnutri_upreg), "phospho_Rap_nutri_sig_upreg.csv")

#normalized_Rapnutri_counts <- rownames_to_column(normalized_Rapnutri_counts, var = "Phospho_Site") 
rapnutri_upreg <- rownames_to_column(rapnutri_upreg, var = "Phospho_Site") 

sig_norm_counts_upreg <- normalized_Rapnutri_counts[rapnutri_upreg$Phospho_Site, ]

sig_norm_counts_upreg <- as.data.frame(sig_norm_counts_upreg)
sig_norm_counts_upreg <- rownames_to_column(sig_norm_counts_upreg, var = "Phospho_Site")
sig_norm_counts_upreg <- left_join(x = sig_norm_counts_upreg,y = grch38[,c("ensgene","symbol","description")], by = "symbol")
sig_norm_counts_upreg <- sig_norm_counts_upreg %>% mutate_all(na_if,"")

#sig_norm_counts_upreg <- filter(sig_norm_counts_upreg, !is.na(ensgene)) #filter out if phosphopeptide corresponding gene doesn't have ensemble gene IDs

#sig_norm_counts_downreg <- filter(sig_norm_counts_downreg, !is.na(description))
#sig_norm_counts_downreg <- dplyr::select(sig_norm_counts_downreg, -(symbol:description))

#gather top 20 significant genes from normalized counts
top_20_upreg <- data.frame(sig_norm_counts_upreg)[1:20, ]
#top_20_upreg <- data.frame(sig_norm_counts_upreg)[1:7, ]

#save CSV file with top 20 upreg
write.csv(top_20_upreg, row.names = TRUE, "Rap_nutri_top_20_upreg.csv")

#merge the metadata to color the plot by sample group
top_20_upreg <- gather(top_20_upreg, key = "sample_name", value = "normalized_counts", 3:9)
Rapnutri_meta_new <- rownames_to_column(Rap_v_Control_nutri_metadata, var = "sample_name")
top_20_upreg <- inner_join(top_20_upreg, Rapnutri_meta_new, by = "sample_name")

#create the expression plot
d=ggplot(top_20_upreg)+
  geom_point(aes(x=reorder(Phospho_Site,normalized_counts), y=normalized_counts, colour= all_comb))+
  scale_y_log10()+
  xlab("Phosphopeptides")+
  ylab("Normalized Counts")+
  ggtitle("Top 20 Significant Upregulated Phosphopeptides")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.05))
d+scale_color_manual(values=c("deepskyblue1", "green1"))

ggsave("rapnutriDEupreg4.png", width=8.5, height=7, dpi=300)





#LABS FOR VOLCANOS
#upreg
volc_upreg <- data.frame(sig_norm_counts_upreg)[1:10, ]
volc_upreg <- gather(volc_upreg, key = "sample_name", value = "normalized_counts", 2:9)

volc_upreg <- inner_join(volc_upreg, Rapnutri_meta_new, by = "sample_name")
volc_upreg <- volc_upreg %>% mutate_all(na_if,"")
volc_upreg <- filter(volc_upreg, !is.na(description))

#downreg
volc_downreg <- data.frame(sig_norm_counts_downreg)[1:10, ]
volc_downreg <- gather(volc_downreg, key = "sample_name", value = "normalized_counts", 2:9)

volc_downreg <- inner_join(volc_downreg, Rapnutri_meta_new, by = "sample_name")




#VOLCANO PLOT
library(EnhancedVolcano)
allDE<- rbind(top_20_downreg,top_20_upreg)

Rapnutri_res_all <- rownames_to_column(Rapnutri_res_all, var = "Phospho_Site")
Rapnutri_res_all <- arrange(Rapnutri_res_all, padj)


EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$Phospho_Site,
                selectLab = Rapnutri_res_all$Phospho_Site[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #xlim = c(-20, 20),
                #ylim = c(0, 275),
                pointSize = 1.25,
                labSize = 1.5,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                shape=15,
                arrowheads = FALSE,
                col=c("black","black","black","green1"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri',
                legendPosition = "none")
ggsave("RapnutriEVolc_highres1.png", width=7, height=7, dpi=300)

EnhancedVolcano(Rapnutri_res_all,
                lab = NA,
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #xlim = c(-20, 20),
                #ylim = c(0, 275),
                shape=15,
                pointSize = 1.25,
                labSize = 3,
                col=c("black","black","black","green1"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri',
                legendPosition = "none")
ggsave("RapnutriEVolcnolab_highres1.png", width=7, height=7, dpi=300)


EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$Phospho_Site,
                selectLab = 'CDCP1 Ser797',
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #xlim = c(-20, 20),
                #ylim = c(0, 275),
                shape=15,
                pointSize = 1.25,
                labSize = 1.5,
                col=c("black","black","black","green1"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri',
                legendPosition = "none")

ggsave("RapnutriEVolcnolab_CDCP1_FCcutoff_1.0_p_val_0.05.png", width=7, height=7, dpi=300)

#to determine which of the phosphopeptides in the mtor pathway are statistically significant (FC>1.5 or <1.5, and p-val <= 0.01)
rapnutri_all_sig <- rbind(rapnutri_downreg, rapnutri_upreg)
mTOR_intermediates <- c('AKT1', 'AKT1S1', 'AKT2', 'AKT3', 'ATP6V1A', 'ATP6V1B1', 'ATP6V1B2', 'ATP6V1C1', 'ATP6V1C2', 'ATP6V1D', 'ATP6V1E1', 'ATP6V1E2', 'ATP6V1F', 'ATP6V1G1', 'ATP6V1G2', 'ATP6V1G3', 'ATP6V1H', 'BRAF', 'CAB39', 'CAB39L', 'CHUK', 'CLIP1', 'DDIT4', 'DEPDC5', 'DEPTOR', 'DVL1', 'DVL2', 'DVL3', 'EIF4B', 'EIF4E', 'EIF4E1B', 'EIF4E2', 'EIF4EBP1', 'FLCN', 'FNIP1', 'FNIP2', 'FZD1', 'FZD10', 'FZD2', 'FZD3', 'FZD4', 'FZD5', 'FZD6', 'FZD7', 'FZD8', 'FZD9', 'GRB10', 'GRB2', 'GSK3B', 'HRAS', 'IGF1', 'IGF1R', 'IKBKB', 'INS', 'INSR', 'IRS1', 'KRAS', 'LAMTOR1', 'LAMTOR2', 'LAMTOR3', 'LAMTOR4', 'LAMTOR5', 'LPIN1', 'LRP5', 'LRP6', 'MAP2K1', 'MAP2K2', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MIOS', 'MLST8', 'MTOR', 'NPRL2', 'NPRL3', 'NRAS', 'PDPK1', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PRKAA1', 'PRKAA2', 'PRKCA', 'PRKCB', 'PRKCG', 'PRR5', 'PTEN', 'RAF1', 'RHEB', 'RHOA', 'RICTOR', 'RNF152', 'RPS6', 'RPS6KA1', 'RPS6KA2', 'RPS6KA3', 'RPS6KA6', 'RPS6KB1', 'RPS6KB2', 'RPTOR', 'RRAGA', 'RRAGB', 'RRAGC', 'RRAGD', 'SEC13', 'SEH1L', 'SESN2', 'SGK1', 'SKP2', 'SLC38A9', 'SLC3A2', 'SLC7A5', 'SOS1', 'SOS2', 'STK11', 'STRADA', 'STRADB', 'TELO2', 'TNF', 'TNFRSF1A', 'TSC1', 'TSC2', 'TTI1', 'ULK1', 'ULK2', 'WDR24', 'WDR59', 'WNT1', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16', 'WNT2', 'WNT2B', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B')
HIPPO_intermediates <- c('ACTB', 'ACTG1', 'AFP', 'AJUBA', 'AMH', 'AMOT', 'APC', 'APC2', 'AREG', 'AXIN1', 'AXIN2', 'BBC3', 'BIRC2', 'BIRC5', 'BMP2', 'BMP4', 'BMP5', 'BMP6', 'BMP7', 'BMP8A', 'BMP8B', 'BMPR1A', 'BMPR2', 'BTRC', 'CCND1', 'CCND3', 'CDH1', 'CRB1', 'CRB2', 'CSNK1D', 'CTGF', 'CTNNA1', 'CTNNA2', 'CTNNA3', 'CTNNB1', 'DLG1', 'DLG2', 'DLG3', 'DLG4', 'DVL1', 'DVL3', 'FBXW11', 'FGF1', 'FRMD1', 'FRMD6', 'FZD1', 'FZD10', 'FZD2', 'FZD3', 'FZD4', 'FZD5', 'FZD6', 'FZD7', 'FZD8', 'FZD9', 'GDF5', 'GDF6', 'GDF7', 'GLI2', 'GSK3B', 'ID1', 'ID2', 'ITGB2', 'LATS1', 'LATS2', 'LEF1', 'LIMD1', 'LLGL1', 'LLGL2', 'MOB1A', 'MOB1B', 'MPP5', 'MYC', 'NF2', 'NKD1', 'PARD3', 'PARD6A', 'PARD6B', 'PARD6G', 'PPP1CA', 'PPP1CB', 'PPP1CC', 'PPP2CA', 'PPP2CB', 'PPP2R1A', 'PPP2R1B', 'PPP2R2A', 'PPP2R2B', 'PPP2R2C', 'PPP2R2D', 'PRKCI', 'PRKCZ', 'RASSF1', 'RASSF6', 'SAV1', 'SCRIB', 'SERPINE1', 'SMAD1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD7', 'SNAI2', 'SOX2', 'STK3', 'TCF7', 'TCF7L1', 'TCF7L2', 'TEAD1', 'TEAD2', 'TEAD3', 'TEAD4', 'TGFB1', 'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TP53BP2', 'TP73', 'WNT1', 'WNT10A', 'WNT10B', 'WNT11', 'WNT16', 'WNT2', 'WNT2B', 'WNT3', 'WNT3A', 'WNT4', 'WNT5A', 'WNT5B', 'WNT6', 'WNT7A', 'WNT7B', 'WNT8A', 'WNT8B', 'WNT9A', 'WNT9B', 'WTIP', 'WWC1', 'WWTR1', 'YAP1', 'YWHAB', 'YWHAE', 'YWHAG', 'YWHAH', 'YWHAQ', 'YWHAZ')

intersection_sig_v_mtor <- intersect(mTOR_intermediates, rapnutri_all_sig$symbol)
intersection_sig_v_hippo <- intersect(HIPPO_intermediates, rapnutri_all_sig$symbol)

write.table(intersection_sig_v_mtor, file="phospho_Rap_nutri_v_control_mTOR_intermediates.csv")
write.table(intersection_sig_v_hippo, file="phospho_Rap_nutri_v_control_HIPPO_intermediates.csv")


#MTOR intermediates
EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$Phospho_Site,
                selectLab = c('EIF4EBP1 Ser65', 'EIF4B Ser422'),
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #xlim = c(-20, 20),
                #ylim = c(0, 270),
                shape=15,
                pointSize = 1.25,
                labSize = 0.9,
                labCol = "red",
                col=c("black","black","black","green1"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, mTOR intermediates',
                legendPosition = "none")
ggsave("Rapnutri_Volcnolab_FCcutoff_1_p_val_0.05_mTOR_intermediates.png", width=7, height=7, dpi=300)

#HIPPO intermediates
EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$Phospho_Site,
                selectLab = intersection_sig_v_hippo,
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #xlim = c(-20, 20),
                #ylim = c(0, 270),
                shape=15,
                pointSize = 1.25,
                labSize = 0.9,
                labCol = "red",
                col=c("black","black","black","green1"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, HIPPO intermediates',
                legendPosition = "none")
ggsave("Rapnutri_Volcnolab_FCcutoff_1_p_val_0.05_HIPPO_intermediates.png", width=7, height=7, dpi=300)


#############################################################################################
#Figure 2 volcano plots (labels of mTOR- and HIPPO-associated genes)
library(data.table) # Not necessary but super-useful
library(ggrepel)    # This is to avoid labels overalapping each other
library(ggplot2)

#import data
mtorc1_tx_down_reg_mtor_hippo_rankings <- read.csv('~/Desktop/Desktop Files/Phosphoproteomics/Keith_pubtator_analysis/mtorc1_tx_down-reg_mtor_hippo_rankings.csv', row.names = 1)
mtorc1_tx_up_reg_mtor_hippo_rankings <- read.csv('~/Desktop/Desktop Files/Phosphoproteomics/Keith_pubtator_analysis/mtorc1_tx_up-reg_mtor_hippo_rankings.csv', row.names = 1)

mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor <- arrange(mtorc1_tx_down_reg_mtor_hippo_rankings, mtor) #sort the mTOR rankings
mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor <- select(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor, symbol, mtor) #select just the desired columns
mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor <- filter(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor, !is.na(mtor)) #get rid of the rows that don't have a ranking
top_10_downreg_mtor <- data.frame(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor)[1:10, ] #select the top 10-ranked genes

mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO <- arrange(mtorc1_tx_down_reg_mtor_hippo_rankings, hippo) #sort the HIPPO rankings
mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO <- select(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO, symbol, hippo) #see above for method
mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO<-filter(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO, !is.na(hippo))
top_10_downreg_HIPPO <- data.frame(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO)[1:10, ]

#repeat with upreg genes
mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor <- arrange(mtorc1_tx_up_reg_mtor_hippo_rankings, mtor) #sort the mTOR rankings
mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor <- select(mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor, symbol, mtor) #select just the desired columns
mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor <- filter(mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor, !is.na(mtor)) #get rid of the rows that don't have a ranking
top_10_up_mtor <- data.frame(mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor)[1:10, ] #select the top 10-ranked genes

mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO <- arrange(mtorc1_tx_up_reg_mtor_hippo_rankings, hippo) #sort the HIPPO rankings
mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO <- select(mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO, symbol, hippo) #see above for method
mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO<-filter(mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO, !is.na(hippo))
top_10_upreg_HIPPO <- data.frame(mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO)[1:10, ]

#marge the upreg and downreg lists
mtorc1_tx_merged_mtor_genes <- rbind(top_10_downreg_mtor, top_10_up_mtor)
mtorc1_tx_merged_mtor_genes_2 <- arrange(mtorc1_tx_merged_mtor_genes, mtor) #sort the mTOR rankings

mtorc1_tx_merged_HIPPO_genes <- rbind(top_10_downreg_HIPPO, top_10_upreg_HIPPO)
mtorc1_tx_merged_HIPPO_genes_2 <- arrange(mtorc1_tx_merged_HIPPO_genes, hippo) #sort the HIPPO rankings

#alternative (merge all the upreg and downreg genes, rather than just the top 10)
mtorc1_tx_merged_mtor_genes_all <- rbind(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_mtor, mtorc1_tx_up_reg_mtor_hippo_rankings_sort_mtor)
mtorc1_tx_merged_mtor_genes_all_2 <- arrange(mtorc1_tx_merged_mtor_genes_all, mtor) #sort the mTOR rankings

mtorc1_tx_merged_HIPPO_genes_all <- rbind(mtorc1_tx_down_reg_mtor_hippo_rankings_sort_HIPPO, mtorc1_tx_up_reg_mtor_hippo_rankings_sort_HIPPO)
mtorc1_tx_merged_HIPPO_genes_all_2 <- arrange(mtorc1_tx_merged_HIPPO_genes_all, hippo) #sort the HIPPO rankings

#append the log2FC and pvalues
mtorc1_tx_merged_mtor_genes_all_3 <- left_join(x = mtorc1_tx_merged_mtor_genes_all_2, y = Rapnutri_res_all, by = "symbol")
mtorc1_tx_merged_HIPPO_genes_all_3 <- left_join(x = mtorc1_tx_merged_HIPPO_genes_all_2, y = Rapnutri_res_all, by = "symbol")


# keyvals <- ifelse(
#   Rapnutri_res_all$symbol %in% mtorc1_tx_merged_mtor_genes_2$symbol, 'purple',
#   ifelse(abs(Rapnutri_res_all$log2FoldChange) >= 2 & Rapnutri_res_all$padj <= 0.05, 'grey',
#          'black')

# keyvals[is.na(keyvals)] <- 'black'
# names(keyvals)[keyvals == 'purple'] <- 'mTOR_gene'
# names(keyvals)[keyvals == 'grey'] <- 'sig_gene'
# names(keyvals)[keyvals == 'black'] <- 'gene'

#sort the genes that I want to label
#colnames(keyvals) <- c("gene_type, ")
#keyvals_sorted <- Rapnutri_res_all[order(unlist(keyvals), decreasing = TRUE),] #sort by the type of gene

## mTOR ranked gene list
EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$symbol,
                selectLab = mtorc1_tx_merged_mtor_genes_2$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                #colCustom = keyvals,
                xlim = c(-3, 3),
                ylim = c(0, 60),
                colAlpha = 1,
                pointSize = 1,
                labSize = 1.5,
                labCol = "red",
                col=c("black","black","black","grey"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, custom color over-ride ',
                legendPosition = "none")
ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_ONLY_mTOR_associated_genes.png", width=7, height=7, dpi=300)

#just plot the mTOR-associated genes to overlay with plot above
EnhancedVolcano(mtorc1_tx_merged_mtor_genes_all_3,
                lab = mtorc1_tx_merged_mtor_genes_all_3$symbol,
                selectLab = mtorc1_tx_merged_mtor_genes$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.6,
                labCol = "red",
                xlim = c(-3, 3),
                ylim = c(0, 60),
                col=c("black","black","black","purple"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, mTOR-associated genes ',
                legendPosition = "none",
                #drawConnectors = TRUE,
                widthConnectors = 0.5,
                typeConnectors = "closed", endsConnectors = "first", lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10")
ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_mTOR_associated_genes_small_labels.png", width=7, height=7, dpi=300)

EnhancedVolcano(mtorc1_tx_merged_mtor_genes_all_3,
                lab = mtorc1_tx_merged_mtor_genes_all_3$symbol,
                selectLab = '', #mtorc1_tx_merged_mtor_genes$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.6,
                labCol = "red",
                xlim = c(-3, 3),
                ylim = c(0, 60),
                col=c("black","black","black","purple"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, mTOR-associated genes ',
                legendPosition = "none")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_mTOR_associated_genes_unlabeled.png", width=7, height=7, dpi=300)

EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$symbol,
                selectLab = '', #mtorc1_tx_merged_mtor_genes$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.6,
                labCol = "red",
                xlim = c(-3, 3),
                ylim = c(0, 60),
                col=c("grey","grey","grey","grey"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, mTOR-associated genes ',
                legendPosition = "none")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_mTOR_associated_genes_background.png", width=7, height=7, dpi=300)


## HIPPO ranked gene list
EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$symbol,
                selectLab = mtorc1_tx_merged_HIPPO_genes_2$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                #colCustom = keyvals,
                colAlpha = 1,
                pointSize = 1,
                labSize = 1.5,
                xlim = c(-3, 3),
                ylim = c(0, 60),
                labCol = "red",
                col=c("black","black","black","grey"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, HIPPO custom color over-ride',
                legendPosition = "none")
ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_ONLY_HIPPO_associated_genes.png", width=7, height=7, dpi=300)

#just plot the mTOR-associated genes to overlay with plot above
EnhancedVolcano(mtorc1_tx_merged_HIPPO_genes_all_3,
                lab = mtorc1_tx_merged_HIPPO_genes_all_3$symbol,
                selectLab = mtorc1_tx_merged_HIPPO_genes$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.7,
                labCol = "red",
                xlim = c(-3, 3),
                ylim = c(0, 60),
                col=c("black","black","black","maroon3"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, HIPPO-associated genes ',
                legendPosition = "none",
                #drawConnectors = TRUE,
                widthConnectors = 0.5,
                typeConnectors = "closed", endsConnectors = "first", lengthConnectors = unit(0.01, "npc"), colConnectors = "grey10")
ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_HIPPO_associated_genes_small_labels.png", width=7, height=7, dpi=300)

EnhancedVolcano(mtorc1_tx_merged_HIPPO_genes_all_3,
                lab = mtorc1_tx_merged_HIPPO_genes_all_3$symbol,
                selectLab = '', #mtorc1_tx_merged_mtor_genes$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.6,
                labCol = "red",
                xlim = c(-3, 3),
                ylim = c(0, 60),
                col=c("black","black","black","maroon3"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, HIPPO-associated genes ',
                legendPosition = "none")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_HIPPO_associated_genes_unlabeled.png", width=7, height=7, dpi=300)

EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$symbol,
                selectLab = '', #mtorc1_tx_merged_mtor_genes$symbol[1:20],
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 1,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.6,
                labCol = "red",
                xlim = c(-3, 3),
                ylim = c(0, 60),
                col=c("grey","grey","grey","grey"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, HIPPO-associated genes ',
                legendPosition = "none")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
ggsave("Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_HIPPO_associated_genes_background.png", width=7, height=7, dpi=300)


write.csv(mtorc1_tx_merged_HIPPO_genes_all_3, "mtorc1_tx_merged_HIPPO_genes.csv")
write.csv(mtorc1_tx_merged_mtor_genes_all_3, "mtorc1_tx_merged_mTOR_genes.csv")


#volcano plot with the known mTOR intermediates (intersection_sig_v_mtor)
EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$symbol,
                selectLab = intersection_sig_v_mtor,
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 2,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 1,
                labCol = "black",
                xlim = c(-20, 20),
                ylim = c(0, 270),
                col=c("black","black","black","grey"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, mTOR-associated genes ',
                legendPosition = "none")+
  theme( #add transparent background to do the merging on Adobe Illustrator
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
ggsave("Phospho_Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_mTOR_genes.png", width=7, height=7, dpi=300)

#volcano plot with the known HIPPO intermediates (intersection_sig_v_hippo)
EnhancedVolcano(Rapnutri_res_all,
                lab = Rapnutri_res_all$symbol,
                selectLab = intersection_sig_v_hippo,
                x="log2FoldChange",
                y="pvalue",
                pCutoff = 0.05,
                FCcutoff = 2,
                shape=15,
                colAlpha = 1,
                pointSize = 0.9,
                labSize = 0.9,
                labCol = "black",
                xlim = c(-20, 20),
                ylim = c(0, 270),
                col=c("black","black","black","grey"),
                title = 'DESeq2 results',
                subtitle = 'Raptor nutri vs vector nutri, HIPPO-associated genes ',
                legendPosition = "none")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

ggsave("Phospho_Rapnutri_Volcnolab_FCcutoff_2_p_val_0.05_HIPPO_genes.png", width=7, height=7, dpi=300)


#############################################################################################

#write.table(Rapnutri_res_all, file="/Users/holzergc/Desktop/RapnutriALL.csv")
#write.table(Rapnutri_res_sig, file="/Users/holzergc/Desktop/RapnutriSIG.csv")


#GO ANALYSIS
library(goseq)

Rapnutri_res_ALL <- arrange(Rapnutri_res_all, padj)
Rapnutri_res_ALL <- left_join(x = Rapnutri_res_ALL,y = grch38[,c("ensgene","symbol","description")], by = "symbol")

#add ensemble gene IDs that weren't converted
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='H0YIS7'] <- 'ENSG00000161939'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='P16402'] <- 'ENSG00000124575'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q02952'] <- 'ENSG00000131016.'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q147U1'] <- 'ENSG00000196605.'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='QARS1'] <- 'ENSG00000172053'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='B7Z6N1'] <- 'ENSP00000222329'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='MACROH2A1'] <- 'ENSG00000113648'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='O43379-4'] <- 'ENSP00000384792'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='H2AX'] <- 'ENSG00000188486'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='H13'] <- 'ENSG00000124575'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='P16402'] <- 'ENSG00000124575'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q59GA1'] <- 'ENSP00000371626'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q59FP3'] <- 'ENSP00000307908'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q59GA1'] <- 'ENSP00000371626'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q02952'] <- 'ENSP00000253332' #swiss-prot ID
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='O43379'] <- 'ENSG00000075702' #uniprotKB
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q59GA1'] <- 'ENSP00000371626' #trembl ID
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='GAGE4'] <- 'ENSP00000371133'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='PABIR1'] <- 'ENSG00000187866'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='A0A384N6H1'] <- 'ENSG00000167658'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='Q53GD7'] <- 'ENSP00000344149'
Rapnutri_res_ALL$ensgene[Rapnutri_res_ALL$symbol=='P07910'] <- 'ENSP00000451291'

#filter all the phosphopeptides that don't have a valid ENSEMBLE ID
Rapnutri_res_ALL <- filter(Rapnutri_res_ALL, !is.na(x = Rapnutri_res_ALL$ensgene)) #losing some phosphopeptides that are not convertable to ensgene IDs

#top 50 that have ensemble gene id's
Rapnutri_res_ALL <- Rapnutri_res_ALL[1:50,]

Rapnutri_res_ALL<-Rapnutri_res_ALL[!duplicated(Rapnutri_res_ALL$ensgene), ]
Rapnutri_res_ALL<-data_frame(Rapnutri_res_ALL)


diff<-Rapnutri_res_ALL
diff<-mutate(diff, de= padj<=0.05 & abs(log2FoldChange)>=0)

diff<-dplyr::select(diff, ensgene, de)
diff1<-diff
diff1<-diff1 %>% mutate_all(na_if,"TRUE")
diff1[is.na(diff1)] <- 1

#remove duplicate rownames
#diff2<-diff1[!duplicated(diff1$ensgene), ]

#important step
genes <- as.integer(Rapnutri_res_ALL$padj <= 0.05 & abs(Rapnutri_res_ALL$log2FoldChange)>=0 & !is.na(Rapnutri_res_ALL$padj)) 
names(genes)<-diff1$ensgene

#genes<-data_frame(genes)
#filter(genes, !is.na())
#rownames(genes)<-diff1$ensgene
#genes<-rownames_to_column(genes, var="ensgene")
#genes<-mutate(genes, ensgene=diff1$ensgene)
#genes<-dplyr::select(genes, -ensgene)
#genes<-data_frame(genes)

#not super useful but there nonetheless
genlens<-getlength(names(genes), "hg19", "ensGene")
names(genlens)<-diff1$ensgene

#useful again
pwf=nullp(genes,"hg19","ensGene")
head(pwf)

#GO analysis 
go.wall=goseq(pwf,"hg19","ensGene")
head(go.wall)

library(ggplot2)
colnames(go.wall)<-c("ID","over_represented_pvalue","under_represented_pvalue","Gene_Number","numInCat","GO_Term","Category")
b<-ggplot(go.wall[1:30,], aes(x=Gene_Number, y=GO_Term, size = Gene_Number, col = Category)) + geom_point(alpha=0.7)
b+ggtitle("Rap Nutri Gene Ontology")+xlab("Gene Number")+ylab("GO Term")+theme_bw()
ggsave("RapnutriGO_highres1.png", width=10, height=7, dpi=300)
ggsave("RapnutriGO_highres1_wide.png", width=12.5, height=7, dpi=300)


#EXPORTING DATA FRAMES FOR SUPPLEMENTAL FIGURES

#ACTsignormcounts <- ACTsignormcounts %>% relocate(Phospho_site, .before = symbol)
#ACTsignormcounts <- ACTsignormcounts %>% relocate(ensgene, .after = symbol)
#ACTsignormcounts$description <- NULL
#write.csv(ACTsignormcounts,"Rap_v_control_normalized_counts.csv", row.names = FALSE)

