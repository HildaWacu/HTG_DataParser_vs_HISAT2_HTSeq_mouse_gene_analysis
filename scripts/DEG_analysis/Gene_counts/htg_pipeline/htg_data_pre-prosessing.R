## ----load results from count matrix input----------------------------------------------------------------------------
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/count_matrix_input.Rdata")

#Load libraries
# setwd("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/RNAseqMEDBIOINFO/")
library(dplyr)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(limma)
library(edgeR)
library(clusterProfiler)

#packageurl <- "https://cran.r-project.org/src/contrib/Archive/rvcheck/rvcheck_0.1.8.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#remove.packages("clusterProfiler")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("clusterProfiler")


#library(org.Hs.eg.db)## this is for Homo sapiens (Hs). Need one for Mus musculus (Mm)
library(enrichplot)
# source("biostudies.R")
#Download data
# if(!dir.exists("data")) dir.create("data")
# files <- getBio(accession = "E-MTAB-2523", path = "data")
# list.files("data")

#Import data
# countTable <- read.table("data/EGAD00001000831_Total_Exon_Reads_The_paired_FF-FFPE_colon_set.txt",
#                          as.is = TRUE, header = TRUE, sep = "\t", row.names = 1)  ### this "row.name =1" option makes the first column used as row.names in the genecount matrix
# sampleTable <- read.table("data/E-MTAB-2523.sdrf.txt", as.is = TRUE, header = TRUE, sep = "\t")

## inspect the data set
# htg_countTable <- htg_gene_counts_sh1

dim(htg_countTable)

table(coldata$condition1)

table(coldata$condition2) ## most impotant classification

table(coldata$condition3)

table(coldata$condition4)



#Inspect dataset
# dim(countTable)
# 
# table(sampleTable$Characteristics.organism.part.)
# 
# table(sampleTable$Characteristics.individual.)
# 
# table(sampleTable$Characteristics.sex.)
# 
# table(sampleTable$Characteristics.disease.)

# specify study design
#Prepare sample annotation

# sampleTable <- sampleTable %>%
#   filter(
#     Characteristics.organism.part. == "colon",
#     grepl("FF$", Comment.ENA_ALIAS.),
#     !duplicated(Comment.ENA_ALIAS.)) %>%
#   transmute(
#     sample = Comment.ENA_ALIAS.,
#     individual = factor(Characteristics.individual.),
#     sex = factor(Characteristics.sex.),
#     disease = factor(recode(Characteristics.disease., `colon carcinoma` = "carcinoma")))
# knitr::kable(sampleTable)
knitr::kable(coldata)

#Synchonize count data with sample table
# countTable <- countTable[, -grep("FFPE", colnames(countTable))] ## exlude FFPE, only include FF, hence "-grep"
# countTable <- countTable[, pmatch(sampleTable$sample, colnames(countTable))] ##
# colnames(countTable) <- sampleTable$sample
# countTable[1:5, 1:5]
# 
knitr::kable(unique(coldata))
unique(coldata)

# plot the hist of htg gene count to observe the distribution characteristics of my data
## firsly, convert to matrix
htg_TT <- as.matrix(htg_countTable)
## plot hist
hist(htg_TT, 
     xlim = c(-10, 200000),
     ylim = c(0,120000), 
     xlab = "Raw gene expression counts", 
     ylab = "No./frequency of all genes",
     main = "Distribution of all HTG raw gene counts")

## ploting the dist of gene counts for a single sample --> sample_25
ggplot(htg_countTable) +
  geom_histogram(aes(x = htg_countTable$sample_25), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes in sample_25")+
  labs(title = "Distribution of HTG raw gene count in sample 25")



## check the mean vs variance of the data
## for all
mean_htgcounts <- apply(htg_countTable, 1,mean)
length(mean_htgcounts)
var_htgcounts <- apply(htg_countTable,1,var) ;length(var_htgcounts)
plot(x = mean_htgcounts,y = var_htgcounts)
plot(x = log10(mean_htgcounts),y = log10(var_htgcounts))

df <- data.frame(mean_htgcounts,var_htgcounts)
ggplot(df) +
  geom_point(aes(x=mean_htgcounts, y=var_htgcounts)) + 
  scale_y_log10(limits = c(1,1e9)) + # uses a logarithmic scale to allow efficient spacing out of values
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

## for colitis&dysplasia+ve
mean_htgcounts <- apply(htg_countTable[,row.names(coldata)[1:10]], 1,mean)
length(mean_htgcounts)
var_htgcounts <- apply(htg_countTable[,row.names(coldata)[1:10]],1,var) ;length(var_htgcounts)
plot(x = mean_htgcounts,y = var_htgcounts)
plot(x = log10(mean_htgcounts),y = log10(var_htgcounts))

df <- data.frame(mean_htgcounts,var_htgcounts)
ggplot(df) +
  geom_point(aes(x=mean_htgcounts, y=var_htgcounts)) + 
  scale_y_log10(limits = c(1,1e11)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")


## for colitis&dysplasia-ve
mean_htgcounts <- apply(htg_countTable[,row.names(coldata)[11:20]], 1,mean)
length(mean_htgcounts)
var_htgcounts <- apply(htg_countTable[,row.names(coldata)[11:20]],1,var) ;length(var_htgcounts)


## for colitisOnly
mean_htgcounts <- apply(htg_countTable[,row.names(coldata)[21:42]], 1,mean)
length(mean_htgcounts)
var_htgcounts <- apply(htg_countTable[,row.names(coldata)[21:42]],1,var) ;length(var_htgcounts)

## for noColitis
mean_htgcounts <- apply(htg_countTable[,row.names(coldata)[43:64]], 1,mean)
length(mean_htgcounts)
var_htgcounts <- apply(htg_countTable[,row.names(coldata)[43:64]],1,var) ;length(var_htgcounts)



## Prefilter using log2(cpm)

#Filter low counts
head(htg_countTable,6)
# hist(htg_countTable)
# hist(rowMeans(htg_countTable))
# htg_cpm <- (cpm(htg_countTable))
# hist(htg_cpm)
# hist(rowMeans(htg_cpm)) 
# hist(rowMeans(log2(cpm(htg_countTable))),xlim = c(-5,16))

htg_meanLog2CPM <- rowMeans(log2(cpm(htg_countTable) + 1)) ## what is recommended in the script

#meanLog2CPM <- rowMeans(log2(cpm(htg_countTable) + 1))
hist(htg_meanLog2CPM, xlim = c(-5,20))
sum(htg_meanLog2CPM <= 1) ## [1] 6 --> only 6 lowly expressed genes. may not make a big difference

## therefore skipped this steps of filtering the 6 genes
cpm_filteredcountTable <- htg_countTable[htg_meanLog2CPM > 1, ]
dim(cpm_filteredcountTable)
cpm_filtered_htg_meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1))
hist(cpm_filtered_htg_meanLog2CPM, xlim = c(-5,20))

# countTable <- htg_countTable
# htg_countTable <- countTable
# dim(htg_countTable)

dds_cpm_filtered <- DESeqDataSetFromMatrix(as.matrix(cpm_filteredcountTable),
                                 design = ~ condition2,
                                 colData = coldata)
print(dds_cpm_filtered)
dds <- dds_cpm_filtered

## prefiltering step -- raw count filtering to retain only counts greater than or equal to 10
### before filtering'
ddsHTG <- DESeqDataSetFromMatrix(as.matrix(htg_countTable),
                                 design = ~ condition2,
                                 colData = coldata)

print(ddsHTG)
keep <- rowSums(counts(ddsHTG)) >= 10

(dds_rowSums_greaterthan_10 <- ddsHTG[keep,])
dim(dds_rowSums_greaterthan_10)  ## [1] 1659   64 --> non-filtered

## another way to perform the prefiltering in a single command
(dds_rowSums_greaterthan_10 <- dds[ rowSums(counts(dds)) >= 10, ])

dds <- dds_rowSums_greaterthan_10

#Prepare data for QC, i.e. have the DESeqDataSet object ready
print(dds)

ddsHTG
dds_cpm_filtered_HTG <- dds_cpm_filtered
dds_rowSums_greaterthan_10_HTG <- dds_rowSums_greaterthan_10

## Note: it is prefered in R that the first level of a factor be the reference level \
#(e.g. control, or untreated samples), so we can B4relevelB4 the condition2 factor like so: \
#(without re-leveling, DESeq2 will choose a reference level for factors based on the alphabetical order)

(dds$condition2 <- relevel(dds$condition2, "noColitis"))

save.image("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/htg_dds_filtering.Rdata")





# #Normalize --> 1. vst
# normCounts_vst <- vst(dds, blind = TRUE)
# assay(normCounts_vst)[1:5, 1:5] ## counts are not whole numbers as they are fractions of gene counts/median of ratios(size factor)
# 
# #Distribution
# hist(assay(normCounts_vst), main = "HTG Pipeline output when VST normalised",xlim = c(0,20))
# 
# hist(assay(normCounts_vst), main = "HTG Pipeline output when VST \
# normalised after log2(cpm) > 1 filtering ", xlim = c(0,20))
# 
# hist(assay(normCounts_vst), main = "HTG Pipeline output when VST \
# normalised after raw counts > 10 filtering ", xlim = c(0,20))
# 
# #Sample heatmap
# sampleDist <- cor(assay(normCounts_vst), method = "spearman")
# sampleColor <- brewer.pal(4, "Accent")[1:4]
# # sampleColor <- cols <- c( "red","orange", "#0080FF","lightgreen")
# # names(sampleColor) <- levels(sampleTable$disease)
# names(sampleColor) <- levels(coldata$condition2)
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = coldata$condition2,
#                                      row.names = row.names(coldata)),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          raw count >10 & VST normalised approach")
# 
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = coldata$condition2,
#                                      row.names = row.names(coldata)),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          log2(cpm) >1 & VST normalised approach")
# 
# #Sample PCA
# pcaRes <- prcomp(t(assay(normCounts_vst)))
# varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
# pcaDF <- data.frame(
#   PC1 = pcaRes$x[, 1],
#   PC2 = pcaRes$x[, 2],
#   Disease = coldata$condition2,
#   Sample = row.names(coldata))
# pcaPlot <- ggplot(
#   data = pcaDF,
#   mapping = aes(x = PC1, y = PC2, color = Disease, label = Sample)) +
#   geom_point(size = 3) +
#   geom_text_repel(size = 0) +
#   labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
#   scale_color_manual(values = cols)
# print(pcaPlot)
# 
# #Remove outlier
# countTable <- subset(countTable, select = -C6_T_FF)
# sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# sampleTable <- droplevels(sampleTable)
# colnames(countTable)
# 
# 
# 
# #Normalize --> 2. using rlog
# library(DESeq2)
# 
# # Estimate size factors
# dds <- ddsHTG
# dds <- estimateSizeFactors(dds)
# 
# # Perform rlog transformation
# normCounts_rlog <- rlog(dds)
# 
# assay(normCounts_rlog)[1:5, 1:5] ## counts are not whole numbers as they are fractions of gene counts/median of ratios(size factor)
# 
# #Distribution
# hist(assay(normCounts_rlog), main = "HTG Pipeline output when rlog normalised", xlim = c(-5,25))
# 
# hist(assay(normCounts_rlog), main = "HTG Pipeline output when rlog \
#      normalised-log2(cpm)>1-filtered counts", xlim = c(-5,25))
# 
# hist(assay(normCounts_rlog), main = "HTG Pipeline output when rlog \
#      normalised-raw counts>10-filtered counts", xlim = c(-5,25))
# 
# #Sample heatmap
# sampleDist <- cor(assay(normCounts_rlog), method = "spearman")
# # sampleColor <- brewer.pal(4, "Accent")[1:4]
# sampleColor <- cols <- c( "red","orange", "#0080FF","lightgreen")
# # names(sampleColor) <- levels(sampleTable$disease)
# names(sampleColor) <- levels(coldata$condition2)
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = coldata$condition2,
#                                      row.names = row.names(coldata)),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          raw count >10 & rlog normalised approach")
# 
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = coldata$condition2,
#                                      row.names = row.names(coldata)),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          log2(cpm) >1 & rlog normalised approach")
# 
# 
# #Sample PCA
# pcaRes <- prcomp(t(assay(normCounts_rlog)))
# varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
# pcaDF <- data.frame(
#   PC1 = pcaRes$x[, 1],
#   PC2 = pcaRes$x[, 2],
#   Disease = coldata$condition2,
#   Sample = row.names(coldata))
# pcaPlot <- ggplot(
#   data = pcaDF,
#   mapping = aes(x = PC1, y = PC2, color = Disease, label = Sample)) +
#   geom_point(size = 3) +
#   geom_text_repel(size = 0) +
#   labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
#   scale_color_manual(values = cols)
# print(pcaPlot)
# 
# #Remove outlier
# countTable <- subset(countTable, select = -C6_T_FF)
# sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# sampleTable <- droplevels(sampleTable)
# colnames(countTable)
# 
# 
# save.image(file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/htg_dds_QC.Rdata")
