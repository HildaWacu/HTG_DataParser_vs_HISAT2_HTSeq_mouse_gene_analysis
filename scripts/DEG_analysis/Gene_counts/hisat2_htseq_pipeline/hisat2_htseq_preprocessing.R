# Perform data preprocessing on the HISAT2-HTSeq dds
# Load the dds 
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/hisat2htseq_deseqdataset.Rdata")



## load packages
library(dplyr)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(limma)
library(edgeR)
library(clusterProfiler)


## inspect my data
table(sampleTable$sampleName)
table(sampleTable$fileName)
table(sampleTable$condition)
dim(counts(ddsHTSeq))

# plot the hist of htg gene count to observe the distribution characteristics of my data
## firsly, convert to matrix
htseq_TT <- as.matrix(counts(ddsHTSeq))
## plot hist
hist(htseq_TT, 
     xlim = c(-10, 200000),
     ylim = c(0,120000), 
     xlab = "Raw gene expression counts", 
     ylab = "No./frequency of all genes",
     main = "Distribution of all HISAT2-HTSeq raw gene counts")
# 
# htseq_TT <- as.matrix(counts(dds_rowSums_greaterthan_10))
# ## plot hist
# hist(htseq_TT, 
#      xlim = c(-10, 200000),
#      ylim = c(0,120000), 
#      xlab = "Raw gene expression counts", 
#      ylab = "No./frequency of all genes",
#      main = "Distribution of all HISAT2-HTSeq raw gene counts     \
#      when counts >10 filtered")



# ## ploting the dist of gene counts for a single sample --> sample_25 !!not working yet
# htseq_counts <- as.matrix(counts(ddsHTSeq))
# ggplot(counts(ddsHTSeq)) +
#   geom_histogram(aes(x = ...), stat = "bin", bins = 200) +
#   xlab("Raw expression counts") +
#   ylab("Number of genes in sample_25")+
#   labs(title = "Distribution of HTG raw gene count in sample 25")


## prefiltering step -- raw count filtering to retain only counts greater than or equal to 10
keep <- rowSums(counts(ddsHTSeq)) >= 10

(dds_rowSums_greaterthan_10 <- ddsHTSeq[keep,])

## another way to perform the prefiltering
(dds_rowSums_greaterthan_10 <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) >= 10, ])

## new dim after prefiltering
dim(counts(dds_rowSums_greaterthan_10))

#Prepare data for QC after method 1 prefiltering --> input for QC processed is dds
print(dds_rowSums_greaterthan_10)
dds <- dds_rowSums_greaterthan_10

## Prefilter using log2(cpm)
#Filter low counts

htseq_countTable <- counts(ddsHTSeq)
head(htseq_countTable,6)

# hist(htseq_countTable)
# hist(rowMeans(htseq_countTable))
# htseq_cpm <- (cpm(htseq_countTable))
# hist(htseq_cpm)
# hist(rowMeans(htseq_cpm)) 
# 
# hist(rowMeans(log2(cpm(htseq_countTable))),xlim = c(-5,16))

htseq_meanLog2CPM <- rowMeans(log2(cpm(htseq_countTable) + 1)) ## what is recommended in the script

#meanLog2CPM <- rowMeans(log2(cpm(htg_countTable) + 1))
hist(htseq_meanLog2CPM, xlim = c(0,20), ylim = c(0,800))
sum(htseq_meanLog2CPM <= 1)

cpm_filteredcountTable <- htseq_countTable[htseq_meanLog2CPM > 1, ]
dim(cpm_filteredcountTable)    


cpm_filtered_htseq_meanLog2CPM <- rowMeans(log2(cpm(cpm_filteredcountTable) + 1)) ## what is recommended in the script

#meanLog2CPM <- rowMeans(log2(cpm(htg_countTable) + 1))
hist(cpm_filtered_htseq_meanLog2CPM, xlim = c(-5,20))
sum(cpm_filtered_htseq_meanLog2CPM <= 1)

#Prepare data for QC after method 2 prefiltering
dds_cpm_filtered <- DESeqDataSetFromMatrix(as.matrix(cpm_filteredcountTable),
                              design = ~ condition,
                              colData = sampleTable)


print(dds_cpm_filtered)
dds <- dds_cpm_filtered

## for comparison purposed, rename the dds
ddsHTSeq
dds_cpm_filtered_HTSeq <- dds_cpm_filtered
dds_rowSums_greaterthan_10_HTSeq <- dds_rowSums_greaterthan_10

dds
## Note: it is prefered in R that the first level of a factor be the reference level \
#(e.g. control, or untreated samples), so we can B4relevelB4 the condition2 factor like so: \
#(without re-leveling, DESeq2 will choose a reference level for factors based on the alphabetical order)

(dds$condition <- relevel(dds$condition, "noColitis"))


save.image("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/htseq_dds_filtering.Rdata")

# 
# #Normalize --> 1. vst
# normCounts_vst <- vst(dds, blind = TRUE)
# assay(normCounts_vst)[1:5, 1:5] ## counts are not whole numbers as they are fractions of gene counts/median of ratios(size factor)
# 
# #Distribution
# hist(assay(normCounts_vst), main = "HISAT2-HTSeq Pipeline output when VST normalised", xlim = c(0,20))
# 
# hist(assay(normCounts_vst), main = "HISAT2-HTSeq Pipeline output when VST \
# normalised after log2(cpm) > 1 filtering ", xlim = c(0,20))
# 
# hist(assay(normCounts_vst), main = "HISAT2-HTSeq Pipeline output when VST \
# normalised after raw counts > 10 filtering ", xlim = c(0,20))
# 
# 
# #Sample heatmap
# sampleDist <- cor(assay(normCounts_vst), method = "spearman")
# # sampleColor <- brewer.pal(4, "Accent")[1:4]
# sampleColor <- cols <- c( "red","orange", "#0080FF","lightgreen")
# names(sampleColor) <- levels(sampleTable$condition)
# # sample_names <- (sub("(sample_[0-9]+).*","\\1",sampleTable$sampleName))
# sampleTable$condition <- factor(sampleTable$condition)
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = sampleTable$condition,
#                                      row.names = sampleTable$sampleName),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          log2(cpm) >1 & VST normalised approach")
# 
# 
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = sampleTable$condition,
#                                      row.names = sampleTable$sampleName),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          raw count >10 & VST normalised approach")
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
# # #Remove outlier
# # countTable <- subset(countTable, select = -C6_T_FF)
# # sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# # sampleTable <- droplevels(sampleTable)
# # colnames(countTable)
# 
# 
# 
# #Normalize --> 2. using rlog
# library(DESeq2)
# 
# # Estimate size factors -#### may be wrong since rlog incorporates the log transformation
# dds_for_rlog <- estimateSizeFactors(dds)
# 
# # Perform rlog transformation
# normCounts_rlog <- rlog(dds_for_rlog)
# 
# assay(normCounts_rlog)[1:5, 1:5] ## counts are not whole numbers as they are fractions of gene counts/median of ratios(size factor)
# 
# #Distribution
# hist(assay(normCounts_rlog), main = "HISAT2-HTSeq Pipeline output when rlog \
#      normalised-log2(cpm)>1-filtered counts", xlim = c(-5,25))
# 
# hist(assay(normCounts_rlog), main = "HISAT2-HTSeq Pipeline output when rlog \
#      normalised-raw counts>10-filtered counts", xlim = c(-5,25))
# 
# #Sample heatmap
# sampleDist <- cor(assay(normCounts_rlog), method = "spearman")
# # sampleColor <- brewer.pal(4, "Accent")[1:4]
# sampleColor <- cols <- c( "red","orange", "#0080FF","lightgreen")
# # names(sampleColor) <- levels(sampleTable$disease)
# names(sampleColor) <- levels(sampleTable$condition)
# sample_names <- (sub("(sample_[0-9]+).*","\\1",sampleTable$sampleName))
# sampleTable$condition <- factor(sampleTable$condition)
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = sampleTable$condition,
#                                      row.names = sampleTable$sampleName),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          log2(cpm) >1 & rlog normalised approach")
# 
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(Disease = sampleTable$condition,
#                                      row.names = sampleTable$sampleName),
#          annotation_colors = list(Disease = sampleColor),
#          main = "Heatmap after filtering by \
#          raw count >10 & rlog normalised approach")
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
# 
# # #Normalize
# # normCounts <- vst(dds, blind = TRUE)
# # assay(normCounts)[1:5, 1:5] ## counts are not whole numbers as they are fractions of gene counts/median of ratios(size factor)
# # 
# # #Distribution
# # hist(assay(normCounts), main = "HISAT2-HTSeq Pipeline output")
# # 
# # 
# # #Sample heatmap
# # sampleDist <- cor(assay(normCounts), method = "spearman")
# # sampleColor <- brewer.pal(4, "Accent")[1:4]
# # # names(sampleColor) <- levels(sampleTable$disease)
# # names(sampleColor) <- levels(sampleTable$condition)
# # sample_names <- (sub("(sample_[0-9]+).*","\\1",sampleTable$sampleName))
# # sampleTable$condition <- factor(sampleTable$condition)
# # pheatmap(sampleDist,
# #          clustering_distance_rows = as.dist(1 - sampleDist),
# #          clustering_distance_cols = as.dist(1 - sampleDist),
# #          annotation_col = data.frame(Disease = sampleTable$condition,
# #                                      row.names = sampleTable$sampleName),
# #          annotation_colors = list(Disease = sampleColor))
# # 
# # #Sample PCA
# # pcaRes <- prcomp(t(assay(normCounts)))
# # varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
# # pcaDF <- data.frame(
# #   PC1 = pcaRes$x[, 1],
# #   PC2 = pcaRes$x[, 2],
# #   Disease = sampleTable$condition,
# #   Sample = sample_names)
# # pcaPlot <- ggplot(
# #   data = pcaDF,
# #   mapping = aes(x = PC1, y = PC2, color = Disease, label = Sample )) +
# #   geom_point(size = 3) +
# #   geom_text_repel(size = 0) +
# #   labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
# #   theme_minimal() +
# #   theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
# #   scale_color_manual(values = brewer.pal(4, "Accent"))
# # print(pcaPlot)
# # 
# # #Remove outlier
# # countTable <- subset(countTable, select = -C6_T_FF)
# # sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# # sampleTable <- droplevels(sampleTable)
# # colnames(countTable)
# # 
# # ### filtered method 2 (log2(cpm))
# # #Normalize
# # normCounts <- vst(ddsHTSeq, blind = TRUE)
# # assay(normCounts)[1:5, 1:5] ## counts are not whole numbers as they are fractions of gene counts/median of ratios(size factor)
# # 
# # #Distribution
# # hist(assay(normCounts), main = "HISAT2-HTSeq Pipeline output")
# # 
# # 
# # #Sample heatmap
# # sampleDist <- cor(assay(normCounts), method = "spearman")
# # # sampleColor <- brewer.pal(4, "Accent")[1:4]
# # sampleColor <- cols <- c( "red","orange", "#0080FF","lightgreen")
# # 
# # # sampleColor <- cols <- c("#FF5500","#FFAA00","#FFFF00", "#AAFF00")
# # # "#0080FF","#0000FF" "#8000FF" "#FF00FF" "#FF0080""
# # # 
# # # "#FF0000" 
# # # "#FF5500" "#FFAA00" "#FFFF00" "#AAFF00" "#55FF00" "#00FF00" "#00FF55"
# # # "#00FFAA" "#00FFFF" "#00AAFF" "#0055FF" "#0000FF" "#5500FF" "#AA00FF" "#FF00FF"
# # # "#FF00AA" "#FF0055"
# # # names(sampleColor) <- levels(sampleTable$disease)
# # names(sampleColor) <- levels(coldata$condition2)
# # pheatmap(sampleDist,
# #          clustering_distance_rows = as.dist(1 - sampleDist),
# #          clustering_distance_cols = as.dist(1 - sampleDist),
# #          annotation_col = data.frame(Disease = coldata$condition2,
# #                                      row.names = row.names(coldata))
# #          ,         annotation_colors = list(Disease = sampleColor))
# # 
# # #Sample PCA
# # pcaRes <- prcomp(t(assay(normCounts)))
# # 
# # varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
# # pcaDF <- data.frame(
# #   PC1 = pcaRes$x[, 1],
# #   PC2 = pcaRes$x[, 2],
# #   Disease = coldata$condition2,
# #   Sample = row.names(coldata))
# # 
# # pcaPlot <- ggplot(
# #   data = pcaDF,
# #   mapping = aes(x = PC1, y = PC2, color = Disease, label = Sample)) +
# #   geom_point(size = 2.5) +
# #   geom_text_repel(size = 0.1, max.overlaps = 100) +
# #   labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
# #   theme_minimal() +
# #   theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
# #   scale_color_manual(values = cols)
# # print(pcaPlot)
# # 
# 
# #Remove outlier
# # countTable <- subset(countTable, select = -C6_T_FF)
# # sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# # sampleTable <- droplevels(sampleTable)
# # colnames(countTable)
# 
# save.image(file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/hisat2htseq_dds_QC.Rdata")
