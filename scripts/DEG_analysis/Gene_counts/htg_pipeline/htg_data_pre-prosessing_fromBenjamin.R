#Load libraries
setwd("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/RNAseqMEDBIOINFO/")
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


getwd()
library(org.Hs.eg.db)
library(enrichplot)
source("biostudies.R")
#Download data
if(!dir.exists("data")) dir.create("data")
files <- getBio(accession = "E-MTAB-2523", path = "data")
list.files("data")


#Import data
countTable <- read.table("data/EGAD00001000831_Total_Exon_Reads_The_paired_FF-FFPE_colon_set.txt",
                         as.is = TRUE, header = TRUE, sep = "\t", row.names = 1)  ### this "row.name =1" option makes the first column used as row.names in the genecount matrix
sampleTable <- read.table("data/E-MTAB-2523.sdrf.txt", as.is = TRUE, header = TRUE, sep = "\t")
#Inspect dataset
dim(countTable)

table(sampleTable$Characteristics.organism.part.)

table(sampleTable$Characteristics.individual.)

table(sampleTable$Characteristics.sex.)

table(sampleTable$Characteristics.disease.)

# specify study design
#Prepare sample annotation
sampleTable <- sampleTable %>%
  filter(
    Characteristics.organism.part. == "colon",
    grepl("FF$", Comment.ENA_ALIAS.),
    !duplicated(Comment.ENA_ALIAS.)) %>%
  transmute(
    sample = Comment.ENA_ALIAS.,
    individual = factor(Characteristics.individual.),
    sex = factor(Characteristics.sex.),
    disease = factor(recode(Characteristics.disease., `colon carcinoma` = "carcinoma")))
knitr::kable(sampleTable)

#Synchonize count data with sample table
countTable <- countTable[, -grep("FFPE", colnames(countTable))] ## exlude FFPE, only include FF, hence "-grep"
countTable <- countTable[, pmatch(sampleTable$sample, colnames(countTable))] ##
colnames(countTable) <- sampleTable$sample
countTable[1:5, 1:5]


#Filter low counts
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1))
hist(meanLog2CPM)
sum(meanLog2CPM <= 1)

countTable <- countTable[meanLog2CPM > 1, ]
dim(countTable)


#Filter low counts
head(countTable,6)
hist(countTable)
hist(rowMeans(countTable))
cT_cpm <- (cpm(countTable))
hist(cT_cpm)
hist(rowMeans(cT_cpm))

hist(rowMeans(log2(cpm(countTable))),xlim = c(-5,16)) ## QS: why no peak here
hist(rowMeans(log2(cpm(countTable) + 1)),xlim = c(-5,14))  ## ## but peak here (--< because we add one to all these
#all these cpm values of 0.)
## try to asnwer QS above
sum(-5<log2(cpm(countTable) <=1)) ## [1] 99763
sum(-5<log2(cpm(countTable) <=0))  ##[1] 53500
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1)) 
hist(meanLog2CPM) ## judging from the data, filter values less than 1
sum(meanLog2CPM <= 1) 

filtered_countTable <- countTable[meanLog2CPM > 1, ]
dim(filtered_countTable) ; dim(countTable)


#Prepare data for QC
dds <- DESeqDataSetFromMatrix(as.matrix(countTable),
                              design = ~ 0 + disease + individual + sex,
                              colData = sampleTable)


dds <- DESeqDataSetFromMatrix(as.matrix(countTable),
                              design = ~ 0 + disease + individual,
                              colData = sampleTable)
print(dds)

#Normalize
normCounts <- vst(dds, blind = TRUE)
assay(normCounts)[1:5, 1:5]

#Distribution
hist(assay(normCounts), main = "E-MTAB-2523")
# hist(assay(normCounts), main = "E-MTAB-2523",  breaks = seq(0, 20, by = 0.1))


#Sample heatmap
sampleDist <- cor(assay(normCounts), method = "spearman")
sampleColor <- brewer.pal(3, "Accent")[1:2]
names(sampleColor) <- levels(sampleTable$disease)
pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(Disease = sampleTable$disease,
                                     row.names = sampleTable$sample),
         annotation_colors = list(Disease = sampleColor))

#Sample PCA
pcaRes <- prcomp(t(assay(normCounts)))
varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
pcaDF <- data.frame(
  PC1 = pcaRes$x[, 1],
  PC2 = pcaRes$x[, 2],
  Disease = sampleTable$disease,
  Sample = sampleTable$sample)
pcaPlot <- ggplot(
  data = pcaDF,
  mapping = aes(x = PC1, y = PC2, color = Disease, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4) +
  labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
  scale_color_manual(values = brewer.pal(3, "Accent"))
print(pcaPlot)

#Remove outlier
countTable <- subset(countTable, select = -C6_T_FF)
sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
sampleTable <- droplevels(sampleTable)
colnames(countTable)


#Step 1: Define design matrix
designMatrix <- model.matrix(~ 0 + disease + individual, data = sampleTable)
designMatrix[1:5, 1:5]


#Step 2: Define contrast matrix
contrastMatrix <- makeContrasts(diseasecarcinoma - diseasenormal,
                                levels = designMatrix)
head(contrastMatrix)


#Step 3: Fit model
dge <- DGEList(countTable)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, designMatrix, robust = TRUE)
fit <- glmQLFit(dge, designMatrix, robust = TRUE)


#Step 4: Perform hypothesis testing
res <- glmQLFTest(fit, contrast = contrastMatrix)
res <- topTags(res, n = nrow(countTable))
sigRes <- subset(res$table, FDR < 0.05 & abs(logFC) > 1)
knitr::kable(head(sigRes))
nrow(sigRes)

#Visualize results
volcanoPlot <- ggplot(res$table,
                      aes(x = logFC, y = -log10(FDR),
                          color = ifelse(FDR < 0.05 & abs(logFC) > 1, "darkred", "grey"))) +
  geom_point() +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("Adjusted P value, Log"[10]*"")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("darkred", "grey", "steelblue")) +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = rownames(res$table)[1:10],
                      size = 2, color = "steelblue"),
                  data = res$table[1:10, ])
print(volcanoPlot)


#GO SEA
goSEA <- enrichGO(
  gene = rownames(sigRes),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)
#Convert symbols to Entrez IDs
sigRes$Entrez <- bitr(
  rownames(sigRes),
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db,
  drop = FALSE)$ENTREZID


#KEGG SEA
keggSEA <- enrichKEGG(
  gene = sigRes$Entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)


#GO GSEA
allFC <- sort(res$table$logFC, decreasing = TRUE)
names(allFC) <- rownames(res$table)
goGSEA <- gseGO(
  gene = allFC,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)


#Visulization
cnetplot(goSEA, colorEdge = TRUE, cex_label_gene = 0.5)
dotplot(goSEA)

# these do not work:
library(enrichplot)
goSEA <- pairwise_termsim(goSEA)
treeplot(goSEA)





