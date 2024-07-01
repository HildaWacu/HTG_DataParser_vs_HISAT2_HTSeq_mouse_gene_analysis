set.seed(1)

sample_df <- data.frame(
  group = factor(rep(letters[1:3], each = 10)),
  value = rnorm(30)
)

group_means_df <- setNames(
  aggregate(value ~ group, sample_df, mean),
  c("group", "group_mean")
)

# The following three code blocks create the same graphic, each using one
# of the three patterns specified above. In each graphic, the sample data
# are plotted in the first layer and the group means data frame is used to
# plot larger red points on top of the sample data in the second layer.

# Pattern 1
# Both the `data` and `mapping` arguments are passed into the `ggplot()`
# call. Those arguments are omitted in the first `geom_point()` layer
# because they get passed along from the `ggplot()` call. Note that the
# second `geom_point()` layer re-uses the `x = group` aesthetic through
# that mechanism but overrides the y-position aesthetic.
ggplot(data = sample_df, mapping = aes(x = group, y = value)) +
  geom_point() +
  geom_point(
    mapping = aes(y = group_mean), data = group_means_df,
    colour = 'red', size = 3
  )

# Pattern 2
# Same plot as above, passing only the `data` argument into the `ggplot()`
# call. The `mapping` arguments are now required in each `geom_point()`
# layer because there is no `mapping` argument passed along from the
# `ggplot()` call.
ggplot(data = sample_df) +
  geom_point(mapping = aes(x = group, y = value)) +
  geom_point(
    mapping = aes(x = group, y = group_mean), data = group_means_df,
    colour = 'red', size = 3
  )

# Pattern 3
# Same plot as above, passing neither the `data` or `mapping` arguments
# into the `ggplot()` call. Both those arguments are now required in
# each `geom_point()` layer. This pattern can be particularly useful when
# creating more complex graphics with many layers using data from multiple
# data frames.
ggplot() +
  geom_point(mapping = aes(x = group, y = value), data = sample_df) +
  geom_point(
    mapping = aes(x = group, y = group_mean), data = group_means_df,
    colour = 'red', size = 3
  )


####======================================================================================================================
---
  title: "GSEA"
author: "Hildah Njoroge"
date: "2023-08-30"
output: html_document
---
  ## https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#enrichplot 
  # Step 1
  
  ### Load the SigRes data
  
  ```{r}
#HTseq data
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/deseq2_statistical_analysis_HTSeq.Rdata")

#HTG data
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/deseq2_statistical_analysis_HTG.Rdata")

```

# Step 2

### Install and load the necessary packages

```{r}
if (!require(clusterProfiler, quietly = T)) {install.packages("clusterProfiler")}
library(clusterProfiler)
if (!require(org.Mm.eg.db, quietly = T)) {BiocManager::install("org.Mm.eg.db")}
library(org.Mm.eg.db)
library(org.Hs.eg.db)
```

# Step 3

### Run GO-SEA

##### A. WT vs colOnly

##### HTSeq

```{r}
#GO SEA
goSEA_HTseq <- enrichGO(
  gene = rownames(sigRes_HTSeq),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)
```

# Step 4 Visualisation

```{r}
# cnetplot(goSEA, colorEdge = TRUE, cex_label_gene = 1)
# dotplot(goSEA)
# 

library(enrichplot)
cnetplot(goSEA_HTseq, 
         colorEdge = TRUE, 
         cex_label_gene = 0.5, 
         node_label = "all",
         )
dotplot(goSEA_HTseq)

barplot(goSEA_HTseq, 
        title = "Enrichment in the Colitis Only Group - HTSeq", 
        showCategory = 10,
        font.size = 16) ## to display the top 10 pathways


goSEA_HTseq1 <- pairwise_termsim(goSEA_HTseq)
treeplot(goSEA_HTseq1)
         # ,
         # label_format = NULL,
         # label_format_cladelab = 30,
         # label_format_tiplab = NULL,
         # offset_tiplab = rel(1),
         # offset = rel(5))

# goSEA_HTseq[1:10,1:7]

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTseq)
paste("Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTseq[1:10,1:7]
results_goSEA_HTseq <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTseq <- t(results_goSEA_HTseq)
colnames(results_goSEA_HTseq) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTseq <- results_goSEA_HTseq[-1,]
results_goSEA_HTseq <- goSEA_HTseq[1:10,1:7]
```


##### HTG WT vs colOnly

```{r}
#GO SEA
goSEA_HTG <- enrichGO(
  gene = rownames(sigRes_HTG),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)
```

# Step 4 Visualisation

```{r}
# cnetplot(goSEA, colorEdge = TRUE, cex_label_gene = 0.9)
# dotplot(goSEA)


cnetplot(goSEA_HTG, 
         colorEdge = T, 
         cex_label_gene = 0.5,
         node_label = "all"
         )

cnetplot(goSEA_HTG, 
         colorEdge = T, 
         cex_label_gene = 0.5,
         node_label = "all",
         showCategory = 10
)


dotplot(goSEA_HTG)

barplot(goSEA_HTG, 
        title = "Enrichment in the Colitis Only Group - HTG", 
        showCategory = 10,
        font.size = 16) ## to display the top 10 pathways


library(enrichplot)
goSEA_HTG1 <- pairwise_termsim(goSEA_HTG)
treeplot(goSEA_HTG1,
         label_format = NULL,
         label_format_cladelab = 30,
         label_format_tiplab = NULL,
         offset_tiplab = rel(1),
         offset = rel(5))

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTG)
paste("WT vs ColOnly, HTG - Number of pathways differentially regulated is", nrow(goSEA_HTG))
paste("WT vs ColOnly, HTSeq - Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTG[1:10,1:7]
results_goSEA_HTG <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTG <- t(results_goSEA_HTG)
colnames(results_goSEA_HTG) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTG <- results_goSEA_HTG[-1,]
results_goSEA_HTG <- goSEA_HTG[1:10,1:7]
```

```{r}
results1 <- data.frame(top10pathways_HTG=results_goSEA_HTG[,1:2],top10pathways_HTSeq=results_goSEA_HTseq[,1:2])
results1
write.csv(results1,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/top10_diff_reg_pathways_nocol_colonly.csv")
setdiff(results_goSEA_HTG[,2],results_goSEA_HTseq[,2]) ## meaning what is in the first that isn't in the second set
setdiff(results_goSEA_HTseq[,2],results_goSEA_HTG[,2])
```

### ggplot for visualisation

```{r}
library(ggplot2)
library(ggrepel)
library(stringr)
# df <- goSEA_HTG@result[1:10,] 
df <- results_goSEA_HTG[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)

### Add the theme section to y axes label more visible - adjust this accordingly 
my_theme <- theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(df, 
       aes(x = Negative_log10_adjusted_PValue, 
           y = reorder(Pathway_Description, 
                       Negative_log10_adjusted_PValue), 
           label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(colour="black",
                  size=4,
                  min.segment.length = 10,
                  hjust = 0.5, 
                  vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", 
       y = "Pathway Description (HTG)") +
  my_theme


df <- results_goSEA_HTseq[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTSeq)")+
  my_theme

### Add the theme sction when you want to plot the graphs with visible axes label
# theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold"))

```

# Step 3

### Run GO-SEA

##### HTSeq WT vs Colnodysp

```{r}
#GO SEA
goSEA_HTseq_Colnodysp <- enrichGO(
  gene = rownames(sigRes_HTSeq_Colnodysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

goSEA_HTseq <- goSEA_HTseq_Colnodysp
```

# Step 4 Visualisation

```{r}
# cnetplot(goSEA, colorEdge = TRUE, cex_label_gene = 0.9)
# dotplot(goSEA)


cnetplot(goSEA_HTseq_Colnodysp, colorEdge = T, cex_label_gene = 0.5)
cnetplot(goSEA_HTseq_Colnodysp, colorEdge = T, cex_label_gene = 0.5, showCategory = 10)
dotplot(goSEA_HTseq_Colnodysp)

barplot(goSEA_HTseq_Colnodysp, 
        title = "Enrichment in the Colitis and Dysplasia (-) Group - HTSeq", 
        showCategory = 10,
        font.size = 16) ## to display the top 10 pathways


library(enrichplot)
goSEA_HTseq_Colnodysp1 <- pairwise_termsim(goSEA_HTseq_Colnodysp)
treeplot(goSEA_HTseq_Colnodysp1)

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTseq)
paste("Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTseq[1:10,1:7]
results_goSEA_HTseq <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTseq <- t(results_goSEA_HTseq)
colnames(results_goSEA_HTseq) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTseq <- results_goSEA_HTseq[-1,]
results_goSEA_HTseq <- goSEA_HTseq[1:10,1:7]
```

##### HTG WT vs Colnodysp

```{r}
#GO SEA
goSEA_HTG_Colnodysp <- enrichGO(
  gene = rownames(sigRes_HTG_Colnodysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

goSEA_HTG <- goSEA_HTG_Colnodysp

```

cnetplot(goSEA_HTG_Colnodysp, colorEdge = T, cex_label_gene = 0.5)
cnetplot(goSEA_HTG_Colnodysp, colorEdge = T, cex_label_gene = 0.5, showCategory = 10)
dotplot(goSEA_HTG_Colnodysp)

barplot(goSEA_HTG_Colnodysp, 
        title = "Enrichment in the Colitis and Dysplasia (-) Group - HTG", 
        showCategory = 10,
        font.size = 16) ## to display the top 10 pathways


library(enrichplot)
goSEA_HTG_Colnodysp1 <- pairwise_termsim(goSEA_HTG_Colnodysp)
treeplot(goSEA_HTG_Colnodysp1)


## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTG)
paste("WT vs Colitis without dysp, HTG - Number of pathways differentially regulated is", nrow(goSEA_HTG))
paste("WT vs Colitis without dysp, HTSeq - Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTG[1:10,1:7]
results_goSEA_HTG <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTG <- t(results_goSEA_HTG)
colnames(results_goSEA_HTG) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTG <- results_goSEA_HTG[-1,]
results_goSEA_HTG <- goSEA_HTG[1:10,1:7]
```

### Then, compare the top 10 differentially regulated pathways in the nocol vs colnodysp groups

```{r}
results2 <- data.frame(top10pathways_HTG=results_goSEA_HTG[,1:2],top10pathways_HTSeq=results_goSEA_HTseq[,1:2])
results2
write.csv(results2,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/top10_diff_reg_pathways_nocol_colnodysp.csv")
setdiff(results_goSEA_HTG[,2],results_goSEA_HTseq[,2])
setdiff(results_goSEA_HTseq[,2],results_goSEA_HTG[,2])
```

### ggplot for visualisation

```{r}
library(stringr)
# df <- goSEA_HTG@result[1:10,] 
df <- results_goSEA_HTG[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTG)") +
  my_theme


df <- results_goSEA_HTseq[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTSeq)") +
  my_theme

```

# Step 3

### Run GO-SEA

##### HTSeq WT vs Coldysp

```{r}
#GO SEA
goSEA_HTseq_Coldysp <- enrichGO(
  gene = rownames(sigRes_HTSeq_Coldysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)


goSEA_HTseq <- goSEA_HTseq_Coldysp

```

cnetplot(goSEA_HTseq_Coldysp, colorEdge = T, cex_label_gene = 0.5)
dotplot(goSEA_HTseq_Coldysp)

barplot(goSEA_HTseq_Coldysp, 
        title = "Enrichment in the Colitis and Dysplasia (+) Group - HTSeq", 
        showCategory = 10,
        font.size = 16) ## to display the top 10 pathways


library(enrichplot)
goSEA_HTseq_Coldysp1 <- pairwise_termsim(goSEA_HTseq_Coldysp)
treeplot(goSEA_HTseq_Coldysp1)


## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTseq)
paste("Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTseq[1:10,1:7]
results_goSEA_HTseq <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTseq <- t(results_goSEA_HTseq)
colnames(results_goSEA_HTseq) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTseq <- results_goSEA_HTseq[-1,]
results_goSEA_HTseq <- goSEA_HTseq[1:10,1:7]
```

##### HTG WT vs Coldysp

```{r}
#GO SEA
goSEA_HTG_Coldysp <- enrichGO(
  gene = rownames(sigRes_HTG_Coldysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

goSEA_HTG <- goSEA_HTG_Coldysp

```

cnetplot(goSEA_HTG_Coldysp, colorEdge = T, cex_label_gene = 0.5)
dotplot(goSEA_HTG_Coldysp)

barplot(goSEA_HTG_Coldysp, 
        title = "Enrichment in the Colitis and Dysplasia (+) Group - HTG", 
        showCategory = 10,
        font.size = 16) ## to display the top 10 pathways

library(enrichplot)
goSEA_HTG_Coldysp1 <- pairwise_termsim(goSEA_HTG_Coldysp)
treeplot(goSEA_HTG_Coldysp1)



## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTG)
paste("WT vs Colitis with dysp, HTG - Number of pathways differentially regulated is", nrow(goSEA_HTG))
paste("WT vs Colitis with dysp, HTSeq - Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTG[1:10,1:7]
results_goSEA_HTG <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTG <- t(results_goSEA_HTG)
colnames(results_goSEA_HTG) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTG <- results_goSEA_HTG[-1,]
results_goSEA_HTG <- goSEA_HTG[1:10,1:7]



```

### Then, compare the top 10 differentially regulated pathways in the nocol vs colnodysp groups

```{r}
results3 <- data.frame(top10pathways_HTG=results_goSEA_HTG[,1:2],top10pathways_HTSeq=results_goSEA_HTseq[,1:2])
results3
write.csv(results3,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/top10_diff_reg_pathways_nocol_coldysp.csv")
setdiff(results_goSEA_HTG[,2],results_goSEA_HTseq[,2])
setdiff(results_goSEA_HTseq[,2],results_goSEA_HTG[,2])
```

### ggplot for visualisation

```{r}

library(stringr)
# df <- goSEA_HTG@result[1:10,] 
df <- results_goSEA_HTG[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTG)") +
  my_theme


df <- results_goSEA_HTseq[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTSeq)") +
  my_theme

```

# Step 3

### Run GO-SEA

##### HTSeq

#### ColOnly vs Colnodysp

```{r}
#GO SEA
goSEA_HTseq_2_Colnodysp <- enrichGO(
  gene = rownames(sigRes_HTSeq_2_Colnodysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)


goSEA_HTseq <- goSEA_HTseq_2_Colnodysp

```

```{r}
cnetplot(goSEA_HTseq_2_Colnodysp, colorEdge = T, cex_label_gene = 0.5)
dotplot(goSEA_HTseq_2_Colnodysp)

barplot(goSEA_HTseq_2_Colnodysp, title = "Enrichment in the Colitis and Dysplasia (-) Group against Colitis Only - HTSeq", showCategory = 10) ## to display the top 10 pathways

library(enrichplot)
goSEA_HTseq_2_Colnodysp1 <- pairwise_termsim(goSEA_HTseq_2_Colnodysp)
treeplot(goSEA_HTseq_2_Colnodysp1)
```




## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTseq)
paste("Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTseq[1:10,1:7]
results_goSEA_HTseq <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTseq <- t(results_goSEA_HTseq)
colnames(results_goSEA_HTseq) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTseq <- results_goSEA_HTseq[-1,]
results_goSEA_HTseq <- goSEA_HTseq[1:10,1:7]
```

##### HTG ColOnly vs Colnodysp

```{r}
#GO SEA
goSEA_HTG_2_Colnodysp <- enrichGO(
  gene = rownames(sigRes_HTG_2_Colnodysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

goSEA_HTG <- goSEA_HTG_2_Colnodysp

```


```{r}
cnetplot(goSEA_HTG_2_Colnodysp, colorEdge = T, cex_label_gene = 0.5, cex_label_category = 0.5)
dotplot(goSEA_HTG_2_Colnodysp)

barplot(goSEA_HTG_2_Colnodysp, title = "Enrichment in the Colitis and Dysplasia (-) Group against Colitis Only - HTG", showCategory = 10) ## to display the top 10 pathways

library(enrichplot)
goSEA_HTG_2_Colnodysp1 <- pairwise_termsim(goSEA_HTG_2_Colnodysp)
treeplot(goSEA_HTG_2_Colnodysp1)
```


## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTG)
paste("ColOnly vs Colitis without dysp,, HTG - Number of pathways differentially regulated is", nrow(goSEA_HTG))
paste("ColOnly vs Colitis without dysp, HTSeq - Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTG[1:10,1:7]
results_goSEA_HTG <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTG <- t(results_goSEA_HTG)
colnames(results_goSEA_HTG) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTG <- results_goSEA_HTG[-1,]
results_goSEA_HTG <- goSEA_HTG[1:10,1:7]
```

### Then, compare the top 10 differentially regulated pathways in the nocol vs colnodysp groups

```{r}
results3 <- data.frame(top10pathways_HTG=results_goSEA_HTG[,2],top10pathways_HTSeq=results_goSEA_HTseq[,2])
results3
write.csv(results3,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/top10_diff_reg_pathways_nocol_coldysp.csv")
setdiff(results_goSEA_HTG[,2],results_goSEA_HTseq[,2])
setdiff(results_goSEA_HTseq[,2],results_goSEA_HTG[,2])
```

### ggplot for visualisation

```{r}

library(stringr)
# df <- goSEA_HTG@result[1:10,] 
df <- results_goSEA_HTG[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=3,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTG)") +
  my_theme


df <- results_goSEA_HTseq[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTSeq)") +
  my_theme

```

# Step 3

### Run GO-SEA

##### HTSeq

#### ColOnly vs Coldysp

```{r}
#GO SEA
goSEA_HTseq_2_Coldysp <- enrichGO(
  gene = rownames(sigRes_HTSeq_2_Coldysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)


goSEA_HTseq <- goSEA_HTseq_2_Coldysp

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTseq)
paste("Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTseq[1:10,1:7]
results_goSEA_HTseq <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTseq <- t(results_goSEA_HTseq)
colnames(results_goSEA_HTseq) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTseq <- results_goSEA_HTseq[-1,]
results_goSEA_HTseq <- goSEA_HTseq[1:10,1:7]
```

##### HTG ColOnly vs Coldysp

```{r}
#GO SEA
goSEA_HTG_2_Coldysp <- enrichGO(
  gene = rownames(sigRes_HTG_2_Coldysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

goSEA_HTG <- goSEA_HTG_2_Coldysp

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTG)
paste("ColOnly vs Colitis with dysp,, HTG - Number of pathways differentially regulated is", nrow(goSEA_HTG))
paste("ColOnly vs Colitis with dysp, HTSeq - Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTG[1:10,1:7]
results_goSEA_HTG <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTG <- t(results_goSEA_HTG)
colnames(results_goSEA_HTG) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTG <- results_goSEA_HTG[-1,]
results_goSEA_HTG <- goSEA_HTG[1:10,1:7]
```

### Then, compare the top 10 differentially regulated pathways in the nocol vs colnodysp groups

```{r}
results3 <- data.frame(top10pathways_HTG=results_goSEA_HTG[,2],top10pathways_HTSeq=results_goSEA_HTseq[,2])
results3
write.csv(results3,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/top10_diff_reg_pathways_nocol_coldysp.csv")
setdiff(results_goSEA_HTG[,2],results_goSEA_HTseq[,2])
setdiff(results_goSEA_HTseq[,2],results_goSEA_HTG[,2])
```

### ggplot for visualisation

```{r}

library(stringr)
# df <- goSEA_HTG@result[1:10,] 
df <- results_goSEA_HTG[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTG)") +
  my_theme


df <- results_goSEA_HTseq[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTSeq)") +
  my_theme

```

# Step 3

### Run GO-SEA

##### HTSeq Colnodysp vs Coldysp

```{r}
#GO SEA
goSEA_HTseq_3_Colnodysp <- enrichGO(
  gene = rownames(sigRes_HTSeq_3_Colnodysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)


goSEA_HTseq <- goSEA_HTseq_3_Colnodysp

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTseq)
paste("Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTseq[1:10,1:7]
results_goSEA_HTseq <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTseq <- t(results_goSEA_HTseq)
colnames(results_goSEA_HTseq) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTseq <- results_goSEA_HTseq[-1,]
results_goSEA_HTseq <- goSEA_HTseq[1:10,1:7]
```

##### HTG WT vs Coldysp

```{r}
#GO SEA
goSEA_HTG_3_Coldysp <- enrichGO(
  gene = rownames(sigRes_HTG_3_Colnodysp),
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

goSEA_HTG <- goSEA_HTG_3_Coldysp

```

## Save the top 10 differentially regulated pathways

```{r}
# class(goSEA_HTG)
# dim(goSEA_HTG)
paste("Col without dysp vs Colitis with dysp, HTG - Number of pathways differentially regulated is", nrow(goSEA_HTG))
paste("Col without dysp vs Colitis with dysp, HTSeq - Number of pathways differentially regulated is", nrow(goSEA_HTseq))
# goSEA_HTG[1:10,1:7]
results_goSEA_HTG <- data.frame(c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"))
results_goSEA_HTG <- t(results_goSEA_HTG)
colnames(results_goSEA_HTG) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue"); results_goSEA_HTG <- results_goSEA_HTG[-1,]
results_goSEA_HTG <- goSEA_HTG[1:10,1:7]
```

### Then, compare the top 10 differentially regulated pathways in the nocol vs colnodysp groups

```{r}
results3 <- data.frame(top10pathways_HTG=results_goSEA_HTG[,2],top10pathways_HTSeq=results_goSEA_HTseq[,2])
results3
write.csv(results3,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/top10_diff_reg_pathways_nocol_coldysp.csv")
setdiff(results_goSEA_HTG[,2],results_goSEA_HTseq[,2])
setdiff(results_goSEA_HTseq[,2],results_goSEA_HTG[,2])
```

### ggplot for visualisation

```{r}

library(stringr)
# df <- goSEA_HTG@result[1:10,] 
df <- results_goSEA_HTG[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTG)") +
  my_theme


df <- results_goSEA_HTseq[1:10,]
head(df)
df$Negative_log10_adjusted_PValue <- -log10(df$p.adjust)
df$Pathway_Description <- str_wrap(df$Description, width = 40)
colnames(df)
ggplot(df, aes(x = Negative_log10_adjusted_PValue, y = reorder(Pathway_Description, Negative_log10_adjusted_PValue), label=ID)) +
  
  geom_bar(stat = "identity")+
  geom_text_repel(
    colour="black",size=4,min.segment.length = 10,
    hjust = 0.5, vjust = 0.5)+
  labs(x = "-Log10(adjusted P-value)", y = "Pathway Description (HTSeq)") +
  my_theme

```

#### Save image

```{r}
save.image("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/GSEA/GSEA.RData")
```

### Attmepting to convert the gene symbols to ENTEREZID first

```{r}

# #Convert symbols to Entrez IDs using bitr (Biological ID translator)
# sigRes$Entrez <- bitr(
#   rownames(sigRes),
#   fromType = "SYMBOL",
#   toType = c("ENTREZID"),
#   OrgDb = org.Mm.eg.db,
#   drop = FALSE)$ENTREZID


gene_names <- data.frame(column1=c("Bcl2a1a_Bcl2a1d","Bcl2a1b","Ccl21_family","Ccl21a","Ccl27_family","Ccl27a","Cd244","Cd244a"))
gene_names <- c("Bcl2a1a_Bcl2a1d","Bcl2a1b","Ccl21_family","Ccl21a","Ccl27_family","Ccl27a","Cd244","Cd244a")

gene_names <- c("Bcl2a1a_Bcl2a1d","Bcl2a1d","Bcl2a1b","Ccl21_family","Ccl21a","Ccl27_family","Ccl27a","Cd244","Cd244a","H2-Ea-ps","H2-Ea","Ifit3_Ifit3b","Ifit3","Ifna_family","Ifna1","Ifna15","Ifna6","Ifnab","Ifnl2_Ifnl3","Ifnl2","Ifnl3","Il22_Iltifb","Il22","Il22b","Klra20","Klra20","Klra21","Klra4_Klra18","Klra7_Klra20","Klra13-ps","Klra4","Klra7","Lyz1_Lyz2","Lyz1","Lyz2","Oas1a_Oas1g","Oas1a","Oas1g","Prame","Pramex1","Skp1a","Skp1")

gene_names$Entrez <- bitr(
  gene_names,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Mm.eg.db,
  drop = FALSE)#$ENTREZID

gene_names2 <- c("BCL2A1A_BCL2A1D","BCL2A1D","BCL2A1B","CCL21_FAMILY","CCL21A","CCL27_FAMILY","CCL27A","CD244","CD244A","H2-EA-PS","H2-EA","IFIT3_IFIT3B","IFIT3","IFNA_FAMILY","IFNA1","IFNA15","IFNA6","IFNAB","IFNL2_IFNL3","IFNL2","IFNL3","IL22_ILTIFB","IL22","IL22B","KLRA20","KLRA20","KLRA21","KLRA4_KLRA18","KLRA7_KLRA20","KLRA13-PS","KLRA4","KLRA7","LYZ1_LYZ2","LYZ1","LYZ2","OAS1A_OAS1G","OAS1A","OAS1G","PRAME","PRAMEX1","SKP1A","SKP1")

gene_names2$Entrez <- bitr(
  gene_names2,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db,
  drop = FALSE)

comparison <- data.frame(gene_names=gene_names$Entrez[,1], EntrezID_Mm=gene_names$Entrez[,2],EntrezID_Hs=gene_names2$Entrez[,2]); comparison
write.csv(comparison,file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/entrezID.csv")

x <- org.Hs.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the entrez gene ID for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}



## then run enrichGO once more
# 
# goSEA_HTseq <- enrichGO(
#   gene = rownames(sigRes_HTSeq),
#   OrgDb = org.Mm.eg.db,
#   keyType = "SYMBOL",
#   ont = "BP", #BP, MF or CC
#   pAdjustMethod = "BH",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05)





# #Convert symbols to Entrez IDs using bitr (Biological ID translator)

trial_sigRes_HTSeq <- sigRes_HTSeq
head(rownames(trial_sigRes_HTSeq))
trial_sigRes_HTSeq$Entrez <- bitr(
  rownames(trial_sigRes_HTSeq),
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Mm.eg.db,
  drop = FALSE)



trial_sigRes_HTG <- sigRes_HTG
head(rownames(trial_sigRes_HTG))
trial_sigRes_HTG$Entrez <- bitr(
  rownames(trial_sigRes_HTG),
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Mm.eg.db,
  drop = FALSE)

# trial_sigRes_HTSeq_entrez <- enrichGO(
#   gene = rownames(trial_sigRes_HTSeq),
#   OrgDb = org.Mm.eg.db,
#   keyType = "ENTREZID",
#   ont = "BP",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   qvalueCutoff = 0.05
#   )

# enrichGO(
#   gene = rownames(sigRes_HTSeq),
#   OrgDb = org.Mm.eg.db,
#   keyType = "ENTREZID",
#   ont = "BP",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH",
#   universe,
#   qvalueCutoff = 0.05,
#   minGSSize = 10,
#   maxGSSize = 500,
#   readable = FALSE,
#   pool = FALSE
# )
```
