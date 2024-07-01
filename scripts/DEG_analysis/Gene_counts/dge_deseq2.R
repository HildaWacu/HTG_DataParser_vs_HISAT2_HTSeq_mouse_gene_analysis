#Step 1: Define design matrix
designMatrix <- model.matrix(~ condition, data = sampleTable)
designMatrix[1:5, 1:4]


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





