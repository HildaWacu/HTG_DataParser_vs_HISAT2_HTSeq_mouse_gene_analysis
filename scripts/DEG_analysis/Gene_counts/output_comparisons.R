#### Comparing log transformed and VST transformed data to stabilise the variace across mean values 


## Load R image to get the filtered data

load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/filtering_criteria.Rdata")

dds <- dds_cpm_filtered_HTSeq1

par( mfrow = c( 2, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.7, main="Log2 of size-factor normalised gene counts")
abline(a = 0, b = 1, col = "red")
# plot(assay(rld)[,1:2], pch=16, cex=0.7, main=("r-log normalised gene counts"))
# abline(a = 0, b = 1, col = "red")
plot(assay(normCounts_vst_htseq[,1:2]),pch=16, cex=0.7, main="VST normalised gene counts")
abline(a = 0, b = 1, col = "red")



dds <- dds_cpm_filtered_HTG1

# par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.7, main="Log2 of size-factor normalised gene counts")
abline(a = 0, b = 1, col = "red")
# plot(assay(rld)[,1:2], pch=16, cex=0.7, main=("r-log normalised gene counts"))
# abline(a = 0, b = 1, col = "red")
plot(assay(normCounts_vst_htg[,1:2]),pch=16, cex=0.7, main="VST normalised gene counts")
abline(a = 0, b = 1, col = "red")



#### Comparing the qqplots 
par(mfrow = c(1, 2))
qqnorm(assay(normCounts_vst_htseq), main = "Normal Q-Q Plot for the HISAT2-HTSeq count data")
qqline(assay(normCounts_vst_htseq), col = "red")

qqnorm(assay(normCounts_vst_htg), main = "Normal Q-Q Plot for the HTG count data")
qqline(assay(normCounts_vst_htg), col = "red")


### comparing the histograms of the distribution
par(mfrow = c(1,2))
hist(assay(normCounts_vst_htseq), main = "HTSeq Count data", xlim = c(0,20), xlab = "VST normalised gene expression values")
hist(assay(normCounts_vst_htg), main = "HTG count data", xlim = c(0,20), xlab = "VST normalised gene expression values")



hist(assay(normCounts_vst_htseq), main = "HTSeq Count data", xlim = c(4,20), xlab = "VST normalised gene expression values")
hist(assay(normCounts_vst_htg), main = "HTG count data", xlim = c(4,20), xlab = "VST normalised gene expression values")




pipeline_rawcounts_pre_post_filtering <- data.frame(HTSeq_raw_counts = dim(ddsHTSeq),
                                             HTSeq_filtered_raw_counts = dim(dds_cpm_filtered_HTSeq1),
                                             HTG_raw_counts = dim(ddsHTG),
                                             HTG_filtered_raw_counts = dim(dds_cpm_filtered_HTG1)) ; knitr::kable(pipeline_rawcounts_pre_post_filtering)



rownames(pipeline_rawcounts_pre_post_filtering) <- c("gene_features", "samples") ; pipeline_rawcounts_pre_post_filtering
write.csv(pipeline_rawcounts_pre_post_filtering, file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/pipeline_rawcounts_pre_post_filtering.csv")


### To comapre the gene names present in both pipelines in the raw count data

#### The Jaccard package

##### Step 1: Install and load package
```{r}
if (!require(devtools)) {install.packages("devtools")}
# library("devtools")
if (!require("jaccard")) {devtools::install_github("ncchung/jaccard")}
```

##### Step 2: Prepare the data sets

I. Load the necessary R image to acquire the genelists : in this case, the filtered raw counts
```{r}

load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/filtering_criteria.Rdata")
```


### prepare dataframes with the genenames from the raw gene count matrix from the two pipelines
genename_htseq <- data.frame("htseq" = (rownames(counts(dds_cpm_filtered_HTSeq1))))
head(genename_htseq)

genename_htg <- data.frame("htg" = (rownames(counts(dds_cpm_filtered_HTG1))))
head(genename_htg)

###### noCol vs colOnly groups
II. Prepare the data.frame with the genes present and absent indicated by 1s and 0s, respectively

```{r}
# create an empty data frame
jaccard_list_raw_genenames <- data.frame()

for (i in 1:nrow(genename_htseq)) {
  # print(i)
  new_row <- data.frame(htg = 0, htseq = 1)
  jaccard_list_raw_genenames <- rbind(jaccard_list_raw_genenames, new_row)
}
head(jaccard_list_raw_genenames)

## to confirm that the number of rows is correct
nrow(jaccard_list_raw_genenames) ; length(genename_htseq$htseq)


#### to assign the rownames in the jaccard_list.... the genenames in genename_htseq df
rownames(jaccard_list_raw_genenames) <- genename_htseq$htseq

head(jaccard_list_raw_genenames)

### Do the same for the HTG raw gene names
jaccard_list_raw_genenames2 <- data.frame()

for (i in 1:nrow(genename_htg)) {
  new_row2 <- data.frame(htg=1, htseq=0)
  jaccard_list_raw_genenames2 <- rbind(jaccard_list_raw_genenames2,new_row2)
}
rownames(jaccard_list_raw_genenames2) <- genename_htg$htg

jaccard_list_raw_genenames <- rbind(jaccard_list_raw_genenames,jaccard_list_raw_genenames2)

## to confirm that the number of rows created is correct
nrow(jaccard_list_raw_genenames) ; sum(nrow(genename_htg), nrow(genename_htseq))

head(jaccard_list_raw_genenames)
```

## to identify the unique gene names in the raw gene coutn amtrix from the two pipelines
in_htg_rawgenes <- setdiff(genename_htg$htg,genename_htseq$htseq) ; head(in_htg_rawgenes) ; length(in_htg_rawgenes)
in_htseq_rawgenes <- setdiff(genename_htseq$htseq,genename_htg$htg) ; head(in_htseq_rawgenes) ; length(in_htseq_rawgenes)

#does not work since the row number is different between the columns being created is different
unique_genes <- data.frame("htseq" = (setdiff(genename_htseq$htseq,genename_htg$htg)), "htg" = (setdiff(genename_htg$htg,genename_htseq$htseq)))

unique_genes_htseq <- data.frame("htseq" = (setdiff(genename_htseq$htseq,genename_htg$htg)))
unique_genes_htg <- data.frame("htg" = (setdiff(genename_htg$htg,genename_htseq$htseq)))
# unique_rows_rbound <- data.frame("htseq" = (setdiff(genename_htseq$htseq,genename_htg$htg)))

unique_rows_rbound <- data.frame(unique_genes_htseq) ; colnames(unique_rows_rbound) <- "gene_names"
unique_rows_rbound2 <- data.frame(unique_genes_htg) ;colnames(unique_rows_rbound2) <- "gene_names"

unique_rows_rbound <- rbind(unique_rows_rbound, unique_rows_rbound2) 
write.csv(unique_rows_rbound, file = "C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/results/pipeline_comparisons/unique_gene_names_rbound.csv")


```{r}
# from the Jaccard list, derive the non-unique
non_unique_rows <- rownames(jaccard_list_raw_genenames)[duplicated(c(genename_htseq$htseq,genename_htg$htg))]

```


```{r}

# create an empty data frame
jaccard_list_raw_genenames <- data.frame()

for (i in 1:length(non_unique_rows)) {
  # print(i)
  new_row <- data.frame(htg = 1, htseq = 1)
  jaccard_list_raw_genenames <- rbind(jaccard_list_raw_genenames, new_row)
}
head(jaccard_list_raw_genenames)

## to confirm that the number of rows is correct
nrow(jaccard_list_raw_genenames) ; length(non_unique_rows)


rownames(jaccard_list_raw_genenames) <- non_unique_rows
head(jaccard_list_raw_genenames)


# unique in HTG pipeline
jaccard_list_raw_genenames2 <- data.frame()
in_htg <- setdiff(genename_htg$htg,genename_htseq$htseq)
for (i in 1:length(in_htg)) {
  new_row2 <- data.frame(htg=1, htseq=0)
  jaccard_list_raw_genenames2 <- rbind(jaccard_list_raw_genenames2,new_row2)
}
rownames(jaccard_list_raw_genenames2) <- in_htg

jaccard_list_raw_genenames <- rbind(jaccard_list_raw_genenames,jaccard_list_raw_genenames2)


## To confirm that the number of rows created is correct, i.e., To confirm that the concatenation has been done accurately = 
nrow(jaccard_list_raw_genenames) ; sum(length(in_htg), length(non_unique_rows))


# unique in HTSeq pipeline
jaccard_list_raw_genenames3 <- data.frame()
in_htseq <- setdiff(genename_htseq$htseq,genename_htg$htg)
for (i in 1:length(in_htseq)) {
  new_row3 <- data.frame(htg=0, htseq=1)
  jaccard_list_raw_genenames3 <- rbind(jaccard_list_raw_genenames3,new_row3)
}
rownames(jaccard_list_raw_genenames3) <- in_htseq

jaccard_list_raw_genenames <- rbind(jaccard_list_raw_genenames,jaccard_list_raw_genenames3)

## to confirm that the number of rows created is correct = 1715
nrow(jaccard_list_raw_genenames) ; sum(length(in_htg),length(in_htseq),length(non_unique_rows))

head(jaccard_list_raw_genenames)
```


III. Run the Jaccard tool - first the Jaccard coeff
```{r}
#first create thethe binary vectors as input for jaccard
binary1 <- jaccard_list_raw_genenames$htg
binary2 <- jaccard_list_raw_genenames$htseq

jaccard::jaccard(binary1,binary2); jaccard::jaccard(binary1,binary2)*100
```
jaccard.test(binary1,binary2, method = "exact") ## to test the significance of the Jaccard Similarity Coeff






              