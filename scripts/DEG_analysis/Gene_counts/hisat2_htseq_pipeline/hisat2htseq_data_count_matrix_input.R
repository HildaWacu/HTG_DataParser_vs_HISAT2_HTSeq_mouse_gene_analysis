# Prepare the DESeqDataSet object for the HISAT2-HTSeq Pipeline


## Load the count_matrix_input.Rdata image to get the coldata

load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/count_matrix_input.Rdata")

## Prepeare a .csv file with the right order of samples. # 
sampleFiles <- data.frame(x = rownames(coldata),
                          y = coldata$condition2)
write.csv(sampleFiles, file="C:/hne_files/MSc_Bioinf/hisat2htseq_dataset/sample_list.csv",quote = FALSE,row.names = FALSE)

##Using a bash script reorder_files.sh, reorder these count files by the order in the .csv sampleFile
## sort these reorderedHTSeq_files according to their timestamp

setwd("C:/hne_files/MSc_Bioinf/hisat2htseq_dataset/reorderedHTSeq_files/")
directory <- "C:/hne_files/MSc_Bioinf/hisat2htseq_dataset/reorderedHTSeq_files/"
list.files()
details <-  file.info(list.files(pattern="*.txt")) ; head(details)

details <-  details[with(details, order(as.POSIXct(mtime))), ]
reodered_sampleFiles <-  rownames(details)
head(sampleFiles) ;head(reodered_sampleFiles)

##Prepare the DESeqDataSet from the ordered HTSeq count files
sampleTable <- data.frame(sampleName = reodered_sampleFiles,
                          fileName = reodered_sampleFiles,
                          condition = coldata$condition2) ;head(sampleTable)
sampleTable$condition <- factor(sampleTable$condition)


## build then we build the DESeqDataSet object

library("DESeq2")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq
dds <- ddsHTSeq

### ---- save results----
save.image(file="C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/hisat2htseq_deseqdataset.Rdata")
