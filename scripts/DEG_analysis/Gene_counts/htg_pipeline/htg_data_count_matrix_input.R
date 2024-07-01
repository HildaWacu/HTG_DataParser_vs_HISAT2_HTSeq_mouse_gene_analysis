#===========================================================================================================================================================
#---- READING FROM THE SECOND SHEET; because it is "absolutely critical that the columns of the count matrix and the rows of the 
#column data (information about samples) are in the same order"

# ## ----installing and loading the necessary packages----
# library(BiocManager)
# if (!require("openxlsx", quietly = TRUE))
#   install.packages("openxlsx")
# #BiocManager::install()
# 
# library(openxlsx)
# 
# ##---- set the working directory----
# setwd("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/")
# 
# 
# ## Load data set --> sheet 1
# path_to_dataset <- "C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH230615.xlsx" #different file-->most recent, check date
# htg_raw_counts_all_sh1 <- read.xlsx(path_to_dataset) 
# 
# # Select the gene features, samples, metadata and controls, and have the first row as sampleID
# htg_raw_counts_with_metainfo_controls_sh1 <- read.xlsx(path_to_dataset, startRow = 8, colNames = TRUE) 
# 
# head(htg_raw_counts_with_metainfo_controls_sh1, 6)
# htg_raw_counts_sh1 <- htg_raw_counts_with_metainfo_controls_sh1[12:nrow(htg_raw_counts_with_metainfo_controls_sh1),] ; head(htg_raw_counts_sh1,3)
# 
# htg_metainfo_controls <- htg_raw_counts_all_sh1[1:17,]; htg_metainfo_controls[,1:7]
# 
# 
# ## assign gene features as the row names
# htg_gene_counts_sh1 <- htg_raw_counts_sh1[,-1]
# 
# row.names(htg_gene_counts_sh1) <- htg_raw_counts_sh1[,1]; head(row.names(htg_gene_counts_sh1))
# 
# dim(htg_gene_counts_sh1)
# 
# 
# # defining the coldata columns one by one (SHOULD INSTEAD HAVE A CSV FILE DEFINING THIS INSTEAD OF DOING IT MANUALLY)
# # read form sheet two as defined by Bettan
# htg_raw_counts_all_sh2 <- read.xlsx(path_to_dataset,"Groupwise - Bettan")
# sampleID_row <- htg_raw_counts_all_sh2[6,]
# 
# sampleID_row_trans <- t(sampleID_row[,-1])
# # rownames(sampleID_row_trans) <- paste(1:nrow(sampleID_row_trans))
# # colnames(sampleID_row_trans) <- "sampleID"
# # coldata<-sampleID_row_trans
# coldata_col1 <- data.frame(sampleID_row_trans,row.names = paste(1:nrow(sampleID_row_trans)))
# colnames(coldata_col1) <- "sampleID"
# coldata_col2 <- data.frame(coldata_col1, condition1=c(rep("colitis",42),rep("healthy",22)))
# coldata_col3 <- data.frame(coldata_col2, condition2=c(rep("colitis&dysplasia",10), rep("colitis&nodysplasia",10), rep("colitisOnly",22),rep("noColitis",22)))
# 
# coldata_col4 <- data.frame(coldata_col3, condition3=c(rep("dysplasia",10), c(rep("nodyplasia",54))))
# coldata_col5 <- data.frame(coldata_col4, condition4=c(rep("paired",20), c(rep("notpaired",44))))
# 
# ## have the sampleID as the  column names 
# # coldata <- coldata_col3[,-1]
# # row.names(coldata) <- coldata_col3[,1]
# 
# coldata <- coldata_col5[,-1]
# row.names(coldata) <- coldata_col5[,1]
# #===========================================================================================================================================================
# ##### check if the naming between counts columns(sampleID) and rows of coldata is consistent
# 
# all(rownames(coldata) %in% colnames(htg_gene_counts_sh1))  ## minus 1 because the first column name is for the features
# 
# if (all(rownames(coldata) %in% colnames(htg_gene_counts_sh1))) {
#   print("same column and row names are present in count and coldata files")
# } else {stop("column and row names present")
# }
# 
# ##### ---- to have the coldata rownames matching exactly with count data colnames
# head(coldata)
# rownames(coldata)
# 
# if (all(rownames(coldata) == colnames(htg_gene_counts_sh1))) {
#   print("column and row names match in count and coldata files")
# } else {stop("column and row names do not match")
# }
# ###if false, run the following command
# htg_gene_counts_sh1 <- htg_gene_counts_sh1[, rownames(coldata)] ## this command orders the columns(sampleID) in htg_counts by the order in the coldata rownames
# 
# if (all(rownames(coldata) == colnames(htg_gene_counts_sh1))) {
#   print("column and row names match in count and coldata files")
# } else {stop("column and row names do not match")
# }
# 
# 
# #####=====================================================================
# #### 
# 
# ## ----installing and loading the necessary packages----
# library(BiocManager)
# if (!require("openxlsx", quietly = TRUE))
#   install.packages("openxlsx")
# #BiocManager::install()
# 
# library(openxlsx)
# # Load the data set --> sheet 2, with samples ordered by disease state
# # path_to_dataset <- "C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH210622.xlsx"
# path_to_dataset <- "C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH230615.xlsx" #different file-->most recent, check date
# htg_raw_counts_all <- read.xlsx(path_to_dataset,"Groupwise - Bettan")
#  
# 
# htg_metainfo_controls <- htg_raw_counts_all[1:17,]; htg_metainfo_controls[,1:7]
# 
# # Select the gene features, samples, metadata and controls, and also label the colnames by sampleID 
# htg_raw_counts_with_metainfo_controls <- read.xlsx(path_to_dataset,"Groupwise - Bettan", startRow = 8, colNames = TRUE) 
# head(htg_raw_counts_with_metainfo_controls, 6)
# htg_raw_counts <- htg_raw_counts_with_metainfo_controls[12:nrow(htg_raw_counts_with_metainfo_controls),] ; head(htg_raw_counts,3)
# 
# 
# ## assign gene features as the row names
# htg_gene_counts <- htg_raw_counts[,-1]
# row.names(htg_gene_counts) <- htg_raw_counts[,1]; head(row.names(htg_gene_counts))
# 
# dim(htg_gene_counts)
# 
# 
# #===========================================================================================================================================================
# 
# # defining the coldata columns one by one
# sampleID_row <- htg_raw_counts_all[6,]
# 
# sampleID_row_trans <- t(sampleID_row[,-1])
# # rownames(sampleID_row_trans) <- paste(1:nrow(sampleID_row_trans))
# # colnames(sampleID_row_trans) <- "sampleID"
# # coldata<-sampleID_row_trans
# coldata_col1 <- data.frame(sampleID_row_trans,row.names = paste(1:nrow(sampleID_row_trans)))
# colnames(coldata_col1) <- "sampleID"
# coldata_col2 <- data.frame(coldata_col1, condition1=c(rep("colitis",42),rep("healthy",22)))
# coldata_col3 <- data.frame(coldata_col2, condition2=c(rep("colitis&dysplasia",10), rep("colitis&nodysplasia",10), rep("colitisOnly",22),rep("noColitis",22)))
# 
# coldata_col4 <- data.frame(coldata_col3, condition3=c(rep("dysplasia",10), c(rep("nodyplasia",54))))
# coldata_col5 <- data.frame(coldata_col4, condition4=c(rep("paired",20), c(rep("notpaired",44))))
# 
# ## have the sampleID as the  column names 
# # coldata <- coldata_col3[,-1]
# # row.names(coldata) <- coldata_col3[,1]
# 
# coldata <- coldata_col5[,-1]
# row.names(coldata) <- coldata_col5[,1]
# #===========================================================================================================================================================
# ##### check if the naming between counts and coldata_col3 is consistent
# 
# all(rownames(coldata) %in% colnames(htg_gene_counts))  ## minus 1 because the first column name is for the features
# 
# if (all(rownames(coldata) %in% colnames(htg_gene_counts))) {
#     print("same column and row names are present in count and coldata files")
#   } else {stop("column and row names present")
#     }
# 
# ##### ---- to have the coldata rownames matching exactly with count data colnames
# head(coldata)
# rownames(coldata)
# 
# if (all(rownames(coldata) == colnames(htg_gene_counts))) {
#   print("column and row names match in count and coldata files")
# } else {stop("column and row names do not match")
# }
# ###if false, run the following command
# htg_gene_counts <- htg_gene_counts[, rownames(coldata)] ## this command orders the columns(sampleID) in htg_counts by the order in the coldata rownames
# 
# if (all(rownames(coldata) == colnames(htg_gene_counts))) {
#   print("column and row names match in count and coldata files")
# } else {stop("column and row names do not match")
# }
# 
# ### finally, convert the conditions column contents to factors 
# coldata$condition1 <- factor(coldata$condition1)
# coldata$condition2 <- factor(coldata$condition2)
# coldata$condition3 <- factor(coldata$condition3)
# coldata$condition4 <- factor(coldata$condition4)


#####=====================================================================
####
###TRY A DIFF. APPROACH

##---- set the working directory----
setwd("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/")


## Load data set --> sheet 1
path_to_dataset <- "C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH230615.xlsx" #different file-->most recent, check date
library(readxl)
column_names <- paste("sample_", sep = "", seq(1:64))

Appendix_1_Raw_mRNA_counts_EHH230615 <- 
  read_excel("C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH230615.xlsx")

htg_gene_counts <- 
  read_excel("C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH230615.xlsx", 
                                                   skip = 19, col_names = FALSE) ## skipping the first 19 rows containing metadata

write.csv(htg_gene_counts,file = "htg_gene_counts", row.names = FALSE)

htg_countTable <- read.table("htg_gene_counts",
                         as.is = TRUE, header = TRUE, sep = ",", row.names = 1)

dim(htg_countTable)
colnames(htg_countTable)
colnames(htg_countTable) <- column_names ;head(htg_countTable)


### defining the coldata columns one by one

# get the order of samples from sheet 2 of the data excel file - as defined by Bettan

sampleID_row <- 
  read_excel("C:/hne_files/MSc_Bioinf/HTG_dataset/raw_gene_count_data/Appendix_1_Raw_mRNA_counts_EHH230615.xlsx", 
             sheet = 2, 
             range = "A8:BM8",col_names = FALSE)  ##range of cells containing the sampleID order
head(sampleID_row)

## begin to create the coldata (i.e., sampleTable) --> condition2 most important classification
sampleID_row_trans <- t(sampleID_row[,-1])
coldata_col1 <- data.frame(sampleID_row_trans,row.names = paste(1:nrow(sampleID_row_trans)))
colnames(coldata_col1) <- "sampleID"
coldata_col2 <- data.frame(coldata_col1, condition1=c(rep("colitis",42),rep("healthy",22)))
coldata_col3 <- data.frame(coldata_col2, condition2=c(rep("colitis_dysplasia",10), 
                                                      rep("colitis_nodysplasia",10), 
                                                      rep("colitisOnly",22),
                                                      rep("noColitis",22)))

coldata_col4 <- data.frame(coldata_col3, condition3=c(rep("dysplasia",10), c(rep("nodyplasia",54))))
coldata_col5 <- data.frame(coldata_col4, condition4=c(rep("paired",20), c(rep("notpaired",44))))


# to rename the rownames in coldata to sample_#
sampleID_labels <- paste("sample_",  sep = "", coldata_col5[,1])
## have the sampleID as the  column names 
# coldata <- coldata_col3[,-1]
# row.names(coldata) <- coldata_col3[,1]
coldata <- coldata_col5[,-1]
row.names(coldata) <- sampleID_labels ; head(coldata)

#===========================================================================================================================================================
##### check if the naming between counts and coldata_col3 is consistent

#all(rownames(coldata) %in% colnames(htg_countTable))  ## minus 1 because the first column name is for the features

if (all(rownames(coldata) %in% colnames(htg_countTable))) {
  print("same column and row names are present in count and coldata files")
} else {stop("column and row names present")
}

##### ---- to have the coldata rownames matching exactly with count data colnames
head(coldata)
rownames(coldata)
colnames(htg_countTable)

if (all(rownames(coldata) == colnames(htg_countTable))) {
  print("column and row names match in count and coldata files")
} else {stop("column and row names do not match")
}
###if false, run the following command
htg_countTable <- htg_countTable[, rownames(coldata)] ## this command orders the columns(sampleID) in htg_counts by the order in the coldata rownames

if (all(rownames(coldata) == colnames(htg_countTable))) {
  print("column and row names match in count and coldata files")
} else {stop("column and row names do not match")
}

### finally, convert the conditions column contents to factors 
coldata$condition1 <- factor(coldata$condition1)
coldata$condition2 <- factor(coldata$condition2)
coldata$condition3 <- factor(coldata$condition3)
coldata$condition4 <- factor(coldata$condition4)



### ---- save results----
save.image(file="count_matrix_input.Rdata")


### all(rownames(coldata_col3[,1]) == colnames(htg_gene_counts)) ## not sure why " 'NULL' is deemed equal to values if columns in htg_gene_count"

#===========================================================================================================================================================
#===========================================================================================================================================================
#===========================================================================================================================================================
#===========================================================================================================================================================



## https://sparkbyexamples.com/r-programming/select-rows-in-r/