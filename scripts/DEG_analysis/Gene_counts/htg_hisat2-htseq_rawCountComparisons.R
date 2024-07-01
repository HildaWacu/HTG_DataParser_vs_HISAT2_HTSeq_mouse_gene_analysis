#1. zero and low expression values filtered
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/htseq_dds_filtering.Rdata")
ddsHTSeq
dds_cpm_filtered_HTSeq 
dds_rowSums_greaterthan_10_HTSeq

###versus

load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/htg_dds_filtering.Rdata")
ddsHTG
dds_cpm_filtered_HTG 
dds_rowSums_greaterthan_10_HTG 


dimensions_of_countTable_HTG <- data.frame(HTG_raw_counts = dim(ddsHTG),
                                       HTG_counts_only_greaterthan10 = dim(dds_rowSums_greaterthan_10_HTG),
                                       HTG_log2_cpm_filtered = dim(dds_cpm_filtered_HTG)) ; knitr::kable(dimensions_of_countTable_HTG)

dimensions_of_countTable_HTSeq <- data.frame(HTSeq_raw_counts = dim(ddsHTSeq),
                                           HTSeq_counts_only_greaterthan10 = dim(dds_rowSums_greaterthan_10_HTSeq),
                                           HTSeq_log2_cpm_filtered = dim(dds_cpm_filtered_HTSeq)) ; knitr::kable(dimensions_of_countTable_HTSeq)

#2. library size following HISAT2 alignment vs HTG data parser software

# Calculate library size for each sample in the HTG counts
##before normalisation
library_sizes_HTGraw <- colSums(counts(ddsHTG))

#after log2(cpm) less than 1 filtering
library_sizes_HTGcpm <- colSums(counts(dds_cpm_filtered_HTG))

#after greater-than-10 filtering --> similar to raw counts beacause none were filtered
library_sizes_HTGgreaterthan10counts <- colSums(counts(dds_rowSums_greaterthan_10_HTG))

### versus
# Calculate library size for each sample in the HISAT2-HTSeq counts
##before normalisation
library_sizes_HTSEQraw <- colSums(counts(ddsHTSeq))

#after log2(cpm) less than 1 filtering
library_sizes_HTSEQcpm <- colSums(counts(dds_cpm_filtered_HTSeq))

#after greater-than-10 filtering --> similar to raw counts beacause none were filtered
library_sizes_HTSEQgreaterthan10counts <- colSums(counts(dds_rowSums_greaterthan_10_HTSeq))


### plot the output

boxplot(library_sizes_HTGraw,library_sizes_HTSEQraw,library_sizes_HTGgreaterthan10counts,library_sizes_HTSEQgreaterthan10counts,library_sizes_HTGcpm,library_sizes_HTSEQcpm)

# Define the datasets and their labels

# Create a data frame with the datasets
datasets <- data.frame(
  HTGraw = library_sizes_HTGraw,
  HTSEQraw = library_sizes_HTSEQraw,
  HTGgreaterthan10counts = library_sizes_HTGgreaterthan10counts,
  HTSEQgreaterthan10counts = library_sizes_HTSEQgreaterthan10counts,
  HTGcpm = library_sizes_HTGcpm,
  HTSEQcpm = library_sizes_HTSEQcpm
)

# Define the datasets and their labels
# datasets_as_list <- list(
#   library_sizes_HTGraw,
#   library_sizes_HTSEQraw,
#   library_sizes_HTGgreaterthan10counts,
#   library_sizes_HTSEQgreaterthan10counts,
#   library_sizes_HTGcpm,
#   library_sizes_HTSEQcpm
# )
# Define the labels
# labels <- c(
#   "HTG Raw",
#   "HTSeq Raw",
#   "HTG >10 Counts",
#   "HTSeq >10 Counts",
#   "HTG CPM",
#   "HTSeq CPM"
# )

labels <- c(
  "HTG",
  "HTSeq",
  "HTG",
  "HTSeq",
  "HTG",
  "HTSeq"
)


# Define the colors for the boxes
colors <- c("red", "#FF5500", "darkorange", "orange", "green", "lightgreen")

# Create the box plot with labeled and colored boxes
boxplot(datasets, names = labels, col = colors)
legend("right", legend = labels, fill = colors)

