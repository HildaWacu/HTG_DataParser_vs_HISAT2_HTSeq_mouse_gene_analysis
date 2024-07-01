library(edgeR)  # To access the cpm function

# A function that takes a matrix and decides which rows (genes) should be kept, based on filtering criteria
filter_function <- function(matrix, group_col_inds_list, perc, threshold)
  {
  # matrix: A data matrix
  # group_col_indices_list: a list where each element is a vector of column indices, for example
  # perc: how many percent of the samples in a group need to have a gene expressed above the threshold value?
  # thresh: a cutoff value below which a gene is considered too lowly expresssed
  
  keep_rows <- c()  # Vector (as yet empty) to hold numbers of the rows that should be kept
  
  # Loop over rows (genes)
  for(i in 1:nrow(matrix))
    {
    keep_row <- 0  # Default is 0 - to not keep the row (gene)
    gene <- as.numeric(matrix[i,])
    
    # Loop over groups of samples
    for(j in 1:length(group_col_inds_list))
      {
      group_inds <- group_col_inds_list[[j]]
      sample_values <- gene[group_inds]  
      
      # How many samples in this group are expressed above the threshold?
      num_above <- length(which(sample_values>=threshold))
      
      # How many percent is that of all samples in the group?
      perc_above <- 100*(num_above/length(sample_values))
      
      # Does this qualify the current gene for being kept in the dataset?
      if(perc_above>=perc)  { keep_row <- 1 }
      }
      
    # If we have now concluded that the current gene should be kept, add its row number to the vector of
    # rows to keep
    if(keep_row==1) { keep_rows <- c(keep_rows, i) }
    }
  
  # All rows done! Return the resulting vector of row numbers to keep
  return(keep_rows)  
  } # end of filter_function



# Make example data
matrix1 <- matrix(sample(1:10000, 5000, replace=T), nrow=100, ncol=50)
matrix1_cpm <- cpm(matrix1)

# Which columns indices represent groups of samples?
gr1 <- 1:20
gr2 <- 21:40
gr3 <- 41:50
group_list <- list(gr1, gr2, gr3)

keep_rows <- filter_function(matrix=matrix1_cpm, group_col_inds_list=group_list, perc=90, threshold=3000)

# Go back to the original matrix (not cpm, but gene counts), and subset it using keep_rows
matrix1_filtered <- matrix1[keep_rows,]


