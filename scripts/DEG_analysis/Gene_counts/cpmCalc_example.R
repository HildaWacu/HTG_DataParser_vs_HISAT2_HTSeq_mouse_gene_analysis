# Example raw counts matrix
raw_counts <- matrix(c(10, 20, 30, 40, 50, 60), nrow = 2, ncol = 3)

# Calculate library size for each sample
library_sizes <- colSums(raw_counts)

# Calculate CPM
cpm <- t(t(raw_counts) / library_sizes) * 1e6
