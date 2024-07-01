## load the sample 1 data
sample1_gene_counts <- read.delim("D:/trials/sample1_gene_counts.txt", header=FALSE)
View(sample1_gene_counts)
tail(sample1_gene_counts)

#subset to get only the gene features
sample_1 <- sample1_gene_counts[1:40879,]  ## this is hard coded. Think of a way to filter it out using the starting double underscores characteristic
head(sample_1)
tail(sample_1)


sum(sample1_gene_counts[,2])
sum(sample_1[,2])
