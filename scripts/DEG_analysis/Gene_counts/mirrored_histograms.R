##### To plot mirrored histograms of the meanLogCPM() values from the two pipelines

##HTSeQ
### filtered using the 0.1 threshold and not 0.05

# hist(assay(normCounts_vst), main = "HISAT2-HTSeq Pipeline output when VST normalised",xlim = c(0,20), breaks = seq(0, 20, by = 0.8), ylim = c(0,15000))

meanLog2CPM <- rowMeans(log2(cpm(counts(ddsHTSeq)) + 1)); htseq_meanLog2CPM <- meanLog2CPM

paste(sum(meanLog2CPM <= 0.1), "genes are filtered")
sum(meanLog2CPM > 0.1)
cpm_filtered_htseq <- counts(ddsHTSeq[meanLog2CPM > 0.1,])
cpm_filtered_htseq_meanLog2CPM <- rowMeans(log2(cpm(cpm_filtered_htseq) + 1))

## HTG

meanLog2CPM <- rowMeans(log2(cpm(counts(ddsHTG)) + 1)) ; htg_meanLog2CPM <- meanLog2CPM

paste(sum(meanLog2CPM <= 0.1), "genes are filtered")
sum(meanLog2CPM > 0.1)
cpm_filtered_htg <- counts(ddsHTG[meanLog2CPM > 0.1,])
cpm_filtered_htg_meanLog2CPM <- rowMeans(log2(cpm(cpm_filtered_htg) + 1))


### before filtering

par(mfrow=c(2,1))

#Make the plot
par(mar=c(0,5,3,3))
hist(htseq_meanLog2CPM, main = "Mirrored histograms showing the distribution of gene expression\
     values (in mean-log2-CPM) for all genes, across all samples",
     ylab="Frequency", 
     xlab="",
     xlim = c(0,20), 
     breaks = seq(0, 20, by = 0.1), 
     xaxt="n", 
     las=1 , 
     col="tomato3")

expr_cutoff=0.1
abline(v = expr_cutoff, col = "orange", lwd = 1)
# hist(x1 , main="" , xlim=c(-2,5), ylab="Frequency for x1", xlab="", ylim=c(0,25) , xaxt="n", las=1 , col="slateblue1", breaks=10)


par(mar=c(5,5,0,3))
hist(htg_meanLog2CPM, main = "",
     xlim = c(0,20), 
     breaks = seq(0, 20, by = 0.1), 
     ylim = c(40,0), 
     ylab="Frequency", 
     xlab="Gene expression values (mean-Log2-CPM)",
     las=1 , 
     col="slateblue1" )

expr_cutoff=0.1
abline(v = expr_cutoff, col = "orange", lwd = 1)
# hist(x2 , main="" , xlim=c(-2,5), ylab="Frequency for x2", xlab="Value of my variable", ylim=c(25,0) , las=1 , col="tomato3"  , breaks=10)



##After filtering

par(mfrow=c(2,1))

#Make the plot
par(mar=c(0,5,3,3))
hist(cpm_filtered_htseq_meanLog2CPM, main = "",
     ylab="Frequency", 
     xlim = c(0,20), 
     breaks = seq(0, 20, by = 0.1), 
     ylim = c(0,40), 
     xaxt="n", 
     las=1 , 
     col="tomato3")

expr_cutoff=0.1
abline(v = expr_cutoff, col = "orange", lwd = 1)
# hist(x1 , main="" , xlim=c(-2,5), ylab="Frequency for x1", xlab="", ylim=c(0,25) , xaxt="n", las=1 , col="slateblue1", breaks=10)


par(mar=c(5,5,0,3))
hist(cpm_filtered_htg_meanLog2CPM, main = "",
     xlim = c(0,20), 
     breaks = seq(0, 20, by = 0.1), 
     ylim = c(40,0), 
     ylab="Frequency", 
     xlab="Gene expression values (mean-Log2-CPM)",
     las=2 , 
     col="slateblue1" )

expr_cutoff=0.1
abline(v = expr_cutoff, col = "orange", lwd = 1)
# hist(x2 , main="" , xlim=c(-2,5), ylab="Frequency for x2", xlab="Value of my variable", ylim=c(25,0) , las=1 , col="tomato3"  , breaks=10)



################################################
## TRYING TO SPLIT THE X AXIS

# Import required libraries  
install.packages("plotrix")
library("plotrix") 

# Create data 
x<-50:1 

# Apply gap.plot function 
gap.barplot(x,gap=c(2,4))
gap.barplot((hist(cpm_filtered_htg_meanLog2CPM, main = "",
              xlim = c(0,20), 
              breaks = seq(0, 20, by = 0.1),
              ylim = c(0,100),
              ylab="Frequency", 
              xlab="Gene expression values (mean-Log2-CPM)",
              las=1 , 
              col="slateblue1" )),gap = c(10,20))




