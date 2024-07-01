## to plot bargraphs of up- and down-regulated DEG

## Load the R.images
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/htg_pipeline/deseq2_statistical_analysis_HTG.Rdata")
load("C:/hne_files/MSc_Bioinf/MScProject-main_C/MScProject-main_10.06.23/scripts/DEG_analysis/Gene_counts/hisat2_htseq_pipeline/deseq2_statistical_analysis_HTSeq.Rdata")


up <- nrow(subset(sigRes_HTG_reg,sigRes_HTG_reg$regulation == "up"))
down <- nrow(subset(sigRes_HTG_reg,sigRes_HTG_reg$regulation == "down"))

col_only_reg_df <- data.frame(htg_colonly = c(up,down))

up <- nrow(subset(sigRes_HTSeq_reg,sigRes_HTSeq_reg$regulation == "up"))
down <- nrow(subset(sigRes_HTSeq_reg,sigRes_HTSeq_reg$regulation == "down"))

col_only_reg_df$htseq_colonly <- c(up,down)

# col_only_reg_df$reg <- c("up","down")

#reorder the columns

#Col no Dysp

up <- nrow(subset(sigRes_HTG_Colnodysp_reg,sigRes_HTG_Colnodysp_reg$regulation == "up"))
down <- nrow(subset(sigRes_HTG_Colnodysp_reg,sigRes_HTG_Colnodysp_reg$regulation == "down"))

col_nodysp_reg_df <- data.frame(htg_col_nodysp = c(up,down))

up <- nrow(subset(sigRes_HTSeq_Colnodysp_reg,sigRes_HTSeq_Colnodysp_reg$regulation == "up"))
down <- nrow(subset(sigRes_HTSeq_Colnodysp_reg,sigRes_HTSeq_Colnodysp_reg$regulation == "down"))

col_nodysp_reg_df$htseq_col_nodysp <- c(up,down)

## Col Dysp

up <- nrow(subset(sigRes_HTG_Coldysp_reg,sigRes_HTG_Coldysp_reg$regulation == "up"))
down <- nrow(subset(sigRes_HTG_Coldysp_reg,sigRes_HTG_Coldysp_reg$regulation == "down"))

col_dysp_reg_df <- data.frame(htg_col_dysp = c(up,down))

up <- nrow(subset(sigRes_HTG_Coldysp_reg,sigRes_HTG_Coldysp_reg$regulation == "up"))
down <- nrow(subset(sigRes_HTG_Coldysp_reg,sigRes_HTG_Coldysp_reg$regulation == "down"))

col_dysp_reg_df$htseq_col_dysp <- c(up,down)


combo_reg <- cbind(col_only_reg_df, c(col_nodysp_reg_df,col_dysp_reg_df))

##Plot the bar graph
par(mfrow = c(1,1))
barplot(as.matrix(combo_reg), col = c("steelblue","coral2"), ylim = c(0,500),legend.text = c("upregulated","downregulated"))

par(mfrow = c(1,3))
barplot(as.matrix(col_only_reg_df), 
        col = c("steelblue","coral2"), 
        ylim = c(0,500),
        axis.lty = 1,
        names.arg = c("HTG","HTSeq"),
        main = "Colitis Only",
        cex.main = 1.6,
        cex.axis = 1.5,
        cex.names = 1.5)
barplot(as.matrix(col_nodysp_reg_df), 
        col = c("steelblue","coral2"),
        ylim = c(0,500),
        axis.lty = 1,
        names.arg = c("HTG","HTSeq"),
        main = "Colitis and Dysplasia (-)",
        cex.main = 1.6,
        cex.axis = 1.5,
        cex.names = 1.5)

barplot(as.matrix(col_dysp_reg_df), 
        col = c("steelblue","coral2"), 
        ylim = c(0,500),
        axis.lty = 1,
        names.arg = c("HTG","HTSeq"),
        main = "Colitis and Dysplasia (+)",
        cex.main = 1.6,
        cex.axis = 1.5,
        cex.names = 1.5)

