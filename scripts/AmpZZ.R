#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
file = args[1]
out = args[2]


sample = sub("_total_*", "",basename(file))

# -------- things ------------
AmpFreq <- read.csv(file,header=T)

nA <- ncol(AmpFreq)

AmpFreq[,3:nA] <- lapply(AmpFreq[,3:nA], function(x) {
  as.numeric(as.character(x))
})
sapply(AmpFreq[,3:nA], class)

AmpFreq <- AmpFreq[ , colSums(is.na(AmpFreq)) == 0]

AmpFreq$Idx <- NULL

rownames(AmpFreq) <- AmpFreq$Ref

AmpFreq$Ref <- NULL

AmpFreq <- AmpFreq[rowSums(AmpFreq)>0,]

AmpFreq <- AmpFreq[,colSums(AmpFreq)>100]

AmpFreq <- t(AmpFreq)

AmpFreqP <- AmpFreq/rowSums(AmpFreq)

AmpFreqP <- t(AmpFreqP)


AmpFreqP_R <- transform(AmpFreqP,SD=apply(AmpFreqP,1, sd, na.rm = TRUE))

AmpFreqP_R <- transform(AmpFreqP_R,Mean=apply(AmpFreqP,1, mean, na.rm = TRUE))

AmpFreqP_R$R <- AmpFreqP_R$SD/AmpFreqP_R$Mean


AmpFreqP_melt <- melt(AmpFreqP)

colnames(AmpFreqP_melt) <- c('Genome','Amplicon','PAbund')

AmpFreqP_melt$Genome <- factor(AmpFreqP_melt$Genome,
                            levels = rownames(AmpFreqP[order(rowMeans(AmpFreqP)),]),
                            ordered = TRUE)



ggheatmap <- ggplot(AmpFreqP_melt, aes(Amplicon, Genome, fill = PAbund))
ggheatmap <- ggheatmap +  geom_tile(color = "white") + scale_fill_gradientn(trans = "sqrt",colours = c('white','blue','green','orange','red'))

ggheatmap <- ggheatmap  + theme(axis.text.x = element_text(size=6,angle = 45, vjust = 0.5, hjust=1))
ggheatmap <- ggheatmap +ggtitle(sample)

pdf(out)
print(ggheatmap)
dev.off()

