#!/usr/bin/env Rscript
library(ggplot2)

# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
file = args[1]
Threshold = as.numeric(args[2])
out = args[3]

# -------- sample name ---------
sample = sub("_total_pi_est.csv", "",basename(file))
# -------- get data ------------
Pi <- as.data.frame(read.delim(file,header=T,sep=","))
Pi$Idx=NULL
colnames(Pi)=c("references","proportion","StdFreq")
Pi = Pi[Pi[,2]>0.00001,]
Pi = Pi[order(-Pi[,2]),]

# some ref are too big
Pi$references = sub(" .*", "", Pi$references)

# color columns differently depending on Thresholding
Pi$color = Pi[,2]>Threshold
Pi$color[Pi$color] = "bold"
Pi$color[Pi$color==FALSE] = "plain"
# order levels so that the manual color scale always colors "plain" in the same color, which is grey
Pi$color = factor(Pi$color,levels=c("bold","plain"))

# -------- plot  ------------
# for visualisation purpose, ignore anything but the 50 first most abundant variant
Pi_small = Pi[1:min(dim(Pi)[1],50),]
p <- ggplot(Pi_small, aes(x=reorder(references,-proportion), y=proportion, fill=color)) + geom_bar(stat="identity")+scale_fill_manual(name="color",values=c("#00BFC4","grey"))
p <- p + theme_bw() + geom_hline(yintercept=Threshold, linetype="dashed", color = "red") 
p <- p + ylab('Predicted variant proportions') + xlab('Reference variants') 
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face =Pi_small$color),legend.position = "none")
p <- p +ggtitle(sample)

# -------- save  ------------
pdf(out)
plot(p)
dev.off()
