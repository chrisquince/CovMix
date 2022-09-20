#!/usr/bin/env Rscript
library(ggplot2)

# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
pi_file = args[1]
ambi_file = args[2]
Threshold = as.numeric(args[3])
out = args[4]

# -------- sample name ---------
sample = sub("_total_pi_est.csv", "",basename(pi_file))

# -------- get data ------------
Pi <- as.data.frame(read.delim(pi_file,header=T,sep=","))
Pi$Idx=NULL
colnames(Pi)=c("references","proportion","StdFreq")
Pi = Pi[Pi[,2]>0.00001,]
Pi = Pi[order(-Pi[,2]),]

# deal with ambiguity
ambi <- as.data.frame(read.delim(ambi_file,header=T,sep=",",row=2))
is_degenerate = ambi[Pi$references,2]>1
# some ref are too big
Pi$references = sub(" .*", "", Pi$references)
Pi$references[is_degenerate] = paste(Pi$references[is_degenerate]," (",ambi[Pi$references[is_degenerate],2],")",sep="")

# color columns differently depending on Thresholding
filtered = Pi[,2]>Threshold
Pi$style = filtered
Pi$style[filtered] = "bold"
Pi$style[!filtered] = "plain"

# order levels so that the manual color scale always colors "plain" in the same color, which is grey
Pi$color = filtered
Pi$color[filtered] = "#00BFC4"
Pi$color[filtered&is_degenerate] = "#00A36C"
Pi$color[!filtered] = "gray"

group.colors <- c("#00BFC4" = "#00BFC4", "#00A36C" = "#00A36C", "gray" ="gray")

# -------- plot  ------------
# for visualisation purpose, ignore anything but the 50 first most abundant variant
Pi_small = Pi[1:min(dim(Pi)[1],50),]
p <- ggplot(Pi_small, aes(x=reorder(references,-proportion), y=proportion, fill=color)) + geom_bar(stat="identity")+scale_fill_manual(name="color",values=group.colors)
p <- p + theme_bw() + geom_hline(yintercept=Threshold, linetype="dashed", color = "red") 
p <- p + ylab('Predicted variant proportions') + xlab('Reference variants') 
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face =Pi_small$style),legend.position = "none")
p <- p +ggtitle(sample)

# -------- save  ------------
pdf(out)
plot(p)
dev.off()
