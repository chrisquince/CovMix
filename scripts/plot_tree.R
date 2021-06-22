#!/usr/bin/env Rscript
library(ggplot2)
library(ggtree)

# -------- get args ------------
args = commandArgs(trailingOnly = TRUE)
file = args[1]
tree_file = args[2]
out = args[3]
sample = sub("_total_.*", "",basename(file))

# -------- get frequencies ------------
Pi <- read.csv(file,header=T)
Pi$MeanFreq <- as.numeric(Pi$MeanFreq)
Pi <- Pi[complete.cases(Pi$MeanFreq),]

# -------- get tree ------------
tree <- read.tree(tree_file)
p <- ggtree(tree)
rownames(Pi) <- Pi$Ref

# -------- plot tree+frequencies ------------
d2 <- data.frame(id=tree$tip.label, value = Pi[tree$tip.label,]$MeanFreq)
p3 <- facet_plot(p, panel='Frequency', data=d2, geom=geom_segment, 
                 aes(x=0, xend=value, y=y, yend=y), size=1, color="#00BFC4")
p3 = p3 + ggtitle(sample)
pdf(out)
plot(p3)
dev.off()