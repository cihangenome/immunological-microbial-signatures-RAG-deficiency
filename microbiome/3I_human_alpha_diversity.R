#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
setDTthreads(4)


#' Read in metadata
A <- read.delim('tables/human_mapfile.txt', row.names = 1)

#' Remove these
sampREMV <- c("FP_26", "FP_46")
A <- subset(A, !SampleID %in% sampREMV)

#' ## Alpha diversity plots
#' - Alpha diversity is already calculated in the metadata table with vegan, so we can plot
cols <- c("#006E5E" , "#FD8D3C")
names(cols) <- c("CTRL", "RAG")



#' ## Stats
#' ### All Ages
#'
#' - CTRL vs all RAG1
kt <- kruskal.test(A$Shannon, A$SampleType)
print(kt)


#' - All samples
ggplot(A, aes(x=SampleType, y=Shannon, color = SampleType, fill=SampleType)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, size=2) +
  theme_bw(base_family = 'Arial') + 
  scale_color_manual(values=cols, guide = "none") +
  scale_fill_manual(values=add.alpha(cols, 0.4, noname = F), guide="none") + 
  xlab('Patient cohort') +
  ylab("Shannon diversity") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

