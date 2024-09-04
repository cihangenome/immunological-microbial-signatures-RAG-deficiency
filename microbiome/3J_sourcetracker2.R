#' Takes in output of SourceTracker2 and compares abundances of ORAL vs other enterotypes
#' For human microbiome data. See *exe_SourceTracker2.sh* in this folder.
#' For figures 3J and supplemental figure 3E


library(data.table)
library(ggplot2)
extrafont::loadfonts('pdf', quiet=T)



#' Read in metadata
A <- read.delim('tables/human_mapfile.txt', row.names = 1)

#' read in proportions
props <- fread('tables/mixing_proportions.txt') |>
  melt(id.vars = 'V1', variable.name = 'SampleID', value.name = 'source_fraction') |>
  merge(A, by='SampleID') |>
  subset(V1 != "Unknown")

#' Stats
#' 
#' - CTRL vs all RAG1 - oral only
oralprops <- subset(props, V1 == "ORAL")
oralkt <- kruskal.test(oralprops$source_fraction, oralprops$SampleType)
print(oralkt)

#' - all enterotypes
statdf <- data.frame()
for (ent in unique(props$V1)) {
  subprops <- subset(props, V1 == ent)
  kt <- kruskal.test(subprops$source_fraction, subprops$SampleType)
  subkt <- data.frame(V1 = ent, statistic = kt$statistic, df = kt$parameter, p.value = kt$p.value, pval=paste("P =",signif(kt$p.value, digits = 2)))
  statdf <- rbind.data.frame(statdf, subkt)
}

print(statdf)



#' boxplot
cols <- c("#006E5E" , "#FD8D3C")
names(cols) <- c("CTRL", "RAG")

ggplot(props, aes(x=SampleType, y=source_fraction, color=SampleType, fill=SampleType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size=2) +
  facet_wrap(vars(V1), scale="free_y", nrow=1) +
  theme_bw(base_family = 'Arial') + 
  scale_color_manual(values=cols, guide = "none") +
  scale_fill_manual(values=add.alpha(cols, 0.4, noname = F), guide="none") + 
  theme(strip.background = element_rect(fill='white')) +
  scale_y_continuous(trans = 'sqrt', expand = expansion(mult = c(0.05, 0.1))) +
  ylab('source fraction') +
  xlab('Patient cohort') 

ggsave("human_st2.pdf", width=5.5, height=4.2)
embed_fonts("human_st2.pdf")


ggplot(oralprops, aes(x=SampleType, y=source_fraction, color=SampleType, fill=SampleType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size=2) +
  theme_bw(base_family = 'Arial') + 
  scale_color_manual(values=cols, guide = "none") +
  scale_fill_manual(values=add.alpha(cols, 0.4, noname = F), guide="none") + 
  theme(strip.background = element_rect(fill='white')) +
  scale_y_continuous(trans = 'sqrt', expand = expansion(mult = c(0.05, 0.1))) +
  ylab('oral source fraction') +
  xlab('Patient cohort') 
ggsave("human_st2_oral.pdf", width=2, height=4)
embed_fonts("human_st2_oral.pdf")

