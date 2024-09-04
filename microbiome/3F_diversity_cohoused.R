#' Alpha diversity is calculated by the Nephele (https://nephele.niaid.nih.gov) DADA2 pipeline.
#' Using the `amp_alphadiv` function from the *ampvis2* R package.  Data is rarefied at 19800 counts.

library(ggplot2)
library(data.table)

#' files
diversitytable <- "tables/diversity_cohoused_study_Rag1_pilot.csv"
shanplotpdf <- "cohoused_shannon_diversity.pdf"

#' Keep only female
divmat <- fread(diversitytable)

#' Stats comparing Genotypes (W mut and wildtype)
#'
#' - Shannon
kt <- kruskal.test(divmat$Shannon, divmat$Genotype)
print(kt)


#' boxplot
shanplot <- ggplot(divmat, aes(x=Genotype, y=Shannon, color=Genotype, fill=Genotype)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.25, size=1.5, stroke=1) + 
  scale_fill_manual(values=add.alpha(cols, 0.35, noname = F), guide="none") + 
  scale_color_manual(values = cols, guide="none") + 
  theme_bw(base_family = 'ArialMT') + 
  theme(strip.background = element_rect(fill="transparent")) + 
  scale_x_discrete(labels = c("+/+", "W/W")) +
  ylab('Shannon index') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(shanplot)


