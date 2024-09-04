extrafont::loadfonts('pdf', quiet=T)
library(ggplot2)
library(ampvis2)
library(lmerTest)

A <- read.delim("tables/timecourse_mapfile_mouse.txt")

otu <- read.delim('tables/timecourse_silva_ASV_table_mouse.txt', row.names = 1, check.names = F)

#' Calculate diversity rarefied to 40k counts
amp <- amp_load(otutable = otu, metadata = A)
ampsub <- amp_subset_samples(ampsub, rarefy=40000)
divmat <- data.frame(Shannon = vegan::diversity(ampsub$abund, MARGIN = 2), t(vegan::estimateR(t(ampsub$abund))))

divmat$Pielou <- divmat$Shannon/log(divmat$S.chao1)
divmat <- cbind.data.frame(divmat, ampsub$metadata)

cols <- c("#ef553b", "#636efa")
names(cols) <-  c("W/W mutant", "LRF WT(control)")


divplot <- ggplot(divmat, aes(x=as.factor(AgeInWeeks), y=Shannon, color=Genotype, fill=Genotype)) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  geom_point(position = position_jitterdodge(jitter.width=0.15)) +
  theme_bw(base_size = 8, base_family = 'Arial') +
  xlab("Age (weeks)") +
  scale_color_manual(values = cols, limits = c("LRF WT(control)", "W/W mutant"), labels = c("+/+", "W/W")) +
  scale_fill_manual(values=add.alpha(cols, 0.35, noname = F), guide="none") + 
  stat_summary(fun = "median", geom = 'line', mapping = aes(group = Genotype), position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
print(divplot)


#' Stats
cat('\n\nDoes the Shannon diversity vary significantly with Genotype?\n\n')
lm <- lmer(Shannon ~ Genotype + AgeInWeeks + (1 + AgeInWeeks | SubjectID), data=divmat)

cat('- model with all samples fails to converge with random slope.  Try without.\n\n')
lm <- lmer(Shannon ~ Genotype + AgeInWeeks + (1 | SubjectID), data=divmat)
lmsumm <- summary(lm)
print(lmsumm)

cat('Taking into account Genotype difference over time \n\n- Is genotype difference over time is significant (does diversity change over time differ by Genotype)? \n\n- yes - weakly significant\n\n')
time.lm <- lmer(Shannon ~ Genotype*AgeInWeeks + (1 | SubjectID), data=divmat)
time.lmsumm <- summary(time.lm)
print(time.lmsumm)


