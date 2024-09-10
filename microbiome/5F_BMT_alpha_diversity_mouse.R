#!/usr/bin/env Rscript
library(ggplot2)
library(data.table)
cols3 <- c("#636efa", "#CC0303", "#FF6363")

import::here(magrittr, `%>%`)

#' diversity is calculated in Nephele DADA2 pipeline using ampvis2 amp_alphadiv See https://nephele.niaid.nih.gov/user_guide_visualization_details/
A <- read.delim('tables/mouse_bmt_mapfile.txt')

div$TreatmentGroup <- 
  gsub(" (no irradiation)", "", div$TreatmentGroup,
       fixed = T)
div <- dplyr::arrange(div, SubjectID)
indexes <- c("Shannon", "Chao1")
meltdiv <- reshape2::melt(div, measure.vars = c("Shannon", "Chao1"))
meltdiv$TreatmentGroup <-
  factor(meltdiv$TreatmentGroup,
         levels = c("wild-type",
                    "Before BMT", "After BMT"))
divplot <- ggplot(meltdiv,
                  aes(
                    x = TreatmentGroup,
                    y = value,
                    color = TreatmentGroup,
                    fill = TreatmentGroup
                  )) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  scale_x_discrete(labels = \(x) relabgroups[x]) +
  scale_color_manual(values = cols3) + 
  scale_fill_manual(values = add.alpha(cols3, 0.35, noname = T), guide = "none") + 
  theme_bw(base_size = 8,base_family = "Arial") + 
  theme(legend.position = "none", strip.background = element_rect(fill = "transparent")) +
  xlab("Group") + ylab("index value") + 
  facet_wrap( ~ variable, scales = "free_y")
divplot1 <- divplot + geom_jitter(width = 0.2, size = 1.5)

print(divplot1)


shandiv <- subset(meltdiv, variable == "Shannon")
shanplot <- ggplot(shandiv, aes(x = TreatmentGroup, y=value, fill=TreatmentGroup, color=TreatmentGroup)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  geom_jitter(width = 0.2, size = 2) +
  scale_x_discrete(labels = \(x) relabgroups[x]) +
  scale_color_manual(values = cols3) + 
  scale_fill_manual(values = add.alpha(cols3, 0.4, noname = T), guide = "none") + 
  theme_bw(base_family = "Arial") + 
  theme(legend.position = "none", strip.background = element_rect(fill = "transparent")) +
  xlab("Group") + ylab("Shannon index") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

print(shanplot)

ggsave('divplot2.pdf', shanplot, width = 3.6, height=5) 
embed_fonts('divplot2.pdf')



## pre vs post
divpvp <- dcast(data = data.table(subset(div, TreatmentGroup !=
                                           "wild-type")), formula = SubjectID ~ TreatmentGroup, value.var = c("Shannon",
                                                                                                              "Chao1"))
testpvp <- coin::wilcoxsign_test(`Shannon_Before BMT` ~ `Shannon_After BMT`,
                                 data = divpvp)
print(testpvp)
testpvpchao1 <- coin::wilcoxsign_test(`Chao1_Before BMT` ~ `Chao1_After BMT`,
                                      data = divpvp)
print(testpvpchao1)

## wt vs pre
divwtvpre <- subset(div, TreatmentGroup != "After BMT")
divwtvpre <- setclasses(divwtvpre, list(factor = "TreatmentGroup"))
print(coin::wilcox_test(Shannon ~ TreatmentGroup, divwtvpre))
print(coin::wilcox_test(Chao1 ~ TreatmentGroup, divwtvpre))


## wt vs post
divwtvpost <- subset(div, TreatmentGroup != "Before BMT")
divwtvpost <- setclasses(divwtvpost, list(factor = "TreatmentGroup"))
print(coin::wilcox_test(Shannon ~ TreatmentGroup, divwtvpost))
print(coin::wilcox_test(Chao1 ~ TreatmentGroup, divwtvpost))
