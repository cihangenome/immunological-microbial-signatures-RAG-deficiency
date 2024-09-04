library(ampvis2)
library(data.table)
library(lmerTest)


#' Uses wrapper around DADA2 for the Nephele pipeline: https://gitlab.com/pooranis/dada2nephele
#' to re-assign ASVs to the Mouse Oral Microbiome Database v5.1 https://momd.org/ftp/16S_rRNA_refseq/MOMD_16S_rRNA_RefSeq/V5.1/
#' Example code:
seqtab <- readRDS('boot60/intermediate_files/seqtab_min75.rds')

dontrun({
  library(dada2nephele)
  taxa <- taxonomy_assign(seqtab = seqtab,
                          refdb = 'MOMD_16S_rRNA_RefSeq_V5.1.train_set.fa',
                          refdb_species = 'MOMD_16S_rRNA_RefSeq_V5.1.species_assignment.fa',
                          nthread = 4,
                          taxmethod = 'rdp',
                          data_type = 'PE',
                          justConcatenate = FALSE,
                          minBoot = 60)
  dada2alloutput(seqtab, taxtab = taxa, outdir = 'momd')
})

A <- read.delim("tables/timecourse_mapfile_mouse.txt")
otu <- read.delim('tables/timecourse_momd_ASV_table_mouse.txt', row.names = 1, check.names = F)

amp <- suppressWarnings(amp_load(otutable = otu, metadata = A))

#' uses functions from: https://gitlab.com/pooranis/datavis16s wrapper around ampvis2
ampmomd <- datavis16s::highertax(ampmomd, 'Genus', keepunclass = T)
ampmomd$tax$Genus[which(ampmomd$tax$Genus %in% c('Muribaculaceae_[G-1]', 'Muribaculaceae_[G-2]'))] <- 'unclassified'
toplotmomd2 <- data.table(cbind(ampmomd$abund, ampmomd$tax))
toplotmomd2$oral <- ifelse(toplotmomd2$Genus == 'unclassified', 'oral_no', 'oral_yes')
toplotmomd2 <- toplotmomd2[,lapply(.SD, sum),by=.(oral), .SDcols=ampmomd$metadata$SampleID]
toplotmomd2[, 2:ncol(toplotmomd2)] <- 100 * vegan::decostand(toplotmomd2[,-1, with=FALSE], method='total', MARGIN = 2)
toplotmomd2 <- transpose(toplotmomd2, keep.names = 'SampleID', make.names = 'oral') |>
  merge(A,  by='SampleID', all.x = T) 
toplotmomd2$Group <- plyr::mapvalues(toplotmomd2$Group, c('wild-type Female', 'W/W Female'), c('+/+', 'W/W'))
toplotmomd2$Group <- factor(toplotmomd2$Group, levels = c('+/+', 'W/W'))

cols1 <- c( "#ef553b", "#636efa")
names(cols1) <-  c("W/W", "+/+")

bp <- ggplot(toplotmomd2, aes(x=as.factor(AgeInWeeks), y=oral_yes, color=Group, fill=Group)) + 
  geom_boxplot(outlier.shape = NA, show.legend = F) + geom_point(size=1.5, position = position_jitterdodge(jitter.width=0)) +
  theme_bw(base_family = 'Arial') +
  scale_color_manual(values = cols1, limits=c('+/+', 'W/W'), name='Genotype') +
  scale_fill_manual(values = add.alpha(cols1, alpha=0.35, noname=F), limits=c('+/+', 'W/W'), name='Genotype') +
  ylab('oral bacteria relative abundance (%)') + 
  xlab('age (weeks)') +
  stat_summary(fun = "median", geom = 'line', mapping = aes(group = Genotype), position = position_dodge(width = 0.9)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))



cat('\n\nDoes the oral bacteria vary significantly with Genotype?\n\n- Yes\n\n')
lm <- lmer(oral_yes ~ Genotype + AgeInWeeks + (1 | SubjectID), data=toplotmomd2)
lmsumm <- summary(lm)
print(lmsumm)

cat('Taking into account Genotype difference over time \n\n- Is genotype difference over time is significant (does the change in oral bacteria abundance over time differ by Genotype)?\n\n- No\n\n')
time.lm <- lmer(oral_yes ~ Genotype*AgeInWeeks + (1 | SubjectID), data=toplotmomd2)

time.lmsumm <- summary(time.lm)
print(time.lmsumm)

