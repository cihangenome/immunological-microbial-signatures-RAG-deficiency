#' Statistical comparisons of taxa for timecourse/longitudinal mouse results.  Used for 
#' Figures 4A, Supp 4A/B

library(Maaslin2)
library(magrittr)
library(ampvis2)
#' From R package:  https://gitlab.com/pooranis/datavis16s
import::from('datavis16s/R/utilities.R', highertax, filterlowabund)


A <- read.delim('tables/timecourse_mapfile_mouse.txt')
otu <- read.delim('tables/timecourse_silva_ASV_table_mouse.txt', row.names = 1, check.names = F)

input <- amp_load(otutable = otu, metadata = A)

#' Filter OTU table
prev <- 5
taxlevel <- 'Family'
# or taxlevel <- 'Genus'
filterfun <- function(amp, sex) {
  ampsub <- amp_subset_samples(amp, Study == "timecourse" & Sex == sex & Genotype != "Q/Q mutant") %>%
    filterlowabund(level = 5, abs=T, persamp=prev)
  if (taxlevel != 'OTU') ampsub <- highertax(ampsub, taxlevel)

  return(ampsub)
}

female <- filterfun(input, 'female')



#' Run Maaslin2
formula = 'expr ~ Genotype + AgeInWeeks + (1 | SubjectID)'
min_variance = 0
transform = 'LOG'
normalization = 'TSS'
analysis_method = 'LM'

maaslinfun <- function(amp) {
  fit.data <- Maaslin2(amp$abund, amp$metadata,
                       fixed_effects = c( "Genotype", "AgeInWeeks"),
                       random_effects = c("SubjectID"),
                       plot_heatmap = F, plot_scatter = F, cores = 6,
                       output = "./maaslin", min_abundance = -1000,
                       min_prevalence = 0, min_variance = min_variance,
                       transform = transform,
                       normalization = normalization,
                       analysis_method = analysis_method,
                       full_formula = as.formula('expr ~ Genotype*AgeInWeeks + (1 + AgeInWeeks | SubjectID)')

                       )
  fit.data$results$feature <- plyr::mapvalues(fit.data$results$feature,
                                              make.names(amp$tax[[taxlevel]]), amp$tax[[taxlevel]])
  fit.data$results <- merge(fit.data$results, amp$tax, by.x='feature',
                             by.y=taxlevel)
  return(fit.data)
}

female.fit.rs <- maaslinfun(female)
saveRDS(female.fit.rs, file=paste0('./maaslin/female.fit.data.', taxlevel, '.rds'))
