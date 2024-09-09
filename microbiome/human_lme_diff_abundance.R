library(data.table)
library(ampvis2)
#' uses functions from: https://gitlab.com/pooranis/datavis16s
library(Maaslin2)
pbopts <- getOption("pboptions")
pbopts$nout <- 1
options(pboptions = pbopts)


#' Read in metadata
A <- read.delim('tables/human_mapfile.txt', row.names = 1) |>
  mixedarrange(SampleType, HouseHoldID)
#' print categorical vars
table(A[,c("PhenoType", "AgeRange")])


#' Read in count data
otu <- fread('tables/human_asv_tab.txt') |>
  setnames('V1', "ASV")
tax <- fread('tables/human_tax_tab.txt', header = T, drop="V1")
amp <- amp_load(otu, A, tax)


#' ## Maaslin2 and LME
#' 
#' - HouseHoldID is a confounder, so we will adjust for it in our model.
#' - $\sim SampleType + (1 | HouseHoldID)
#' 
#' - I used CLR transformation, bc that's what we used in Mouse
#' 
#' qval and filter threshold & Maaslin parameters
#+ class.source='fold-none'
filterlevel <- 3
filterprev <- 20
normalization = 'NONE'
transform = 'NONE'
analysis_method = 'LM'
min_abundance = -1000
inputmat <- 'clrmat'
options(please_run = F)

#' ### Family level
dontrun(please_run = F, {
 famres <- local({
    ba <- amp
    
    
    ## filter prev
    ba <- highertax(ba, "Family")
    ba <- filterlowabund(ba, level = filterlevel, abs=T,  persamp = filterprev)
    ## map to get real names from matrix sanitized row.names that Maaslin2 uses
    map <- data.frame(taxa = ba$tax$Family, sani = make.names(ba$tax$Family))
    
    ## clr
    ba$clrmat <- log1p(t(ba$abund)) - log1p(apply(t(ba$abund), 1, mean))
    
    aaa <- capture.output(fit.data.ba <- Maaslin2(
      input_data = ba[[inputmat]],
      input_metadata = ba$metadata,
      min_abundance = min_abundance,
      min_prevalence = 0,
      min_variance = -1,
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method,
      fixed_effects = "SampleType",
      random_effects = c('HouseHoldID'),
      output = './maaslin',
      cores = 4,
      plot_scatter = F,
      plot_heatmap = F
    ))
    fit.data.ba$results <- subset(fit.data.ba$results, !grepl('unclassified', feature))
    fit.data.ba$results$feature <- map$taxa[match(fit.data.ba$results$feature, map$sani)]
    fit.data.ba$clrmat <- ba[[inputmat]]
    fit.data.ba$amp <- ba
    return(fit.data.ba)
  })
  
  saveRDS(famres, 'rds/maaslin_human_family.rds')
})
famres <- readRDS('rds/maaslin_human_family.rds')


#' ### Genus Level
dontrun(please_run = F, {
  genusres <- local({
    ba <- amp
    
    
    ## filter prev
    ba <- highertax(ba, "Genus")
    ba <- filterlowabund(ba, level = filterlevel, abs=T,  persamp = filterprev)
    ## map to get real names from matrix sanitized row.names that Maaslin2 uses
    map <- data.frame(taxa = ba$tax$Genus, sani = make.names(ba$tax$Genus))
    
    ## clr
    ba$clrmat <- log1p(t(ba$abund)) - log1p(apply(t(ba$abund), 1, mean))
    
    aaa <- capture.output(fit.data.ba <- Maaslin2(
      input_data = ba[[inputmat]],
      input_metadata = ba$metadata,
      min_abundance = min_abundance,
      min_prevalence = 0,
      min_variance = -1,
      normalization = normalization,
      transform = transform,
      analysis_method = analysis_method,
      fixed_effects = "SampleType",
      random_effects = c('HouseHoldID'),
      output = './maaslin',
      cores = 4,
      plot_scatter = F,
      plot_heatmap = F
    ))
    fit.data.ba$results <- subset(fit.data.ba$results, !grepl('unclassified', feature))
    fit.data.ba$results$feature <- map$taxa[match(fit.data.ba$results$feature, map$sani)]
    fit.data.ba$clrmat <- ba[[inputmat]]
    fit.data.ba$amp <- ba
    return(fit.data.ba)
  })
  
  saveRDS(genusres, 'rds/maaslin_human_genus.rds')
})
genusres <- readRDS('rds/maaslin_human_genus.rds')
