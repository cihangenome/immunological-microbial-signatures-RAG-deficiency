#' PICRUSt was run in Nephele: https://nephele.niaid.nih.gov/pipeline_details/picrust2/
#' This script finds the pathways significantly different between the W/W and +/+ based
#' on abundance.  

library(data.table)

A <- readmapfile("tables/merged_metadata_mouse.txt")

P <- fread('tables/path_abun_unstrat_descrip.tsv.gz')
names(P)[3:17] <- gsub('.', '_', names(P)[3:17], fixed=T)
pmap <- P[,1:2]
pmap[,sani := make.names(pathway)]
P$description <- NULL

#' Maaslin parameters
normalization = 'NONE'
transform = 'NONE'
analysis_method = 'LM'
min_abundance = -1000
min_prevalence = 0.3
threshold <- 0.11
inputmat <- 'clrmat'
options(please_run = F)


# maaslin

#' ## Compare W mut vs wild-type
#'
dontrun({
  library(Maaslin2)
  pbopts <- getOption("pboptions")
  pbopts$nout <- 1
  options(pboptions = pbopts)
  
  ba <- list()
  ba$metadata <- setDF(A)
  printdebug(dim(P))
  ba$abund <- setDF(P[,ba$metadata$SampleID, with=FALSE], rownames = P$pathway)
  ## remove zero rows, if any
  v <- which(rowSums(ba$abund) > 0)
  ba$abund <- ba$abund[v,]
  printdebug(dim(ba$abund))
  w <- which(apply(ba$abund, 1, \(x) length(which(x>0))) > min_prevalence * nrow(ba$metadata))
  clr <- log1p(t(ba$abund))
  ba$clrmat <- t(clr - apply(clr, 1, mean))
  ba$clrmat <- ba$clrmat[w,]
  printdebug(dim(ba$clrmat))

  ba$metadata$Genotype <- factor(as.character(ba$metadata$Genotype))
  row.names(ba$metadata) <- ba$metadata$SampleID

  fit.data.ba <- Maaslin2(input_data = ba[[inputmat]],
                          input_metadata = ba$metadata,
                          min_abundance = min_abundance,
                          min_prevalence = min_prevalence,
                          min_variance = -1,
                          normalization = normalization,
                          transform = transform,
                          analysis_method = analysis_method,
                          fixed_effects = "Genotype",
                          output = './maaslin', cores = 4, plot_scatter = F)
  fit.data.ba$results <- merge(pmap,fit.data.ba$results, by.y='feature', by.x = 'sani', all.y = T)
  fit.data.ba$results$sani <- NULL
  saveRDS(namedlist(fit.data.ba, ba), "ivan_result/picrustrevcomp_PipelineResults.d40d29072d31/pwy_fit.data.rds")
})


#' We only look at significant pathways: 
#' sigpwys <- c('CENTFERM-PWY', 'P163-PWY', 'PWY-5676')
#' signames <- paste(sigpwys, c('pyruvate fermentation to butanoate', 'L-lysine fermentation to acetate and butanoate', 'acetyl-CoA fermentation to butanoate II'), sep=': ')
#' The results of these stats are in tables/picrust2_mouse_stats.xlsx
