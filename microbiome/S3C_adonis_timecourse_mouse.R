extrafont::loadfonts('pdf', quiet=T)
library(ggplot2)
library(ampvis2)
library(vegan)



#' Read in ASVs and tree
A <- read.delim( "tables/timecourse_mapfile_mouse.txt")
datafile <- tools::file_path_as_absolute("tables/timecourse_silva_ASV_table_mouse.txt")
otu <- read.delim(datafile, check.names = F)
row.names(otu) <- otu[,1]
otu[,1] <- NULL

#' Tree is made with QIIME2 align-to-tree-mafft-fasttree with default parameters
amp <- suppressWarnings(amp_load(otutable = otu, metadata = A, tree=ape::read.tree("tables/timecourse_mouse_tree.nwk")))

#' Filter taxa below 5 counts, only look at timecourse study, and remove Q/Q Female
ampll <- datavis16s:::filterlowabund(amp, level = 5, abs=T)
ampall <- amp_subset_taxa(amp, tax_vector = ampll$tax$OTU)


#' ## Compute Weighted UniFrac distance and plot beta diversity
colors <- c("#ef553b", "#00cc96", "#636efa",  "#19d3f3")
names(colors) <- c("W/W Female", "W/W Male", "wild-type Female", "wild-type Male")
labels <- setNames(c("W/W female", "W/W male", "+/+ female", "+/+ male"), nm = names(colors))

#' **Update** remove female from labels since we only include female in the paper
labels <- setNames(c("W/W", "W/W male", "+/+", "+/+ male"), nm = names(colors))

makepcoa <- function(amptc, age=T, subtitle=NULL, print=T, size=1.5, weighted = TRUE, wunitc=NULL, hull=F) {
  if (weighted) {
    if (is.null(wunitc)) wunitc <- amp_poorani(amptc, type="PCOA", distmeasure="wunifrac", transform = "none", detailed_output = T, sample_color_by = "Group", gunifrac_alpha = 1)
    # title <- "Weighted UniFrac beta diversity PCoA (2020)"
    title <- subtitle
  } else {
    if (is.null(wunitc)) wunitc <- amp_poorani(amptc, type="PCOA", distmeasure="unifrac", transform = "none", detailed_output = T, sample_color_by = "Group")
    # title <- "UniFrac beta diversity PCoA (2020)"
    title = subtitle
  }

  pcoawuni <- ggplot(wunitc$dsites, aes(x=PCo1, y=PCo2, color=Group))
  
  if (hull) {
    hulls <- wunitc$plot$data |>
      dplyr::group_by(Group) |>
      dplyr::slice(chull(PCo1, PCo2))
    pcoawuni <- pcoawuni + geom_polygon(data=hulls, aes(fill=Group), show.legend = F) 
  }
  
  
  if (age) {
    pcoawuni <- pcoawuni + geom_point(aes(shape=as.factor(AgeInWeeks)), size=size)
  } else {
    pcoawuni <- pcoawuni + geom_point(size=size)
  }
  mycolors <- colors[unique(wunitc$dsites$Group)]
  mylabels <- labels[unique(wunitc$dsites$Group)]
  pcoawuni <- pcoawuni + theme_bw(base_family = 'Arial') +
    scale_color_manual(values=mycolors, labels = mylabels, name='Genotype') +
    scale_shape_discrete(name="Age (weeks)") + xlab(paste0("PCo1 [", round(wunitc$screeplot$data$eigenvalues[1], 1), "%]")) + ylab(paste0("PCo2 [", round(wunitc$screeplot$data$eigenvalues[2], 1), "%]")) + ggtitle(title)
  
  if (hull) {
    pcoawuni <- pcoawuni + scale_fill_manual(values = add.alpha(mycolors, alpha=0.1, noname=F), guide='none')
  }
  
  if (print) print(pcoawuni)
  invisible(pcoawuni)
}

printadonis <- function(aaa, dat, distm="Jaccard") {
  cat("\n\n")
  cap <- paste(gsub('\n', '<br />', attr(aaa, 'heading'), fixed=T), collapse='')
  # add variables being compared to caption
  cap <- paste(cap, paste(levels(dat$metadata$Genotype), collapse = ', '), sep = '<br />') |>
    paste(paste0("Dist: ", distm), sep = '<br />')
  
  pander:::pander.anova(aaa, add.significance.stars = T, caption = cap, justify='lrrrrrl', missing='', width=120)
  cat("\n\n")
  
}

#' ### 3 weeks
amp3wf <- amp_subset_samples(ampall, Group %in% c("W/W Female", "wild-type Female") & AgeInWeeks == 3)
uni3wf <- amp_poorani(amp3wf, type="PCOA", distmeasure="wunifrac", transform = "none", detailed_output = T, sample_color_by = "Group", gunifrac_alpha = 1)
plot3wf <- makepcoa(amp3wf, age=F, subtitle = "3 weeks", size=2.5, print=T, wunitc = uni3wf)

plot3wfhull <- makepcoa(amp3wf, age=F, subtitle = "3 weeks", size=2.5, print=F, wunitc = uni3wf, hull=T)

#' Adonis
set.seed(50)
perm <- how(nperm = 99999)
set.seed(50)
adonis_3w <- adonis2(formula = uni3wf$inputmatrix ~ Genotype, data = amp3wf$metadata, permutations = perm, parallel = 4)
printadonis(adonis_3w, amp3wf)

#' ### 14 weeks
amp14wf <- amp_subset_samples(ampall, Group %in% c("W/W Female", "wild-type Female") & AgeInWeeks == 14)
uni14wf <- amp_poorani(amp14wf, type="PCOA", distmeasure="wunifrac", transform = "none", detailed_output = T, sample_color_by = "Group", gunifrac_alpha = 1)
plot14wf <- makepcoa(amp14wf, age=F, subtitle = "14 weeks", size=2.5, print=T, wunitc = uni14wf)


plot14wfhull <- makepcoa(amp14wf, age=F, subtitle = "14 weeks", size=2.5, print=F, wunitc = uni14wf, hull = T)



#' Adonis
set.seed(50)
perm <- how(nperm = 99999)
set.seed(50)
adonis_14w <- adonis2(formula = uni14wf$inputmatrix ~ Genotype, data = amp14wf$metadata, permutations = perm, parallel = 4)
printadonis(adonis_14w, amp14wf)



#' ### References
#'
#' adonis: McArdle, B.H. and M.J. Anderson. 2001. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. *Ecology*, **82**: 290–297.
#'
#' Oksanen J, Blanchet FG, Friendly M, Kindt R, Legendre P, McGlinn D, Minchin PR, O'Hara RB, Simpson GL, Solymos P, Stevens MHH, Szoecs E, Wagner H (2020).
#' <em>vegan: Community Ecology Package</em>.
#' R package version 2.5-7, <a href="https://CRAN.R-project.org/package=vegan">https://CRAN.R-project.org/package=vegan</a>.
#'
#' PCoA `amp_ordinate` function in *ampvis2* R package v2.7.9: Andersen KS, Kirkegaard RH, Karst SM, Albertsen M (2018).
#' &ldquo;ampvis2: an R package to analyse and visualise 16S rRNA amplicon data.&rdquo;
#' <em>bioRxiv</em>.
#' <a href="https://doi.org/10.1101/299537">https://doi.org/10.1101/299537</a>.
#'
#' UniFrac distance: Lozupone C, Knight R. UniFrac: a New Phylogenetic Method for Comparing Microbial Communities. <i>Appl Environ Microbiol</i> 2005;<b>71</b>:8228–35. doi:<a href="https://doi.org/10.1128/AEM.71.12.8228-8235.2005">10.1128/AEM.71.12.8228-8235.2005</a>
#'
