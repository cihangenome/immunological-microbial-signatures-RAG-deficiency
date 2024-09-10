#setwd("current_directory")

library(dplyr)
library(Seurat)
library(patchwork)
library(scran)

library(scDblFinder)

library(Matrix)
library(writexl)
library(readxl)

library(scater)
library(loomR)
library(destiny)
library(ggplot2)
library(ggthemes)
library(SingleCellExperiment)
library(ggpubr)


matrix_dir_HTO="./HTO_counts_rep6/"
barcode.path <- paste0(matrix_dir_HTO, "barcodes.tsv")
features.path <- paste0(matrix_dir_HTO, "features.tsv")
matrix.path <- paste0(matrix_dir_HTO, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


LP_CD4_rep6=Read10X_h5("./RNA_counts_rep6/raw_feature_bc_matrix.h5")

colnames(mat) <- paste(colnames(mat),"-1", sep ="")

LP_CD4_rep6<-LP_CD4_rep6[,order(colnames(LP_CD4_rep6))]
mat<-mat[,order(colnames(mat))]

joint_bcs <- intersect(colnames(LP_CD4_rep6),colnames(mat))

length(colnames(LP_CD4_rep6))
length(colnames(mat))
length(joint_bcs)
##################################

LP_CD4_rep6 <- LP_CD4_rep6[,joint_bcs]
hto_LP_CD4_rep6<- as.matrix(mat)
hto_LP_CD4_rep6<-hto_LP_CD4_rep6[,joint_bcs]


##################################
#!#################################
hto_LP_CD4_rep6<-hto_LP_CD4_rep6[1:5,]##
#!#################################
rownames(hto_LP_CD4_rep6)<-c("HTO_1","HTO_2","HTO_3","HTO_4","HTO_5")
##################################
#!#################################

LP_CD4_rep6 <- CreateSeuratObject(counts = LP_CD4_rep6, min.cells = 100, min.features = 200)


LP_CD4_rep6[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep6, pattern = "^mt-")


LP_CD4_rep6_sc <- as.SingleCellExperiment(LP_CD4_rep6)
LP_CD4_rep6_sc <- scDblFinder(LP_CD4_rep6_sc)

ind_doublets_LP_CD4_rep6<-which(LP_CD4_rep6_sc$scDblFinder.class=="doublet")
LP_CD4_rep6_sc <- LP_CD4_rep6_sc[,-ind_doublets_LP_CD4_rep6]

high_mitoexp_LP_CD4_rep6<-isOutlier(LP_CD4_rep6_sc@colData@listData$percent.mt, nmads=4, type="higher")
high_ncount_LP_CD4_rep6<-isOutlier(LP_CD4_rep6_sc@colData@listData$nCount_RNA, nmads=2.5, type="higher")
low_ncount_LP_CD4_rep6<-isOutlier(LP_CD4_rep6_sc@colData@listData$nCount_RNA, nmads=3, type="lower")
high_nfeature_LP_CD4_rep6<-isOutlier(LP_CD4_rep6_sc@colData@listData$nFeature_RNA, nmads=2.5, type="higher")
low_nfeature_LP_CD4_rep6<-isOutlier(LP_CD4_rep6_sc@colData@listData$nFeature_RNA, nmads=3, type="lower")

ind_outliers_LP_CD4_rep6<-Reduce("|", list(high_mitoexp_LP_CD4_rep6,high_ncount_LP_CD4_rep6,low_ncount_LP_CD4_rep6,high_nfeature_LP_CD4_rep6,low_nfeature_LP_CD4_rep6))
table(ind_outliers_LP_CD4_rep6)


ind_outliers_LP_CD4_rep6<-which(ind_outliers_LP_CD4_rep6=="TRUE")


LP_CD4_rep6_sc <- LP_CD4_rep6_sc[,-ind_outliers_LP_CD4_rep6]

plot1 <- FeatureScatter(LP_CD4_rep6, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by="orig.ident")
plot2 <- FeatureScatter(LP_CD4_rep6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by="orig.ident")
CombinePlots(plots = list(plot1, plot2))

LP_CD4_rep6 <- subset(LP_CD4_rep6, cells=LP_CD4_rep6_sc@colData@rownames)


# Normalize RNA data with log normalization
LP_CD4_rep6 <- NormalizeData(LP_CD4_rep6)
# Find and scale variable features
LP_CD4_rep6 <- FindVariableFeatures(LP_CD4_rep6,selection.method = "vst",mean.function = ExpMean)
LP_CD4_rep6 <- ScaleData(LP_CD4_rep6, features = rownames(LP_CD4_rep6))

# Add HTO data as a new assay independent from RNA
LP_CD4_rep6[["HTO"]] <- CreateAssayObject(counts = hto_LP_CD4_rep6[, colnames(x = LP_CD4_rep6)])
# Normalize HTO data,here we use centered log-ratio (CLR) transformation
LP_CD4_rep6 <- NormalizeData(LP_CD4_rep6, assay = "HTO", normalization.method = "CLR" )

LP_CD4_rep6_demux <- HTODemux(LP_CD4_rep6, assay = "HTO", positive.quantile = 0.95)

table(LP_CD4_rep6_demux$HTO_classification.global)

table(LP_CD4_rep6_demux@meta.data[["HTO_maxID"]])


Idents(LP_CD4_rep6_demux) <- "HTO_classification.global"


# First, we will remove negative cells from the object
LP_CD4_rep6_demux.subset <- subset(LP_CD4_rep6_demux, idents = "Negative", invert = TRUE)

DefaultAssay(LP_CD4_rep6_demux.subset) <- "HTO"
LP_CD4_rep6_demux.subset <- ScaleData(LP_CD4_rep6_demux.subset, features = rownames(LP_CD4_rep6_demux.subset),
    verbose = FALSE)
LP_CD4_rep6_demux.subset <- RunPCA(LP_CD4_rep6_demux.subset, features = rownames(LP_CD4_rep6_demux.subset), approx = FALSE)


# Extract the singlets
LP_CD4_rep6.singlet <- subset(LP_CD4_rep6_demux, idents = "Singlet")
LP_CD4_rep6.singlet[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep6.singlet, pattern = "^mt-")

LP_CD4_rep6_filtered<-LP_CD4_rep6.singlet

############################
#!###########################

table(LP_CD4_rep6_filtered@meta.data[["HTO_maxID"]])
table(LP_CD4_rep6_filtered@meta.data[["HTO_classification"]])


LP_CD4_rep6_filtered@meta.data$orig.ident<-as.character(LP_CD4_rep6_filtered@meta.data$orig.ident)

ind_Wmutant_ANMV<-which(LP_CD4_rep6_filtered@meta.data$HTO_maxID %in% c("HTO-1","HTO-2","HTO-3"))
ind_Wmutant<-which(LP_CD4_rep6_filtered@meta.data$HTO_maxID %in% c("HTO-4","HTO-5"))

LP_CD4_rep6_filtered@meta.data$orig.ident[ind_Wmutant_ANMV]<-"WW_ANMV_12wks"
LP_CD4_rep6_filtered@meta.data$orig.ident[ind_Wmutant]<-"WW_CTRL_12wks"


Idents(LP_CD4_rep6_filtered) <- "orig.ident"

saveRDS(LP_CD4_rep6_filtered,"LP_CD4_rep6_filtered.rds")
