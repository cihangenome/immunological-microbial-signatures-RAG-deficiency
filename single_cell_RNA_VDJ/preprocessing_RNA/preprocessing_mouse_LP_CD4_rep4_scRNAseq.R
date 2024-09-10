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


matrix_dir_HTO="./HTO_counts_rep4/"
barcode.path <- paste0(matrix_dir_HTO, "barcodes.tsv")
features.path <- paste0(matrix_dir_HTO, "features.tsv")
matrix.path <- paste0(matrix_dir_HTO, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


LP_CD4_rep4=Read10X_h5("./RNA_counts_rep4/raw_feature_bc_matrix.h5")

colnames(mat) <- paste(colnames(mat),"-1", sep ="")

LP_CD4_rep4<-LP_CD4_rep4[,order(colnames(LP_CD4_rep4))]
mat<-mat[,order(colnames(mat))]

joint_bcs <- intersect(colnames(LP_CD4_rep4),colnames(mat))

length(colnames(LP_CD4_rep4))
length(colnames(mat))
length(joint_bcs)
##################################

LP_CD4_rep4 <- LP_CD4_rep4[,joint_bcs]
hto_LP_CD4_rep4<- as.matrix(mat)
hto_LP_CD4_rep4<-hto_LP_CD4_rep4[,joint_bcs]

##################################
#!#################################
hto_LP_CD4_rep4<-hto_LP_CD4_rep4[2:7,]##this step removes HTO1
#!#################################
rownames(hto_LP_CD4_rep4)<-c("HTO_2","HTO_3","HTO_4","HTO_5","HTO_6","HTO_7") 
#!#################################
##################################

LP_CD4_rep4 <- CreateSeuratObject(counts = LP_CD4_rep4, min.cells = 100, min.features = 200)


LP_CD4_rep4[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep4, pattern = "^mt-")


LP_CD4_rep4_sc <- as.SingleCellExperiment(LP_CD4_rep4)
LP_CD4_rep4_sc <- scDblFinder(LP_CD4_rep4_sc)

ind_doublets_LP_CD4_rep4<-which(LP_CD4_rep4_sc$scDblFinder.class=="doublet")
LP_CD4_rep4_sc <- LP_CD4_rep4_sc[,-ind_doublets_LP_CD4_rep4]

high_mitoexp_LP_CD4_rep4<-isOutlier(LP_CD4_rep4_sc@colData@listData$percent.mt, nmads=4, type="higher")
high_ncount_LP_CD4_rep4<-isOutlier(LP_CD4_rep4_sc@colData@listData$nCount_RNA, nmads=2.5, type="higher")
low_ncount_LP_CD4_rep4<-isOutlier(LP_CD4_rep4_sc@colData@listData$nCount_RNA, nmads=3, type="lower")
high_nfeature_LP_CD4_rep4<-isOutlier(LP_CD4_rep4_sc@colData@listData$nFeature_RNA, nmads=2.5, type="higher")
low_nfeature_LP_CD4_rep4<-isOutlier(LP_CD4_rep4_sc@colData@listData$nFeature_RNA, nmads=3, type="lower")

ind_outliers_LP_CD4_rep4<-Reduce("|", list(high_mitoexp_LP_CD4_rep4,high_ncount_LP_CD4_rep4,low_ncount_LP_CD4_rep4,high_nfeature_LP_CD4_rep4,low_nfeature_LP_CD4_rep4))
table(ind_outliers_LP_CD4_rep4)


ind_outliers_LP_CD4_rep4<-which(ind_outliers_LP_CD4_rep4=="TRUE")


LP_CD4_rep4_sc <- LP_CD4_rep4_sc[,-ind_outliers_LP_CD4_rep4]


LP_CD4_rep4 <- subset(LP_CD4_rep4, cells=LP_CD4_rep4_sc@colData@rownames)


# Normalize RNA data with log normalization
LP_CD4_rep4 <- NormalizeData(LP_CD4_rep4)
# Find and scale variable features
LP_CD4_rep4 <- FindVariableFeatures(LP_CD4_rep4,selection.method = "vst",mean.function = ExpMean)
LP_CD4_rep4 <- ScaleData(LP_CD4_rep4, features = rownames(LP_CD4_rep4))

# Add HTO data as a new assay independent from RNA
LP_CD4_rep4[["HTO"]] <- CreateAssayObject(counts = hto_LP_CD4_rep4[, colnames(x = LP_CD4_rep4)])
# Normalize HTO data,here we use centered log-ratio (CLR) transformation
LP_CD4_rep4 <- NormalizeData(LP_CD4_rep4, assay = "HTO", normalization.method = "CLR" )

LP_CD4_rep4_demux <- HTODemux(LP_CD4_rep4, assay = "HTO", positive.quantile = 0.95)

table(LP_CD4_rep4_demux$HTO_classification.global)

table(LP_CD4_rep4_demux@meta.data[["HTO_maxID"]])


Idents(LP_CD4_rep4_demux) <- "HTO_classification.global"


# First, we will remove negative cells from the object
LP_CD4_rep4_demux.subset <- subset(LP_CD4_rep4_demux, idents = "Negative", invert = TRUE)

DefaultAssay(LP_CD4_rep4_demux.subset) <- "HTO"
LP_CD4_rep4_demux.subset <- ScaleData(LP_CD4_rep4_demux.subset, features = rownames(LP_CD4_rep4_demux.subset),
    verbose = FALSE)
LP_CD4_rep4_demux.subset <- RunPCA(LP_CD4_rep4_demux.subset, features = rownames(LP_CD4_rep4_demux.subset), approx = FALSE)


# Extract the singlets
LP_CD4_rep4.singlet <- subset(LP_CD4_rep4_demux, idents = "Singlet")
LP_CD4_rep4.singlet[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep4.singlet, pattern = "^mt-")

LP_CD4_rep4_filtered<-LP_CD4_rep4.singlet

table(LP_CD4_rep4_filtered@meta.data[["HTO_maxID"]])
table(LP_CD4_rep4_filtered@meta.data[["HTO_classification"]])

############################
#!###########################
LP_CD4_rep4_filtered@meta.data$orig.ident<-as.character(LP_CD4_rep4_filtered@meta.data$orig.ident)
ind_Wmutant<-which(LP_CD4_rep4_filtered@meta.data$HTO_maxID %in% c("HTO-3","HTO-5","HTO-7")) #"HTO-1", this is an exception
ind_WT<-which(LP_CD4_rep4_filtered@meta.data$HTO_maxID %in% c("HTO-2","HTO-4","HTO-6"))
LP_CD4_rep4_filtered@meta.data$orig.ident[ind_WT]<-"WT_12wks"
LP_CD4_rep4_filtered@meta.data$orig.ident[ind_Wmutant]<-"WW_12wks"
#!###########################
############################

Idents(LP_CD4_rep4_filtered) <- "orig.ident"


saveRDS(LP_CD4_rep4_filtered,"LP_CD4_rep4_filtered.rds")
