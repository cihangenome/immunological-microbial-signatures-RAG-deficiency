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

matrix_dir_rep1="./HTO_counts_rep1/"
barcode.path <- paste0(matrix_dir_rep1, "barcodes.tsv")
features.path <- paste0(matrix_dir_rep1, "features.tsv")
matrix.path <- paste0(matrix_dir_rep1, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


LP_CD4_rep1=Read10X_h5("./RNA_counts_rep1/raw_feature_bc_matrix.h5")

colnames(mat) <- paste(colnames(mat),"-1", sep ="")

head(colnames(LP_CD4_rep1))


head(colnames(mat))


LP_CD4_rep1<-LP_CD4_rep1[,order(colnames(LP_CD4_rep1))]
mat<-mat[,order(colnames(mat))]

joint_bcs <- intersect(colnames(LP_CD4_rep1),colnames(mat))

length(colnames(LP_CD4_rep1))
length(colnames(mat))
length(joint_bcs)
##################################

LP_CD4_rep1 <- LP_CD4_rep1[,joint_bcs]
hto_LP_CD4_rep1<- as.matrix(mat)
hto_LP_CD4_rep1<-hto_LP_CD4_rep1[,joint_bcs]



hto_LP_CD4_rep1<-hto_LP_CD4_rep1[1:6,]#


rownames(hto_LP_CD4_rep1)<-c("HTO_1","HTO_2","HTO_3","HTO_4","HTO_5","HTO_6")

LP_CD4_rep1 <- CreateSeuratObject(counts = LP_CD4_rep1, min.cells = 100, min.features = 200)


LP_CD4_rep1[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep1, pattern = "^mt-")


LP_CD4_rep1_sc <- as.SingleCellExperiment(LP_CD4_rep1)
LP_CD4_rep1_sc <- scDblFinder(LP_CD4_rep1_sc)

ind_doublets_LP_CD4_rep1<-which(LP_CD4_rep1_sc$scDblFinder.class=="doublet")
LP_CD4_rep1_sc <- LP_CD4_rep1_sc[,-ind_doublets_LP_CD4_rep1]

high_mitoexp_LP_CD4_rep1<-isOutlier(LP_CD4_rep1_sc@colData@listData$percent.mt, nmads=4, type="higher")
high_ncount_LP_CD4_rep1<-isOutlier(LP_CD4_rep1_sc@colData@listData$nCount_RNA, nmads=2.5, type="higher")
low_ncount_LP_CD4_rep1<-isOutlier(LP_CD4_rep1_sc@colData@listData$nCount_RNA, nmads=3, type="lower")
high_nfeature_LP_CD4_rep1<-isOutlier(LP_CD4_rep1_sc@colData@listData$nFeature_RNA, nmads=2.5, type="higher")
low_nfeature_LP_CD4_rep1<-isOutlier(LP_CD4_rep1_sc@colData@listData$nFeature_RNA, nmads=3, type="lower")

ind_outliers_LP_CD4_rep1<-Reduce("|", list(high_mitoexp_LP_CD4_rep1,high_ncount_LP_CD4_rep1,low_ncount_LP_CD4_rep1,high_nfeature_LP_CD4_rep1,low_nfeature_LP_CD4_rep1))
table(ind_outliers_LP_CD4_rep1)


ind_outliers_LP_CD4_rep1<-which(ind_outliers_LP_CD4_rep1=="TRUE")


LP_CD4_rep1_sc <- LP_CD4_rep1_sc[,-ind_outliers_LP_CD4_rep1]

LP_CD4_rep1 <- subset(LP_CD4_rep1, cells=LP_CD4_rep1_sc@colData@rownames)


# Normalize RNA data with log normalization
LP_CD4_rep1 <- NormalizeData(LP_CD4_rep1)
# Find and scale variable features
LP_CD4_rep1 <- FindVariableFeatures(LP_CD4_rep1,selection.method = "vst",mean.function = ExpMean)
LP_CD4_rep1 <- ScaleData(LP_CD4_rep1, features = rownames(LP_CD4_rep1))#

# Add HTO data as a new assay independent from RNA
LP_CD4_rep1[["HTO"]] <- CreateAssayObject(counts = hto_LP_CD4_rep1[, colnames(x = LP_CD4_rep1)])
# Normalize HTO data,here we use centered log-ratio (CLR) transformation
LP_CD4_rep1 <- NormalizeData(LP_CD4_rep1, assay = "HTO", normalization.method = "CLR" )
#"CLR" "LogNormalize"
LP_CD4_rep1_demux <- HTODemux(LP_CD4_rep1, assay = "HTO", positive.quantile = 0.95)

Idents(LP_CD4_rep1_demux) <- "HTO_classification.global"


LP_CD4_rep1_demux.subset <- LP_CD4_rep1_demux

DefaultAssay(LP_CD4_rep1_demux.subset) <- "HTO"
LP_CD4_rep1_demux.subset <- ScaleData(LP_CD4_rep1_demux.subset, features = rownames(LP_CD4_rep1_demux.subset),
    verbose = FALSE)
LP_CD4_rep1_demux.subset <- RunPCA(LP_CD4_rep1_demux.subset, features = rownames(LP_CD4_rep1_demux.subset), approx = FALSE)


# Extract the singlets
LP_CD4_rep1.singlet <- subset(LP_CD4_rep1_demux, idents = "Singlet")
LP_CD4_rep1.singlet[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep1.singlet, pattern = "^mt-")


LP_CD4_rep1_filtered<-LP_CD4_rep1.singlet


LP_CD4_rep1_filtered@meta.data$orig.ident<-as.character(LP_CD4_rep1_filtered@meta.data$orig.ident)

ind_WT<-which(LP_CD4_rep1_filtered@meta.data$HTO_maxID %in% c("HTO-1","HTO-3","HTO-5"))
ind_Wmutant<-which(LP_CD4_rep1_filtered@meta.data$HTO_maxID %in% c("HTO-2","HTO-4","HTO-6"))

LP_CD4_rep1_filtered@meta.data$orig.ident[ind_WT]<-"WT_3wks"
LP_CD4_rep1_filtered@meta.data$orig.ident[ind_Wmutant]<-"WW_3wks"


Idents(LP_CD4_rep1_filtered) <- "orig.ident"


saveRDS(LP_CD4_rep1_filtered,"LP_CD4_rep1_filtered.rds")
