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


matrix_dir_rep5="./HTO_counts_rep5/"
barcode.path <- paste0(matrix_dir_rep5, "barcodes.tsv")
features.path <- paste0(matrix_dir_rep5, "features.tsv")
matrix.path <- paste0(matrix_dir_rep5, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1


LP_CD4_rep5=Read10X_h5("./RNA_counts_rep5/raw_feature_bc_matrix.h5")



colnames(mat) <- paste(colnames(mat),"-1", sep ="")

head(colnames(LP_CD4_rep5))


head(colnames(mat))


LP_CD4_rep5<-LP_CD4_rep5[,order(colnames(LP_CD4_rep5))]
mat<-mat[,order(colnames(mat))]

joint_bcs <- intersect(colnames(LP_CD4_rep5),colnames(mat))

length(colnames(LP_CD4_rep5))
length(colnames(mat))
length(joint_bcs)
##################################

LP_CD4_rep5 <- LP_CD4_rep5[,joint_bcs]
hto_LP_CD4_rep5<- as.matrix(mat)
hto_LP_CD4_rep5<-hto_LP_CD4_rep5[,joint_bcs]

rownames(hto_LP_CD4_rep5)



hto_LP_CD4_rep5<-hto_LP_CD4_rep5[1:5,]#

rownames(hto_LP_CD4_rep5)<-c("HTO_1","HTO_2","HTO_3","HTO_4","HTO_5")

LP_CD4_rep5 <- CreateSeuratObject(counts = LP_CD4_rep5, min.cells = 100, min.features = 200)

LP_CD4_rep5[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep5, pattern = "^mt-")


LP_CD4_rep5_sc <- as.SingleCellExperiment(LP_CD4_rep5)
LP_CD4_rep5_sc <- scDblFinder(LP_CD4_rep5_sc)

ind_doublets_LP_CD4_rep5<-which(LP_CD4_rep5_sc$scDblFinder.class=="doublet")
LP_CD4_rep5_sc <- LP_CD4_rep5_sc[,-ind_doublets_LP_CD4_rep5]

high_mitoexp_LP_CD4_rep5<-isOutlier(LP_CD4_rep5_sc@colData@listData$percent.mt, nmads=4, type="higher")
high_ncount_LP_CD4_rep5<-isOutlier(LP_CD4_rep5_sc@colData@listData$nCount_RNA, nmads=2.5, type="higher")
low_ncount_LP_CD4_rep5<-isOutlier(LP_CD4_rep5_sc@colData@listData$nCount_RNA, nmads=3, type="lower")
high_nfeature_LP_CD4_rep5<-isOutlier(LP_CD4_rep5_sc@colData@listData$nFeature_RNA, nmads=2.5, type="higher")
low_nfeature_LP_CD4_rep5<-isOutlier(LP_CD4_rep5_sc@colData@listData$nFeature_RNA, nmads=3, type="lower")

ind_outliers_LP_CD4_rep5<-Reduce("|", list(high_mitoexp_LP_CD4_rep5,high_ncount_LP_CD4_rep5,low_ncount_LP_CD4_rep5,high_nfeature_LP_CD4_rep5,low_nfeature_LP_CD4_rep5))
table(ind_outliers_LP_CD4_rep5)


ind_outliers_LP_CD4_rep5<-which(ind_outliers_LP_CD4_rep5=="TRUE")


LP_CD4_rep5_sc <- LP_CD4_rep5_sc[,-ind_outliers_LP_CD4_rep5]

LP_CD4_rep5 <- subset(LP_CD4_rep5, cells=LP_CD4_rep5_sc@colData@rownames)


# Normalize RNA data with log normalization
LP_CD4_rep5 <- NormalizeData(LP_CD4_rep5)
# Find and scale variable features
LP_CD4_rep5 <- FindVariableFeatures(LP_CD4_rep5,selection.method = "vst",mean.function = ExpMean)
LP_CD4_rep5 <- ScaleData(LP_CD4_rep5, features = rownames(LP_CD4_rep5))

# Add HTO data as a new assay independent from RNA
LP_CD4_rep5[["HTO"]] <- CreateAssayObject(counts = hto_LP_CD4_rep5[, colnames(x = LP_CD4_rep5)])
# Normalize HTO data,here we use centered log-ratio (CLR) transformation
LP_CD4_rep5 <- NormalizeData(LP_CD4_rep5, assay = "HTO", normalization.method = "CLR" )
#"CLR" "LogNormalize"
LP_CD4_rep5_demux <- HTODemux(LP_CD4_rep5, assay = "HTO", positive.quantile = 0.95)

Idents(LP_CD4_rep5_demux) <- "HTO_classification.global"


# First, we will remove negative cells from the object
LP_CD4_rep5_demux.subset <- subset(LP_CD4_rep5_demux, idents = "Negative", invert = TRUE)

DefaultAssay(LP_CD4_rep5_demux.subset) <- "HTO"
LP_CD4_rep5_demux.subset <- ScaleData(LP_CD4_rep5_demux.subset, features = rownames(LP_CD4_rep5_demux.subset),
    verbose = FALSE)
LP_CD4_rep5_demux.subset <- RunPCA(LP_CD4_rep5_demux.subset, features = rownames(LP_CD4_rep5_demux.subset), approx = FALSE)


# Extract the singlets
LP_CD4_rep5.singlet <- subset(LP_CD4_rep5_demux, idents = "Singlet")
LP_CD4_rep5.singlet[["percent.mt"]] <- PercentageFeatureSet(LP_CD4_rep5.singlet, pattern = "^mt-")

LP_CD4_rep5_filtered<-LP_CD4_rep5.singlet

LP_CD4_rep5_filtered@meta.data$orig.ident<-as.character(LP_CD4_rep5_filtered@meta.data$orig.ident)

ind_Wmutant_ANMV<-which(LP_CD4_rep5_filtered@meta.data$HTO_maxID %in% c("HTO-1","HTO-2","HTO-3"))
ind_Wmutant<-which(LP_CD4_rep5_filtered@meta.data$HTO_maxID %in% c("HTO-4","HTO-5"))

LP_CD4_rep5_filtered@meta.data$orig.ident[ind_Wmutant_ANMV]<-"WW_ANMV_12wks"
LP_CD4_rep5_filtered@meta.data$orig.ident[ind_Wmutant]<-"WW_CTRL_12wks"


Idents(LP_CD4_rep5_filtered) <- "orig.ident"

saveRDS(LP_CD4_rep5_filtered,"LP_CD4_rep5_filtered.rds")
