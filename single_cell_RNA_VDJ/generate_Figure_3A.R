####################
####################
#setwd("current_directory")

library(dplyr)
library(Seurat)
library(patchwork)

LP_CD4.integrated<-readRDS("sc_LP_CD4.rds")

###############
###############
df_LP_metadata<-LP_CD4.integrated@meta.data


df_cell_sub_type_seurat_clusters<-data.frame(seurat_clusters=LP_CD4.integrated@meta.data$seurat_clusters)

df_cell_sub_type_seurat_clusters$seurat_clusters_new<-df_cell_sub_type_seurat_clusters$seurat_clusters

mapping_cell_sub_type_seurat_clusters<-data.frame(simple_new=c("Treg","Naive","IFN-gamma+","IL-17+","IFN-gamma response (Gpb4+)", "Activated T", "Tfh", "IFN-gamma response", "APC", "Th2", "Proliferative T", "Memory T"),original= sort(unique(LP_CD4.integrated@meta.data$seurat_clusters)))


df_cell_sub_type_seurat_clusters$seurat_clusters_new <- dplyr::recode(
  df_cell_sub_type_seurat_clusters$seurat_clusters_new,
  !!!setNames(mapping_cell_sub_type_seurat_clusters$simple_new, mapping_cell_sub_type_seurat_clusters$original)
)


LP_CD4.integrated@meta.data$seurat_clusters_new<-df_cell_sub_type_seurat_clusters$seurat_clusters_new
####################
####################

Idents(LP_CD4.integrated)<-"seurat_clusters_new"

LP_CD4.integrated<-subset(LP_CD4.integrated,idents=c("Memory T"),invert=TRUE)

#unique(LP_CD4.integrated@meta.data[["orig.ident"]])
#[1] "WW_ANMV_12wks"

LP_CD4.integrated$orig.ident.new<-LP_CD4.integrated$orig.ident
LP_CD4.integrated$orig.ident.new<-as.character(LP_CD4.integrated$orig.ident.new)

LP_CD4.integrated$orig.ident.new[grep("ANMV",LP_CD4.integrated$orig.ident.new)]<-"WW_12wks"


#####################
#####################
LP_CD4.integrated$orig.ident.new<-as.factor(LP_CD4.integrated$orig.ident.new)
LP_CD4.integrated$orig.ident.new<-factor(LP_CD4.integrated$orig.ident.new,levels=c("WT_3wks","WW_3wks","WT_12wks","WW_12wks"))
#####################
#####################

##################
##################
Idents(object = LP_CD4.integrated) <- "seurat_clusters_new"
pdf(file = "Figure_3A.pdf",width=9, height=6)
DimPlot(LP_CD4.integrated,split.by="orig.ident.new",reduction = "umap", label=FALSE, label.size = 5,pt.size=0.5,ncol=2)
dev.off()
##################
##################
