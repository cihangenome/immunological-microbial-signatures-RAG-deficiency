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

idents_selec<-grep("ANMV",unique(df_LP_metadata$orig.ident), value = TRUE)

Idents(LP_CD4.integrated)<-"orig.ident"

LP_CD4.integrated<-subset(LP_CD4.integrated,idents=idents_selec,invert=FALSE)

Idents(LP_CD4.integrated)<-"seurat_clusters_new"

LP_CD4.integrated<-subset(LP_CD4.integrated,idents=c("Memory T"),invert=TRUE)

#unique(LP_CD4.integrated@meta.data[["orig.ident"]])
#[1] "WW_ANMV_12wks"

##################
##################
Idents(object = LP_CD4.integrated) <- "seurat_clusters_new"
pdf(file = "Figure_5E.pdf",width=9, height=6)
DimPlot(LP_CD4.integrated, reduction = "umap", label=FALSE, label.size = 5,pt.size=0.5)
dev.off()
##################
##################