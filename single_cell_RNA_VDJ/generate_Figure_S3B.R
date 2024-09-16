

####################
####################
#setwd("current_directory")

library(dplyr)
library(Seurat)
library(patchwork)
library(scCustomize)

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


gene_list_new<-read_excel("LP_dotplot_genes_cell_type_specific_markers.xlsx",sheet= "Sheet1")


features_selec<-gene_list_new$Genes



#############
#############
p_DotPlot_LP_CD4.integrated_selected_markers<-DotPlot_scCustom(seurat_object = LP_CD4.integrated, features = features_selec, x_lab_rotate = TRUE)

pdf(file = "Figure_S3B.pdf",width=9, height=6)
print(p_DotPlot_LP_CD4.integrated_selected_markers)
dev.off()
#############
#############
