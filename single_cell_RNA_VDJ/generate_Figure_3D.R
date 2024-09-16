####################
####################
#setwd("current_directory")

library(dplyr)
library(Seurat)
library(patchwork)

library(scales)
library(ggplot2)

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




###########
###########
df_cloneType_props_per_orig.ident.new<-LP_CD4.integrated@meta.data[-which(is.na(LP_CD4.integrated@meta.data$cloneType==TRUE)),] %>%
  group_by(orig.ident.new, cloneType) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n))

df_cloneType_props_per_orig.ident.new<-df_cloneType_props_per_orig.ident.new[order(-df_cloneType_props_per_orig.ident.new$proportion),]
###########
###########


sum(df_cloneType_props_per_orig.ident.new$proportion[which(df_cloneType_props_per_orig.ident.new$orig.ident.new=="WT_3wks")])

###########
###########

###########
###########
p_df_cloneType_props_per_orig.ident.new<-ggplot(df_cloneType_props_per_orig.ident.new, aes(x = orig.ident.new, y = proportion, fill = cloneType))+geom_bar(stat = "identity")+scale_fill_manual(values = rev(c("#29285F","#613077","#AC4B6B","#DE9155","#E6E45D")))+ labs(x = "",y="Fraction of cells") 
###########
###########
pdf(file = "Figure_3D.pdf",width=9, height=6)
print(p_df_cloneType_props_per_orig.ident.new)
dev.off()
###########
###########