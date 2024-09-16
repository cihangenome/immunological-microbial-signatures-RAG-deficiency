#setwd("current_directory")

library(ggplot2)
library(Seurat)

sc_LP_CD4<-readRDS("sc_LP_CD4.rds")
sc_Treg_WT_WW_12wks<-readRDS("sc_Treg_WT_WW_12wks.rds")

Idents(sc_LP_CD4)<-"seurat_clusters"

#subsetting for Tregs
sc_Treg_LP_CD4<-subset(sc_LP_CD4,idents=c("0"))

#projecting the Treg ANMV on the UMAP of the 12wks WT+WW cells

Idents(sc_Treg_LP_CD4)<-"orig.ident"
sc_Treg_WW_ANMV_12wks<-subset(sc_Treg_LP_CD4,idents=c("WW_ANMV_12wks"))

#sc_Treg_WW_ANMV_12wks$orig.ident<-as.factor(sc_Treg_WW_ANMV_12wks$orig.ident)

levels(sc_Treg_WT_WW_12wks$orig.ident)
sc_Treg_WT_WW_12wks$orig.ident<-as.factor(sc_Treg_WT_WW_12wks$orig.ident)

Idents(sc_Treg_WT_WW_12wks) <- "RNA_snn_res.0.4"

sc_Treg_LP_CD4.anchors <- FindTransferAnchors(reference = sc_Treg_WT_WW_12wks, query = sc_Treg_WW_ANMV_12wks,dims = 1:30, reference.reduction = "pca")

predictions <- TransferData(anchorset = sc_Treg_LP_CD4.anchors, refdata = sc_Treg_WT_WW_12wks$RNA_snn_res.0.4,dims = 1:30)

sc_Treg_WW_ANMV_12wks <- AddMetaData(sc_Treg_WW_ANMV_12wks, metadata = predictions)

sc_Treg_WT_WW_12wks <- RunUMAP(sc_Treg_WT_WW_12wks, dims = 1:30, reduction = "pca", return.model = TRUE)

sc_Treg_WW_ANMV_12wks <- MapQuery(anchorset = sc_Treg_LP_CD4.anchors, reference = sc_Treg_WT_WW_12wks, query = sc_Treg_WW_ANMV_12wks,refdata = sc_Treg_WT_WW_12wks$RNA_snn_res.0.4, reference.reduction = "pca", reduction.model = "umap")

Idents(sc_Treg_WW_ANMV_12wks) <- "predicted.id"

pdf(file = "Figure_S5C.pdf",width=6, height=4)
DimPlot(sc_Treg_WW_ANMV_12wks, reduction = "ref.umap", label=TRUE, label.size = 5,pt.size=0.5)+ggplot2::labs(title="LP CD4+ WW_ANMV_12wks projected onto reclustered Tregs in WT_WW_12wks \n All cells are Tregs in the original cluster 0 with VDJ data (cluster IDs: RNA_snn_res.0.4))")+ theme(plot.title = element_text(size=8))
dev.off()