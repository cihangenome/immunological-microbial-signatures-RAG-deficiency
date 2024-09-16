############
############
#setwd("current_directory")

library(readxl)
library(ggplot2)

############
############
df_dist_TRA_LP_Treg_cluster_concise<-read_excel("LP_CD4_Treg_cluster_clonal_overlap_TRA_Ikzf2_Rorc_Il10.xlsx",sheet= "Overlap_Jaccard_index")
############
############
pops_to_include<-data.frame(Sample=c("LP_Treg_cluster_WW_12wks_Ikzf2+_Rorc-_Il10-",
"LP_Treg_cluster_WW_ANMV_12wks_Ikzf2+_Rorc-_Il10-",
"LP_Treg_cluster_WT_3wks_Ikzf2+_Rorc-_Il10-",
"LP_Treg_cluster_WT_12wks_Ikzf2+_Rorc-_Il10-",
"LP_Treg_cluster_WT_12wks_Ikzf2-_Rorc+_Il10+",
"LP_Treg_cluster_WW_3wks_Ikzf2+_Rorc-_Il10-",
"LP_Treg_cluster_WW_ANMV_12wks_Ikzf2-_Rorc+_Il10+",
"LP_Treg_cluster_WT_3wks_Ikzf2-_Rorc+_Il10+",
"LP_Treg_cluster_WW_3wks_Ikzf2-_Rorc+_Il10+",
"LP_Treg_cluster_WW_12wks_Ikzf2-_Rorc+_Il10+"))



########################
########################

ind_relev1<-which(df_dist_TRA_LP_Treg_cluster_concise$sample1 %in% pops_to_include$Sample)

ind_relev2<-which(df_dist_TRA_LP_Treg_cluster_concise$sample2 %in% pops_to_include$Sample)

df_dist_TRA_LP_Treg_cluster_concise_12and3wks<-df_dist_TRA_LP_Treg_cluster_concise[(intersect(ind_relev1,ind_relev2)),]
##################
##################

############
############
colnames(df_dist_TRA_LP_Treg_cluster_concise_12and3wks)[3]<-"TRA_Overlap"

df_dist_TRA_LP_Treg_cluster_concise_12and3wks_add<-df_dist_TRA_LP_Treg_cluster_concise_12and3wks

dum1<-df_dist_TRA_LP_Treg_cluster_concise_12and3wks_add$sample1
dum2<-df_dist_TRA_LP_Treg_cluster_concise_12and3wks_add$sample2

df_dist_TRA_LP_Treg_cluster_concise_12and3wks_add$sample1<-dum2
df_dist_TRA_LP_Treg_cluster_concise_12and3wks_add$sample2<-dum1

df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo<-rbind(df_dist_TRA_LP_Treg_cluster_concise_12and3wks,df_dist_TRA_LP_Treg_cluster_concise_12and3wks_add)

target_order<-sort(unique(as.character(df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo$sample1)))



df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo$sample1<-factor(df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo$sample1,levels=target_order)

df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo$sample2<-factor(df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo$sample2,levels=target_order)

max_val<-max(df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo$TRA_Overlap)

pdf("Figure_S5D.pdf",width=8, height=6)
ggplot(df_dist_TRA_LP_Treg_cluster_concise_12and3wks_combo,aes(x=sample1,y=sample2,fill=TRA_Overlap))+geom_tile()+ ylab("")+ xlab("")+ scale_fill_continuous(limits=c(0,max_val), breaks=seq(0,max_val,by=0.04))+scale_fill_distiller(palette="RdYlBu")+ theme_bw()+
theme(strip.background =element_rect(fill="white"))+theme(axis.text.x = element_text(angle = 90,vjust = 1, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
############
############
