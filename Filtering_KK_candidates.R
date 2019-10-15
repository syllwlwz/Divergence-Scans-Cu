options(java.parameters = "-Xmx12g")
require(xlsx)

SNPeffKK<-read.table("SNPeff_KromKosi.ann.HIGH.table",sep="\t",header=F)
colnames(SNPeffKK)<-c("Scaffold","Pos","Base","Effect","Impact","Gene","Warning")
Gene<-gsub("^.*-","",SNPeffKK$Gene)
SNPeffKK$Gene<-Gene
SNPeffKKexcl<-SNPeffKK[!SNPeffKK$Warning=="",]

#KromKosi directory
KK_FST<-read.csv("FST_KromKosi_01percent.csv",header=T)
KK_DXY<-read.csv("DXY_KromKosi_01percent.csv",header=T)
KK_AFDabs<-read.csv("AFDabs_KromKosi_01percent.csv",header=T)
KK_NIELSEN<-read.csv("NIELSEN_KromKosi_01percent.csv",header=T)
KK_FLK<-read.csv("FLK_KromKosi_01percent.csv",header=T)
KK_DD<-read.csv("DD_KromKosi_01percent.csv",header=T)
KK_VARLD<-read.csv("VarLD_KromKosi_01percent.csv",header=T)
KK_HighSNP<-read.table("GenesHIGHeffect_above40percentAFdiff_withUG.txt",sep="\t",header=T,fill=T,quote="")
KK_HighSNP<-KK_HighSNP[KK_HighSNP$UG_AF_diff!="NO",]
KK_HighSNP<-KK_HighSNP[!KK_HighSNP$Gene%in%SNPeffKKexcl$Gene,]

#change to GROM_analysis directory
KK_INDEL<-read.xlsx2("F://Cu_project//CNVS_Indels//Indel_HC_GROM.xlsx",1,header=T)

KK_FSTr<-cbind(KK_FST,rep("FST",nrow(KK_FST)))
KK_DXYr<-cbind(KK_DXY,rep("DXY",nrow(KK_DXY)))
KK_AFDabsr<-cbind(KK_AFDabs,rep("AFDabs",nrow(KK_AFDabs)))
KK_NIELSENr<-cbind(KK_NIELSEN,rep("NIELSEN",nrow(KK_NIELSEN)))
KK_FLKr<-cbind(KK_FLK,rep("FLK",nrow(KK_FLK)))
KK_DDr<-cbind(KK_DD,rep("DD",nrow(KK_DD)))
KK_VARLDr<-cbind(KK_VARLD,rep("VARLD",nrow(KK_VARLD)))

names(KK_FSTr)[99]<-"Outlier_metric"
names(KK_DXYr)[99]<-"Outlier_metric"
names(KK_AFDabsr)[99]<-"Outlier_metric"
names(KK_NIELSENr)[99]<-"Outlier_metric"
names(KK_FLKr)[99]<-"Outlier_metric"
names(KK_DDr)[99]<-"Outlier_metric"
names(KK_VARLDr)[99]<-"Outlier_metric"

names(KK_FSTr)[93]<-"Rank_in_outlier_metric"
names(KK_DXYr)[93]<-"Rank_in_outlier_metric"
names(KK_AFDabsr)[93]<-"Rank_in_outlier_metric"
names(KK_NIELSENr)[93]<-"Rank_in_outlier_metric"
names(KK_FLKr)[93]<-"Rank_in_outlier_metric"
names(KK_DDr)[93]<-"Rank_in_outlier_metric"
names(KK_VARLDr)[93]<-"Rank_in_outlier_metric"

names(KK_FLKr)[25]<-"Flk_pvalue"
KK<-rbind(KK_FSTr,KK_DXYr,KK_AFDabsr,KK_NIELSENr,KK_FLKr,KK_DDr,KK_VARLDr)
#844 halleri genes

KKFiltered<-KK[(KK$TD_diff_nearest_Krom_Kosi_stand<0|KK$Sweed_Krom_nearest_Likelihood>KK$Sweed_Kosi_nearest_Likelihood)&(KK$TajimasD_Krom_nearest_stand<0|KK$Sweed_Krom_nearest_Likelihood>0),]

KKFiltered2<-KK[KK$TajimasD_Krom_nearest_stand<0|KK$Sweed_Krom_nearest_Likelihood>0,]

KKFiltered3<-KK[KK$TajimasD_Krom_nearest_stand<0|KK$Sweed_Krom_nearest_Likelihood>median(na.omit(KK$Sweed_Krom_nearest_Likelihood)),]

KK_HighSNP_fixed<-KK_HighSNP[KK_HighSNP$AF_Krom<=0.1|KK_HighSNP$AF_Krom>=0.9,]

Krom_cn<-read.xlsx2("F://Cu_project//CNVS_Indels//CNVs_overlap_strict.xlsx",1,header=T)

colnames(KK_HighSNP_fixed)[1:3]<-c("Ah_ID","Scaffold","Pos")



#write.xlsx2(KKFiltered,"Genes_KromKosi.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE)
#write.xlsx2(KKFiltered2,"Genes_KromKosi.xlsx",sheetName="Less_strict_genome_scans",col.names=TRUE,row.names=FALSE,append=T)
#write.xlsx2(KK_HighSNP_fixed,"Genes_KromKosi.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
#write.xlsx2(KK_INDEL,"Genes_KromKosi.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T)
#write.xlsx2(Krom_cn,"Genes_KromKosi.xlsx",sheetName="CopyNumberVariations",col.names=TRUE,row.names=FALSE,append=T)

write.xlsx2(KKFiltered3,"Genes_KromKosi_wVarLD.xlsx",sheetName="Classical_genome_scans",col.names=TRUE,row.names=FALSE,quote=F)
write.xlsx2(KK_HighSNP_fixed,"Genes_KromKosi_wVarLD.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T,quote=F)
write.xlsx2(KK_INDEL,"Genes_KromKosi_wVarLD.xlsx",sheetName="Indels",col.names=TRUE,row.names=FALSE,append=T,quote=F)
write.xlsx2(Krom_cn,"Genes_KromKosi_wVarLD.xlsx",sheetName="CopyNumberVariations",col.names=TRUE,row.names=FALSE,append=T,quote=F)



