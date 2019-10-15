############################################################
#GROM data#
############################################################

Krom_Indels<-read.table("Krom.indels.vcf",sep="\t",header=F)
Krom_CNVs1<-read.table("Krom.cnvs.vcf",sep="\t",header=F)
Kosi_Indels<-read.table("Kosi.indels.vcf",sep="\t",header=F)
Kosi_CNVs1<-read.table("Kosi.cnvs.vcf",sep="\t",header=F)
Lang_Indels<-read.table("Lang.indels.vcf",sep="\t",header=F)
Lang_CNVs1<-read.table("Lang.cnvs.vcf",sep="\t",header=F)
Best_Indels<-read.table("Best.indels.vcf",sep="\t",header=F)
Best_CNVs1<-read.table("Best.cnvs.vcf",sep="\t",header=F)
Noss_Indels<-read.table("Noss.indels.vcf",sep="\t",header=F)
Noss_CNVs1<-read.table("Noss.cnvs.vcf",sep="\t",header=F)
Pais_Indels<-read.table("Pais.indels.vcf",sep="\t",header=F)
Pais_CNVs1<-read.table("Pais.cnvs.vcf",sep="\t",header=F)

names(Krom_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Krom_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Kosi_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Kosi_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Lang_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Lang_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Best_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Best_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Noss_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Noss_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Pais_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Pais_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Noss_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Noss_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Pais_Indels)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")
names(Pais_CNVs1)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT_names","FORMAT")

#GROM filtering: SRD>=24 SPR==0 SPR: probability of start breakpoint evidence ocurring by chance; SRD: physical read depth at start breakpoint 
#need to split into insertions and deletions
Krom_Indels1_filter<-Krom_Indels[Krom_Indels$FORMAT_names=="SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP",]
Krom_Indels1_filter1<-strsplit(as.character(Krom_Indels1_filter$FORMAT),":")
Krom_Indels1_filter2<-lapply(Krom_Indels1_filter1,"[",1:max(unlist(lapply(Krom_Indels1_filter1,length))))
Krom_Indels1_filter3<-data.frame(matrix(unlist(Krom_Indels1_filter2),nrow=length(Krom_Indels1_filter2),byrow=T))
Krom_Indels1<-cbind(Krom_Indels1_filter,Krom_Indels1_filter3)
colnames(Krom_Indels1)[11:19]<-c("SPR","SEV","SRD","SCO","ECO","SOT","EOT","SSC","HP")
Krom_Indelsa<-Krom_Indels1[(as.numeric(as.character(Krom_Indels1$SRD))>=24)&(round(as.numeric(as.character(Krom_Indels1$SPR)),digits=0)==0),]

Krom_Indels2_filter<-Krom_Indels[Krom_Indels$FORMAT_names=="SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP",]
Krom_Indels2_filter1<-strsplit(as.character(Krom_Indels2_filter$FORMAT),":")
Krom_Indels2_filter2<-lapply(Krom_Indels2_filter1,"[",1:max(unlist(lapply(Krom_Indels2_filter1,length))))
Krom_Indels2_filter3<-data.frame(matrix(unlist(Krom_Indels2_filter2),nrow=length(Krom_Indels2_filter2),byrow=T))
Krom_Indels2<-cbind(Krom_Indels2_filter,Krom_Indels2_filter3)
colnames(Krom_Indels2)[11:23]<-c("SPR","EPR","SEV","EEV","SRD","ERD","SCO","ECO","SOT","EOT","SSC","ESC","HP")
Krom_Indelsb<-Krom_Indels2[(as.numeric(as.character(Krom_Indels2$SRD))>=24)&(round(as.numeric(as.character(Krom_Indels2$SPR)),digits=0)==0),]

Krom_Indels<-rbind(Krom_Indelsa[,1:10],Krom_Indelsb[,1:10])

Kosi_Indels1_filter<-Kosi_Indels[Kosi_Indels$FORMAT_names=="SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP",]
Kosi_Indels1_filter1<-strsplit(as.character(Kosi_Indels1_filter$FORMAT),":")
Kosi_Indels1_filter2<-lapply(Kosi_Indels1_filter1,"[",1:max(unlist(lapply(Kosi_Indels1_filter1,length))))
Kosi_Indels1_filter3<-data.frame(matrix(unlist(Kosi_Indels1_filter2),nrow=length(Kosi_Indels1_filter2),byrow=T))
Kosi_Indels1<-cbind(Kosi_Indels1_filter,Kosi_Indels1_filter3)
colnames(Kosi_Indels1)[11:19]<-c("SPR","SEV","SRD","SCO","ECO","SOT","EOT","SSC","HP")
Kosi_Indelsa<-Kosi_Indels1[(as.numeric(as.character(Kosi_Indels1$SRD))>=24)&(round(as.numeric(as.character(Kosi_Indels1$SPR)),digits=0)==0),]

Kosi_Indels2_filter<-Kosi_Indels[Kosi_Indels$FORMAT_names=="SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP",]
Kosi_Indels2_filter1<-strsplit(as.character(Kosi_Indels2_filter$FORMAT),":")
Kosi_Indels2_filter2<-lapply(Kosi_Indels2_filter1,"[",1:max(unlist(lapply(Kosi_Indels2_filter1,length))))
Kosi_Indels2_filter3<-data.frame(matrix(unlist(Kosi_Indels2_filter2),nrow=length(Kosi_Indels2_filter2),byrow=T))
Kosi_Indels2<-cbind(Kosi_Indels2_filter,Kosi_Indels2_filter3)
colnames(Kosi_Indels2)[11:23]<-c("SPR","EPR","SEV","EEV","SRD","ERD","SCO","ECO","SOT","EOT","SSC","ESC","HP")
Kosi_Indelsb<-Kosi_Indels2[(as.numeric(as.character(Kosi_Indels2$SRD))>=24)&(round(as.numeric(as.character(Kosi_Indels2$SPR)),digits=0)==0),]

Kosi_Indels<-rbind(Kosi_Indelsa[,1:10],Kosi_Indelsb[,1:10])

Lang_Indels1_filter<-Lang_Indels[Lang_Indels$FORMAT_names=="SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP",]
Lang_Indels1_filter1<-strsplit(as.character(Lang_Indels1_filter$FORMAT),":")
Lang_Indels1_filter2<-lapply(Lang_Indels1_filter1,"[",1:max(unlist(lapply(Lang_Indels1_filter1,length))))
Lang_Indels1_filter3<-data.frame(matrix(unlist(Lang_Indels1_filter2),nrow=length(Lang_Indels1_filter2),byrow=T))
Lang_Indels1<-cbind(Lang_Indels1_filter,Lang_Indels1_filter3)
colnames(Lang_Indels1)[11:19]<-c("SPR","SEV","SRD","SCO","ECO","SOT","EOT","SSC","HP")
Lang_Indelsa<-Lang_Indels1[(as.numeric(as.character(Lang_Indels1$SRD))>=24)&(round(as.numeric(as.character(Lang_Indels1$SPR)),digits=0)==0),]

Lang_Indels2_filter<-Lang_Indels[Lang_Indels$FORMAT_names=="SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP",]
Lang_Indels2_filter1<-strsplit(as.character(Lang_Indels2_filter$FORMAT),":")
Lang_Indels2_filter2<-lapply(Lang_Indels2_filter1,"[",1:max(unlist(lapply(Lang_Indels2_filter1,length))))
Lang_Indels2_filter3<-data.frame(matrix(unlist(Lang_Indels2_filter2),nrow=length(Lang_Indels2_filter2),byrow=T))
Lang_Indels2<-cbind(Lang_Indels2_filter,Lang_Indels2_filter3)
colnames(Lang_Indels2)[11:23]<-c("SPR","EPR","SEV","EEV","SRD","ERD","SCO","ECO","SOT","EOT","SSC","ESC","HP")
Lang_Indelsb<-Lang_Indels2[(as.numeric(as.character(Lang_Indels2$SRD))>=24)&(round(as.numeric(as.character(Lang_Indels2$SPR)),digits=0)==0),]

Lang_Indels<-rbind(Lang_Indelsa[,1:10],Lang_Indelsb[,1:10])

Best_Indels1_filter<-Best_Indels[Best_Indels$FORMAT_names=="SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP",]
Best_Indels1_filter1<-strsplit(as.character(Best_Indels1_filter$FORMAT),":")
Best_Indels1_filter2<-lapply(Best_Indels1_filter1,"[",1:max(unlist(lapply(Best_Indels1_filter1,length))))
Best_Indels1_filter3<-data.frame(matrix(unlist(Best_Indels1_filter2),nrow=length(Best_Indels1_filter2),byrow=T))
Best_Indels1<-cbind(Best_Indels1_filter,Best_Indels1_filter3)
colnames(Best_Indels1)[11:19]<-c("SPR","SEV","SRD","SCO","ECO","SOT","EOT","SSC","HP")
Best_Indelsa<-Best_Indels1[(as.numeric(as.character(Best_Indels1$SRD))>=24)&(round(as.numeric(as.character(Best_Indels1$SPR)),digits=0)==0),]

Best_Indels2_filter<-Best_Indels[Best_Indels$FORMAT_names=="SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP",]
Best_Indels2_filter1<-strsplit(as.character(Best_Indels2_filter$FORMAT),":")
Best_Indels2_filter2<-lapply(Best_Indels2_filter1,"[",1:max(unlist(lapply(Best_Indels2_filter1,length))))
Best_Indels2_filter3<-data.frame(matrix(unlist(Best_Indels2_filter2),nrow=length(Best_Indels2_filter2),byrow=T))
Best_Indels2<-cbind(Best_Indels2_filter,Best_Indels2_filter3)
colnames(Best_Indels2)[11:23]<-c("SPR","EPR","SEV","EEV","SRD","ERD","SCO","ECO","SOT","EOT","SSC","ESC","HP")
Best_Indelsb<-Best_Indels2[(as.numeric(as.character(Best_Indels2$SRD))>=24)&(round(as.numeric(as.character(Best_Indels2$SPR)),digits=0)==0),]

Best_Indels<-rbind(Best_Indelsa[,1:10],Best_Indelsb[,1:10])

Noss_Indels1_filter<-Noss_Indels[Noss_Indels$FORMAT_names=="SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP",]
Noss_Indels1_filter1<-strsplit(as.character(Noss_Indels1_filter$FORMAT),":")
Noss_Indels1_filter2<-lapply(Noss_Indels1_filter1,"[",1:max(unlist(lapply(Noss_Indels1_filter1,length))))
Noss_Indels1_filter3<-data.frame(matrix(unlist(Noss_Indels1_filter2),nrow=length(Noss_Indels1_filter2),byrow=T))
Noss_Indels1<-cbind(Noss_Indels1_filter,Noss_Indels1_filter3)
colnames(Noss_Indels1)[11:19]<-c("SPR","SEV","SRD","SCO","ECO","SOT","EOT","SSC","HP")
Noss_Indelsa<-Noss_Indels1[(as.numeric(as.character(Noss_Indels1$SRD))>=24)&(round(as.numeric(as.character(Noss_Indels1$SPR)),digits=0)==0),]

Noss_Indels2_filter<-Noss_Indels[Noss_Indels$FORMAT_names=="SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP",]
Noss_Indels2_filter1<-strsplit(as.character(Noss_Indels2_filter$FORMAT),":")
Noss_Indels2_filter2<-lapply(Noss_Indels2_filter1,"[",1:max(unlist(lapply(Noss_Indels2_filter1,length))))
Noss_Indels2_filter3<-data.frame(matrix(unlist(Noss_Indels2_filter2),nrow=length(Noss_Indels2_filter2),byrow=T))
Noss_Indels2<-cbind(Noss_Indels2_filter,Noss_Indels2_filter3)
colnames(Noss_Indels2)[11:23]<-c("SPR","EPR","SEV","EEV","SRD","ERD","SCO","ECO","SOT","EOT","SSC","ESC","HP")
Noss_Indelsb<-Noss_Indels2[(as.numeric(as.character(Noss_Indels2$SRD))>=24)&(round(as.numeric(as.character(Noss_Indels2$SPR)),digits=0)==0),]

Noss_Indels<-rbind(Noss_Indelsa[,1:10],Noss_Indelsb[,1:10])

Pais_Indels1_filter<-Pais_Indels[Pais_Indels$FORMAT_names=="SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP",]
Pais_Indels1_filter1<-strsplit(as.character(Pais_Indels1_filter$FORMAT),":")
Pais_Indels1_filter2<-lapply(Pais_Indels1_filter1,"[",1:max(unlist(lapply(Pais_Indels1_filter1,length))))
Pais_Indels1_filter3<-data.frame(matrix(unlist(Pais_Indels1_filter2),nrow=length(Pais_Indels1_filter2),byrow=T))
Pais_Indels1<-cbind(Pais_Indels1_filter,Pais_Indels1_filter3)
colnames(Pais_Indels1)[11:19]<-c("SPR","SEV","SRD","SCO","ECO","SOT","EOT","SSC","HP")
Pais_Indelsa<-Pais_Indels1[(as.numeric(as.character(Pais_Indels1$SRD))>=24)&(round(as.numeric(as.character(Pais_Indels1$SPR)),digits=0)==0),]

Pais_Indels2_filter<-Pais_Indels[Pais_Indels$FORMAT_names=="SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP",]
Pais_Indels2_filter1<-strsplit(as.character(Pais_Indels2_filter$FORMAT),":")
Pais_Indels2_filter2<-lapply(Pais_Indels2_filter1,"[",1:max(unlist(lapply(Pais_Indels2_filter1,length))))
Pais_Indels2_filter3<-data.frame(matrix(unlist(Pais_Indels2_filter2),nrow=length(Pais_Indels2_filter2),byrow=T))
Pais_Indels2<-cbind(Pais_Indels2_filter,Pais_Indels2_filter3)
colnames(Pais_Indels2)[11:23]<-c("SPR","EPR","SEV","EEV","SRD","ERD","SCO","ECO","SOT","EOT","SSC","ESC","HP")
Pais_Indelsb<-Pais_Indels2[(as.numeric(as.character(Pais_Indels2$SRD))>=24)&(round(as.numeric(as.character(Pais_Indels2$SPR)),digits=0)==0),]

Pais_Indels<-rbind(Pais_Indelsa[,1:10],Pais_Indelsb[,1:10])

Krom_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Krom_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Krom_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Krom_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Krom_CNVs<-data.frame(Krom_CNVs1[,1:2],as.numeric(as.character(Krom_CNVs_END[,2])),Krom_CNVs1[,3:7],Krom_CNVs_info)
names(Krom_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Kosi_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Kosi_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Kosi_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Kosi_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Kosi_CNVs<-data.frame(Kosi_CNVs1[,1:2],as.numeric(as.character(Kosi_CNVs_END[,2])),Kosi_CNVs1[,3:7],Kosi_CNVs_info)
names(Kosi_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Lang_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Lang_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Lang_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Lang_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Lang_CNVs<-data.frame(Lang_CNVs1[,1:2],as.numeric(as.character(Lang_CNVs_END[,2])),Lang_CNVs1[,3:7],Lang_CNVs_info)
names(Lang_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Best_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Best_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Best_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Best_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Best_CNVs<-data.frame(Best_CNVs1[,1:2],as.numeric(as.character(Best_CNVs_END[,2])),Best_CNVs1[,3:7],Best_CNVs_info)
names(Best_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Noss_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Noss_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Noss_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Noss_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Noss_CNVs<-data.frame(Noss_CNVs1[,1:2],as.numeric(as.character(Noss_CNVs_END[,2])),Noss_CNVs1[,3:7],Noss_CNVs_info)
names(Noss_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

Pais_CNVs_info<-as.data.frame(matrix(unlist(strsplit(as.character(Pais_CNVs1$FORMAT),":")),ncol = 4, byrow = TRUE))
Pais_CNVs_END<-as.data.frame(matrix(unlist(strsplit(as.character(Pais_CNVs1[,8]),"=")),ncol = 2, byrow = TRUE))
Pais_CNVs<-data.frame(Pais_CNVs1[,1:2],as.numeric(as.character(Pais_CNVs_END[,2])),Pais_CNVs1[,3:7],Pais_CNVs_info)
names(Pais_CNVs)<-c("Scaffold","Start_pos","End_pos","ID","REF","ALT","QUAL","FILTER","CNV standard deviation","CNV_probability_score","CNV_copy_number","CNV_copy_number_standard_deviation")

#######################################################
#CN.mops data#
#######################################################
Krom_cnvs_HC_strict<-read.csv("Krom_cnvs_only_strict.table",header=T,sep="\t")
Lang_cnvs_HC_strict<-read.csv("Lang_cnvs_only_strict.table",header=T,sep="\t")
Kosi_cnvs_HC_strict<-read.csv("Kosi_cnvs_only_strict.table",header=T,sep="\t")
Best_cnvs_HC_strict<-read.csv("Best_cnvs_only_strict.table",header=T,sep="\t")
Pais_cnvs_HC_strict<-read.csv("Pais_cnvs_only_strict.table",header=T,sep="\t")
Noss_cnvs_HC_strict<-read.csv("Noss_cnvs_only_strict.table",header=T,sep="\t")

require(GenomicRanges)
#find overlapping regions between Krom and Kosi with similar CNs

Krom_CNVs_GRange<-GRanges(seqnames=tolower(Krom_CNVs$Scaffold),ranges=IRanges(start=Krom_CNVs$Start_pos,end=Krom_CNVs$End_pos))
values(Krom_CNVs_GRange)<-Krom_CNVs[,4:12]
Kosi_CNVs_GRange<-GRanges(seqnames=tolower(Kosi_CNVs$Scaffold),ranges=IRanges(start=Kosi_CNVs$Start_pos,end=Kosi_CNVs$End_pos))
values(Kosi_CNVs_GRange)<-Kosi_CNVs[,4:12]
KrKo_merge<-mergeByOverlaps(Krom_CNVs_GRange,Kosi_CNVs_GRange)
KrKo_merge_same<-KrKo_merge[abs(as.numeric(as.character(KrKo_merge@ listData$ Krom_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(KrKo_merge@ listData$ Kosi_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
KrKo_merge_same_df<-as.data.frame(KrKo_merge_same)[,1:3]
names(KrKo_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
KrKo_merge_same_df_GRange<-GRanges(seqnames=tolower(KrKo_merge_same_df$Scaffold),ranges=IRanges(start=KrKo_merge_same_df$Start_pos,end=KrKo_merge_same_df$End_pos))

#subset Krom CNVs to regions without similar CN in Kosi

KrKo_diff<-setdiff(Krom_CNVs_GRange,KrKo_merge_same_df_GRange)
KrKo<-subsetByOverlaps(Krom_CNVs_GRange,KrKo_diff)
KoKr_diff<-setdiff(Kosi_CNVs_GRange,KrKo_merge_same_df_GRange)
KoKr<-subsetByOverlaps(Kosi_CNVs_GRange,KoKr_diff)

#overlap cn.mops and GROM regions

Krom_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Krom_cnvs_HC_strict$Chr),ranges=IRanges(start=Krom_cnvs_HC_strict$gene_start,end=Krom_cnvs_HC_strict$gene_end))
values(Krom_cnvs_HC_strict_GRange)<-Krom_cnvs_HC_strict[,c(1:12,16:70)]
Krom_overlap_strict<-mergeByOverlaps(Krom_cnvs_HC_strict_GRange,KrKo)
Krom_overlap_strict_df<-as.data.frame(Krom_overlap_strict)
Kosi_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Kosi_cnvs_HC_strict$Chr),ranges=IRanges(start=Kosi_cnvs_HC_strict$gene_start,end=Kosi_cnvs_HC_strict$gene_end))
values(Kosi_cnvs_HC_strict_GRange)<-Kosi_cnvs_HC_strict[,c(1:12,16:70)]
Kosi_overlap_strict<-mergeByOverlaps(Kosi_cnvs_HC_strict_GRange,KoKr)
Kosi_overlap_strict_df<-as.data.frame(Kosi_overlap_strict)
Krom_CNVs_overlapping_strict<-Krom_overlap_strict_df[,c(1:5,154:162,73:139)]
names(Krom_CNVs_overlapping_strict)[1:5]<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand")
Kosi_CNVs_overlapping_strict<-Kosi_overlap_strict_df[,c(1:5,154:162,73:139)]
names(Kosi_CNVs_overlapping_strict)[1:5]<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand")
Krom_CNVs_HC_GROM_strict<-rbind(Krom_CNVs_overlapping_strict,Kosi_CNVs_overlapping_strict)

#CNVs_Krom_HC_GROM_cov<-merge(Krom_CNVs_HC_GROM_strict,Coverage_all,by.x="Ahal_Gene",by.y="ID.x",all.x=T)
CNVs_Krom_HC_GROM_cov<-Krom_CNVs_HC_GROM_strict
CNVs_Krom_HC_GROM_cov<-CNVs_Krom_HC_GROM_cov[!((CNVs_Krom_HC_GROM_cov$CN_class=="CN0"|CNVs_Krom_HC_GROM_cov$CN_class=="CN1")&CNVs_Krom_HC_GROM_cov$ALT=="<DUP>"),]
CNVs_Krom_HC_GROM_cov<-CNVs_Krom_HC_GROM_cov[!(CNVs_Krom_HC_GROM_cov$ALT=="<DEL>"&!(CNVs_Krom_HC_GROM_cov$CN_class=="CN0"|CNVs_Krom_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Krom_HC_GROM_cov,"Krom_CNVs_overlap_strict.table",row.names=F,sep="\t")
require(xlsx)
write.xlsx2(CNVs_Krom_HC_GROM_cov,"CNVs_overlap_strict.xlsx",sheetName="KromKosi",row.names=F)

Lang_CNVs_GRange<-GRanges(seqnames=tolower(Lang_CNVs$Scaffold),ranges=IRanges(start=Lang_CNVs$Start_pos,end=Lang_CNVs$End_pos))
values(Lang_CNVs_GRange)<-Lang_CNVs[,4:12]
Best_CNVs_GRange<-GRanges(seqnames=tolower(Best_CNVs$Scaffold),ranges=IRanges(start=Best_CNVs$Start_pos,end=Best_CNVs$End_pos))
values(Best_CNVs_GRange)<-Best_CNVs[,4:12]
LB_merge<-mergeByOverlaps(Lang_CNVs_GRange,Best_CNVs_GRange)
LB_merge_same<-LB_merge[abs(as.numeric(as.character(LB_merge@ listData$ Lang_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(LB_merge@ listData$ Best_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
LB_merge_same_df<-as.data.frame(LB_merge_same)[,1:3]
names(LB_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
LB_merge_same_df_GRange<-GRanges(seqnames=tolower(LB_merge_same_df$Scaffold),ranges=IRanges(start=LB_merge_same_df$Start_pos,end=LB_merge_same_df$End_pos))

#subset Lang CNVs to regions without similar CN in Best

LB_diff<-setdiff(Lang_CNVs_GRange,LB_merge_same_df_GRange)
LB<-subsetByOverlaps(Lang_CNVs_GRange,LB_diff)
BL_diff<-setdiff(Best_CNVs_GRange,LB_merge_same_df_GRange)
BL<-subsetByOverlaps(Best_CNVs_GRange,BL_diff)

#overlap cn.mops and GROM regions

Lang_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Lang_cnvs_HC_strict$Chr),ranges=IRanges(start=Lang_cnvs_HC_strict$gene_start,end=Lang_cnvs_HC_strict$gene_end))
values(Lang_cnvs_HC_strict_GRange)<-Lang_cnvs_HC_strict[,c(1:12,16:70)]
Lang_overlap_strict<-mergeByOverlaps(Lang_cnvs_HC_strict_GRange,LB)
Lang_overlap_strict_df<-as.data.frame(Lang_overlap_strict)
Best_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Best_cnvs_HC_strict$Chr),ranges=IRanges(start=Best_cnvs_HC_strict$gene_start,end=Best_cnvs_HC_strict$gene_end))
values(Best_cnvs_HC_strict_GRange)<-Best_cnvs_HC_strict[,c(1:12,16:70)]
Best_overlap_strict<-mergeByOverlaps(Best_cnvs_HC_strict_GRange,BL)
Best_overlap_strict_df<-as.data.frame(Best_overlap_strict)
Lang_CNVs_overlapping_strict<-Lang_overlap_strict_df[,c(1:5,154:162,73:139)]
names(Lang_CNVs_overlapping_strict)[1:5]<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand")
Best_CNVs_overlapping_strict<-Best_overlap_strict_df[,c(1:5,154:162,73:139)]
names(Best_CNVs_overlapping_strict)[1:5]<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand")
Lang_CNVs_HC_GROM_strict<-rbind(Lang_CNVs_overlapping_strict,Best_CNVs_overlapping_strict)

#CNVs_Lang_HC_GROM_cov<-merge(Lang_CNVs_HC_GROM_strict,Coverage_all,by.x="Ahal_Gene",by.y="ID.x",all.x=T)
CNVs_Lang_HC_GROM_cov<-Lang_CNVs_HC_GROM_strict
CNVs_Lang_HC_GROM_cov<-CNVs_Lang_HC_GROM_cov[!((CNVs_Lang_HC_GROM_cov$CN_class=="CN0"|CNVs_Lang_HC_GROM_cov$CN_class=="CN1")&CNVs_Lang_HC_GROM_cov$ALT=="<DUP>"),]
CNVs_Lang_HC_GROM_cov<-CNVs_Lang_HC_GROM_cov[!(CNVs_Lang_HC_GROM_cov$ALT=="<DEL>"&!(CNVs_Lang_HC_GROM_cov$CN_class=="CN0"|CNVs_Lang_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Lang_HC_GROM_cov,"Lang_CNVs_overlap_strict.table",row.names=F,sep="\t")
require(xlsx)
write.xlsx2(CNVs_Lang_HC_GROM_cov,"CNVs_overlap_strict.xlsx",sheetName="LangBest",row.names=F,append=T)


Noss_CNVs_GRange<-GRanges(seqnames=tolower(Noss_CNVs$Scaffold),ranges=IRanges(start=Noss_CNVs$Start_pos,end=Noss_CNVs$End_pos))
values(Noss_CNVs_GRange)<-Noss_CNVs[,4:12]
Pais_CNVs_GRange<-GRanges(seqnames=tolower(Pais_CNVs$Scaffold),ranges=IRanges(start=Pais_CNVs$Start_pos,end=Pais_CNVs$End_pos))
values(Pais_CNVs_GRange)<-Pais_CNVs[,4:12]
NP_merge<-mergeByOverlaps(Noss_CNVs_GRange,Pais_CNVs_GRange)
NP_merge_same<-NP_merge[abs(as.numeric(as.character(NP_merge@ listData$ Noss_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number))-as.numeric(as.character(NP_merge@ listData$ Pais_CNVs_GRange@ elementMetadata@ listData$ CNV_copy_number)))<0.8,]
NP_merge_same_df<-as.data.frame(NP_merge_same)[,1:3]
names(NP_merge_same_df)<-c("Scaffold","Start_pos","End_pos")
NP_merge_same_df_GRange<-GRanges(seqnames=tolower(NP_merge_same_df$Scaffold),ranges=IRanges(start=NP_merge_same_df$Start_pos,end=NP_merge_same_df$End_pos))

#subset Noss CNVs to regions without similar CN in Pais

NP_diff<-setdiff(Noss_CNVs_GRange,NP_merge_same_df_GRange)
NP<-subsetByOverlaps(Noss_CNVs_GRange,NP_diff)
PN_diff<-setdiff(Pais_CNVs_GRange,NP_merge_same_df_GRange)
PN<-subsetByOverlaps(Pais_CNVs_GRange,PN_diff)

#overlap cn.mops and GROM regions

Noss_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Noss_cnvs_HC_strict$Chr),ranges=IRanges(start=Noss_cnvs_HC_strict$gene_start,end=Noss_cnvs_HC_strict$gene_end))
values(Noss_cnvs_HC_strict_GRange)<-Noss_cnvs_HC_strict[,c(1:12,16:70)]
Noss_overlap_strict<-mergeByOverlaps(Noss_cnvs_HC_strict_GRange,NP)
Noss_overlap_strict_df<-as.data.frame(Noss_overlap_strict)
Pais_cnvs_HC_strict_GRange<-GRanges(seqnames=tolower(Pais_cnvs_HC_strict$Chr),ranges=IRanges(start=Pais_cnvs_HC_strict$gene_start,end=Pais_cnvs_HC_strict$gene_end))
values(Pais_cnvs_HC_strict_GRange)<-Pais_cnvs_HC_strict[,c(1:12,16:70)]
Pais_overlap_strict<-mergeByOverlaps(Pais_cnvs_HC_strict_GRange,PN)
Pais_overlap_strict_df<-as.data.frame(Pais_overlap_strict)
Noss_CNVs_overlapping_strict<-Noss_overlap_strict_df[,c(1:5,154:162,73:139)]
names(Noss_CNVs_overlapping_strict)[1:5]<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand")
Pais_CNVs_overlapping_strict<-Pais_overlap_strict_df[,c(1:5,154:162,73:139)]
names(Pais_CNVs_overlapping_strict)[1:5]<-c("Scaffold","Overlap_start","Overlap_end","Overlap_width","Overlap_strand")
Noss_CNVs_HC_GROM_strict<-rbind(Noss_CNVs_overlapping_strict,Pais_CNVs_overlapping_strict)

#CNVs_Noss_HC_GROM_cov<-merge(Noss_CNVs_HC_GROM_strict,Coverage_all,by.x="Ahal_Gene",by.y="ID.x",all.x=T)
CNVs_Noss_HC_GROM_cov<-Noss_CNVs_HC_GROM_strict
CNVs_Noss_HC_GROM_cov<-CNVs_Noss_HC_GROM_cov[!((CNVs_Noss_HC_GROM_cov$CN_class=="CN0"|CNVs_Noss_HC_GROM_cov$CN_class=="CN1")&CNVs_Noss_HC_GROM_cov$ALT=="<DUP>"),]
CNVs_Noss_HC_GROM_cov<-CNVs_Noss_HC_GROM_cov[!(CNVs_Noss_HC_GROM_cov$ALT=="<DEL>"&!(CNVs_Noss_HC_GROM_cov$CN_class=="CN0"|CNVs_Noss_HC_GROM_cov$CN_class=="CN1")),]
write.table(CNVs_Noss_HC_GROM_cov,"Noss_CNVs_overlap_strict.table",row.names=F,sep="\t")
require(xlsx)
write.xlsx2(CNVs_Noss_HC_GROM_cov,"CNVs_overlap_strict.xlsx",sheetName="NossPais",row.names=F,append=T)


###################################################
#INDEL analysis#
###################################################
require(xlsx)

############################################################
#GROM data#
############################################################

Lang_indel_grom1<-read.table("Lang.GROM.Indel.ann.vcf",sep="\t",header=F,fill=T,quote="")
Best_indel_grom1<-read.table("Best.GROM.Indel.ann.vcf",sep="\t",header=F,fill=T,quote="")
Noss_indel_grom1<-read.table("Noss.GROM.Indel.ann.vcf",sep="\t",header=F,fill=T,quote="")
Pais_indel_grom1<-read.table("Pais.GROM.Indel.ann.vcf",sep="\t",header=F,fill=T,quote="")
Krom_indel_grom1<-read.table("Krom.GROM.Indel.ann.vcf",sep="\t",header=F,fill=T,quote="")
Kosi_indel_grom1<-read.table("Kosi.GROM.Indel.ann.vcf",sep="\t",header=F,fill=T,quote="")
Krom_indel_grom<-Krom_indel_grom1[paste(as.character(Krom_indel_grom1$V1),as.character(Krom_indel_grom1$V2),sep=":")%in%paste(as.character(Krom_Indels$Scaffold),as.character(Krom_Indels$Pos),sep=":"),]
Kosi_indel_grom<-Kosi_indel_grom1[paste(as.character(Kosi_indel_grom1$V1),as.character(Kosi_indel_grom1$V2),sep=":")%in%paste(as.character(Kosi_Indels$Scaffold),as.character(Kosi_Indels$Pos),sep=":"),]
Lang_indel_grom<-Lang_indel_grom1[paste(as.character(Lang_indel_grom1$V1),as.character(Lang_indel_grom1$V2),sep=":")%in%paste(as.character(Lang_Indels$Scaffold),as.character(Lang_Indels$Pos),sep=":"),]
Best_indel_grom<-Best_indel_grom1[paste(as.character(Best_indel_grom1$V1),as.character(Best_indel_grom1$V2),sep=":")%in%paste(as.character(Best_Indels$Scaffold),as.character(Best_Indels$Pos),sep=":"),]
Noss_indel_grom<-Noss_indel_grom1[paste(as.character(Noss_indel_grom1$V1),as.character(Noss_indel_grom1$V2),sep=":")%in%paste(as.character(Noss_Indels$Scaffold),as.character(Noss_Indels$Pos),sep=":"),]
Pais_indel_grom<-Pais_indel_grom1[paste(as.character(Pais_indel_grom1$V1),as.character(Pais_indel_grom1$V2),sep=":")%in%paste(as.character(Pais_Indels$Scaffold),as.character(Pais_Indels$Pos),sep=":"),]

names(Lang_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")
names(Best_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")
names(Noss_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")
names(Pais_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")
names(Krom_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")
names(Kosi_indel_grom)<-c("Scaffold","Pos","ID","REF","ALT","QUAL","FILTER","ANN","FORMAT_names","FORMAT")

Krom_indel_gromAnnsplit<-strsplit(as.character(Krom_indel_grom$ANN),",",fixed=T)
Krom_indel_gromAnnsplit1<-lapply(Krom_indel_gromAnnsplit,"[",1:max(unlist(lapply(Krom_indel_gromAnnsplit,length))))
Krom_indel_gromAnnsplit2<-data.frame(matrix(unlist(Krom_indel_gromAnnsplit1),nrow=length(Krom_indel_gromAnnsplit1),byrow=T))
Krom_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Krom_indel_gromAnnsplit,length))))
	{Krom_indel_grom1<-ifelse(!is.na(Krom_indel_gromAnnsplit2[,i]),strsplit(as.character(Krom_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Krom_indel_gromAnn1<-lapply(Krom_indel_grom1,"[",c(1:3,5))
	Krom_indel_gromAnndf1<-data.frame(matrix(unlist(Krom_indel_gromAnn1),nrow=length(Krom_indel_gromAnnsplit1),byrow=T))
	names(Krom_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Krom_indel_gromcomp1<-cbind(Krom_indel_grom[,1:7],Krom_indel_gromAnndf1)
	Krom_indel_gromdf<-rbind(Krom_indel_gromdf,na.exclude(Krom_indel_gromcomp1))
	}
Krom_indel_gromdf2<-Krom_indel_gromdf
Krom_indel_gromdf<-Krom_indel_gromdf2
Gene<-gsub("^.*-","",Krom_indel_gromdf2$Gene)
Krom_indel_gromdf$Gene<-Gene
Kosi_indel_gromAnnsplit<-strsplit(as.character(Kosi_indel_grom$ANN),",",fixed=T)
Kosi_indel_gromAnnsplit1<-lapply(Kosi_indel_gromAnnsplit,"[",1:max(unlist(lapply(Kosi_indel_gromAnnsplit,length))))
Kosi_indel_gromAnnsplit2<-data.frame(matrix(unlist(Kosi_indel_gromAnnsplit1),nrow=length(Kosi_indel_gromAnnsplit1),byrow=T))
Kosi_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Kosi_indel_gromAnnsplit,length))))
	{Kosi_indel_grom1<-ifelse(!is.na(Kosi_indel_gromAnnsplit2[,i]),strsplit(as.character(Kosi_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Kosi_indel_gromAnn1<-lapply(Kosi_indel_grom1,"[",c(1:3,5))
	Kosi_indel_gromAnndf1<-data.frame(matrix(unlist(Kosi_indel_gromAnn1),nrow=length(Kosi_indel_gromAnnsplit1),byrow=T))
	names(Kosi_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Kosi_indel_gromcomp1<-cbind(Kosi_indel_grom[,1:7],Kosi_indel_gromAnndf1)
	Kosi_indel_gromdf<-rbind(Kosi_indel_gromdf,na.exclude(Kosi_indel_gromcomp1))
	}
require(GenomicRanges)
Kosi_indel_gromdf2<-Kosi_indel_gromdf
Kosi_indel_gromdf<-Kosi_indel_gromdf2
Gene<-gsub("^.*-","",Kosi_indel_gromdf2$Gene)
Kosi_indel_gromdf$Gene<-Gene
KK_INDEL<-read.table("SNP_KK_Indels_diffabove04_hQD.txt",sep="\t",header=T,fill=T,quote="")

#subset to main scaffolds because otherwise NAs

require(GenomicRanges)
KromKosi_Indels2<-rbind(Krom_indel_gromdf,Kosi_indel_gromdf)
KromKosi_Indels2<-KromKosi_Indels2[KromKosi_Indels2$Scaffold %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),]
KK_INDEL<-KK_INDEL[KK_INDEL$Chrom %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),]
KK_INDEL_GR<-GRanges(seqnames=tolower(KK_INDEL$Chrom),ranges=IRanges(start=as.numeric(KK_INDEL$Pos),end=as.numeric(KK_INDEL$Pos)))
KromKosi_Indels_GR<-GRanges(seqnames=tolower(KromKosi_Indels2$Scaffold),ranges=IRanges(start=KromKosi_Indels2$Pos,end=KromKosi_Indels2$Pos))
values(KK_INDEL_GR)<-KK_INDEL[,c(1,4:64)]
values(KromKosi_Indels_GR)<-KromKosi_Indels2[,3:11]
KK_Indels_nearest<-nearest(KK_INDEL_GR,KromKosi_Indels_GR,ignore.strand=T)
KK_nearest<-KromKosi_Indels2[KK_Indels_nearest,1:11]
KK_Indel_table<-cbind(KK_INDEL,KK_nearest)
KK_Indel_table2<-KK_Indel_table[abs(KK_Indel_table[,3]-KK_Indel_table[,66])<=max(nchar(as.character(KK_Indel_table[,8])),nchar(as.character(KK_Indel_table[,9])),nchar(as.character(KK_Indel_table[,68])),nchar(as.character(KK_Indel_table[,69]))),]
SNPeffdf<-KK_Indel_table2

SNPeffdf2<-SNPeffdf[SNPeffdf[,12]=="HIGH"&SNPeffdf[,74]=="HIGH",]
SNPeffdf3<-SNPeffdf[SNPeffdf[,12]=="HIGH"|SNPeffdf[,74]=="HIGH",]
nrow(SNPeffdf2[!duplicated(SNPeffdf2),])
#9
nrow(SNPeffdf3[!duplicated(SNPeffdf3),])
#114
SNPeffdf2<-SNPeffdf2[!duplicated(SNPeffdf2),]
SNPeffdf3<-SNPeffdf3[!duplicated(SNPeffdf3),]

SNPeffdf2<-SNPeffdf2[,c(1:12,65:75,13:64)]
SNPeffdf3<-SNPeffdf3[,c(1:12,65:75,13:64)]

#remove duplicates where different high effect mutations occurred in both populations
dup_genes<-SNPeffdf2$Ah_ID[duplicated(SNPeffdf2$Ah_ID)]
unique_pos<-!(duplicated(SNPeffdf2$Pos))
Gene_tab<-tabulate(SNPeffdf2$Ah_ID)
sum_AF<-sum(SNPeffdf2$AF_Krom>SNPeffdf2$AF_Kosi)

removal_df<-SNPeffdf2[SNPeffdf2$Ah_ID%in%dup_genes&unique_pos=="TRUE",]
removal<-c()
for(i in removal_df$Ah_ID)
	{for (j in unique(removal_df$Ath_ID[removal_df$Ah_ID==i]))
		{if (is.na(j)==TRUE)
			{removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i,])!=sum(SNPeffdf2$AF_Krom[SNPeffdf2$Ah_ID==i]>SNPeffdf2$AF_Kosi[SNPeffdf2$Ah_ID==i])&nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i,])!=sum(SNPeffdf2$AF_Krom[SNPeffdf2$Ah_ID==i]<SNPeffdf2$AF_Kosi[SNPeffdf2$Ah_ID==i]),paste0(i),"keep"))
			}
		else
			{
			removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j,])!=sum(SNPeffdf2$AF_Krom[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]>SNPeffdf2$AF_Kosi[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j])&nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j,])!=sum(SNPeffdf2$AF_Krom[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]<SNPeffdf2$AF_Kosi[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]),paste0(i,"_",j),"keep"))
			}
		}
	}


SNPeffdf2_clean<-SNPeffdf2[!(paste0(SNPeffdf2$Ah_ID,"_",SNPeffdf2$Ath_ID)%in%removal|(is.na(SNPeffdf2$Ath_ID)&SNPeffdf2$Ah_ID%in%removal)),]

write.table(SNPeffdf2_clean,"Indel_KK_HC_GROM.table",row.names=F,sep="\t")
write.xlsx(SNPeffdf2_clean,"Indel_KK_HC_GROM.xlsx",sheetName="Indel_KromKosi_HC_GROM",row.names=F)


Lang_indel_gromAnnsplit<-strsplit(as.character(Lang_indel_grom$ANN),",",fixed=T)
Lang_indel_gromAnnsplit1<-lapply(Lang_indel_gromAnnsplit,"[",1:max(unlist(lapply(Lang_indel_gromAnnsplit,length))))
Lang_indel_gromAnnsplit2<-data.frame(matrix(unlist(Lang_indel_gromAnnsplit1),nrow=length(Lang_indel_gromAnnsplit1),byrow=T))
Lang_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Lang_indel_gromAnnsplit,length))))
	{Lang_indel_grom1<-ifelse(!is.na(Lang_indel_gromAnnsplit2[,i]),strsplit(as.character(Lang_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Lang_indel_gromAnn1<-lapply(Lang_indel_grom1,"[",c(1:3,5))
	Lang_indel_gromAnndf1<-data.frame(matrix(unlist(Lang_indel_gromAnn1),nrow=length(Lang_indel_gromAnnsplit1),byrow=T))
	names(Lang_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Lang_indel_gromcomp1<-cbind(Lang_indel_grom[,1:7],Lang_indel_gromAnndf1)
	Lang_indel_gromdf<-rbind(Lang_indel_gromdf,na.exclude(Lang_indel_gromcomp1))
	}
Lang_indel_gromdf2<-Lang_indel_gromdf
Lang_indel_gromdf<-Lang_indel_gromdf2
Gene<-gsub("^.*-","",Lang_indel_gromdf2$Gene)
Lang_indel_gromdf$Gene<-Gene
Best_indel_gromAnnsplit<-strsplit(as.character(Best_indel_grom$ANN),",",fixed=T)
Best_indel_gromAnnsplit1<-lapply(Best_indel_gromAnnsplit,"[",1:max(unlist(lapply(Best_indel_gromAnnsplit,length))))
Best_indel_gromAnnsplit2<-data.frame(matrix(unlist(Best_indel_gromAnnsplit1),nrow=length(Best_indel_gromAnnsplit1),byrow=T))
Best_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Best_indel_gromAnnsplit,length))))
	{Best_indel_grom1<-ifelse(!is.na(Best_indel_gromAnnsplit2[,i]),strsplit(as.character(Best_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Best_indel_gromAnn1<-lapply(Best_indel_grom1,"[",c(1:3,5))
	Best_indel_gromAnndf1<-data.frame(matrix(unlist(Best_indel_gromAnn1),nrow=length(Best_indel_gromAnnsplit1),byrow=T))
	names(Best_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Best_indel_gromcomp1<-cbind(Best_indel_grom[,1:7],Best_indel_gromAnndf1)
	Best_indel_gromdf<-rbind(Best_indel_gromdf,na.exclude(Best_indel_gromcomp1))
	}
require(GenomicRanges)
Best_indel_gromdf2<-Best_indel_gromdf
Best_indel_gromdf<-Best_indel_gromdf2
Gene<-gsub("^.*-","",Best_indel_gromdf2$Gene)
Best_indel_gromdf$Gene<-Gene
LB_INDEL<-read.table("SNP_LB_Indels_diffabove04_hQD.txt",sep="\t",header=T,fill=T,quote="")

#subset to main scaffolds because otherwise NAs

require(GenomicRanges)
LangBest_Indels2<-rbind(Lang_indel_gromdf,Best_indel_gromdf)
LangBest_Indels2<-LangBest_Indels2[LangBest_Indels2$Scaffold %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),]
LB_INDEL<-LB_INDEL[LB_INDEL$Chrom %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),]
LB_INDEL_GR<-GRanges(seqnames=tolower(LB_INDEL$Chrom),ranges=IRanges(start=as.numeric(LB_INDEL$Pos),end=as.numeric(LB_INDEL$Pos)))
LangBest_Indels_GR<-GRanges(seqnames=tolower(LangBest_Indels2$Scaffold),ranges=IRanges(start=LangBest_Indels2$Pos,end=LangBest_Indels2$Pos))
values(LB_INDEL_GR)<-LB_INDEL[,c(1,4:64)]
values(LangBest_Indels_GR)<-LangBest_Indels2[,3:11]
LB_Indels_nearest<-nearest(LB_INDEL_GR,LangBest_Indels_GR,ignore.strand=T)
LB_nearest<-LangBest_Indels2[LB_Indels_nearest,1:11]
LB_Indel_table<-cbind(LB_INDEL,LB_nearest)
LB_Indel_table2<-LB_Indel_table[abs(LB_Indel_table[,3]-LB_Indel_table[,66])<=max(nchar(as.character(LB_Indel_table[,8])),nchar(as.character(LB_Indel_table[,9])),nchar(as.character(LB_Indel_table[,68])),nchar(as.character(LB_Indel_table[,69]))),]
SNPeffdf<-LB_Indel_table2

SNPeffdf2<-SNPeffdf[SNPeffdf[,12]=="HIGH"&SNPeffdf[,74]=="HIGH",]
SNPeffdf3<-SNPeffdf[SNPeffdf[,12]=="HIGH"|SNPeffdf[,74]=="HIGH",]
nrow(SNPeffdf2[!duplicated(SNPeffdf2),])
#119
nrow(SNPeffdf3[!duplicated(SNPeffdf3),])
#1255
SNPeffdf2<-SNPeffdf2[!duplicated(SNPeffdf2),]
SNPeffdf3<-SNPeffdf3[!duplicated(SNPeffdf3),]

SNPeffdf2<-SNPeffdf2[,c(1:12,65:75,13:64)]
SNPeffdf3<-SNPeffdf3[,c(1:12,65:75,13:64)]

#remove duplicates where different high effect mutations occurred in both populations
dup_genes<-SNPeffdf2$Ah_ID[duplicated(SNPeffdf2$Ah_ID)]
unique_pos<-!(duplicated(SNPeffdf2$Pos))
Gene_tab<-tabulate(SNPeffdf2$Ah_ID)
sum_AF<-sum(SNPeffdf2$AF_Lang>SNPeffdf2$AF_Best)

removal_df<-SNPeffdf2[SNPeffdf2$Ah_ID%in%dup_genes&unique_pos=="TRUE",]
removal<-c()
for(i in removal_df$Ah_ID)
	{for (j in unique(removal_df$Ath_ID[removal_df$Ah_ID==i]))
		{if (is.na(j)==TRUE)
			{removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i,])!=sum(SNPeffdf2$AF_Lang[SNPeffdf2$Ah_ID==i]>SNPeffdf2$AF_Best[SNPeffdf2$Ah_ID==i])&nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i,])!=sum(SNPeffdf2$AF_Lang[SNPeffdf2$Ah_ID==i]<SNPeffdf2$AF_Best[SNPeffdf2$Ah_ID==i]),paste0(i),"keep"))
			}
		else
			{
			removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j,])!=sum(SNPeffdf2$AF_Lang[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]>SNPeffdf2$AF_Best[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j])&nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j,])!=sum(SNPeffdf2$AF_Lang[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]<SNPeffdf2$AF_Best[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]),paste0(i,"_",j),"keep"))
			}
		}
	}


SNPeffdf2_clean<-SNPeffdf2[!(paste0(SNPeffdf2$Ah_ID,"_",SNPeffdf2$Ath_ID)%in%removal|(is.na(SNPeffdf2$Ath_ID)&SNPeffdf2$Ah_ID%in%removal)),]

write.table(SNPeffdf2_clean,"Indel_LB_HC_GROM.table",row.names=F,sep="\t")
write.xlsx(SNPeffdf2_clean,"Indel_LB_HC_GROM.xlsx",sheetName="Indel_LangBest_HC_GROM",row.names=F)


Noss_indel_gromAnnsplit<-strsplit(as.character(Noss_indel_grom$ANN),",",fixed=T)
Noss_indel_gromAnnsplit1<-lapply(Noss_indel_gromAnnsplit,"[",1:max(unlist(lapply(Noss_indel_gromAnnsplit,length))))
Noss_indel_gromAnnsplit2<-data.frame(matrix(unlist(Noss_indel_gromAnnsplit1),nrow=length(Noss_indel_gromAnnsplit1),byrow=T))
Noss_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Noss_indel_gromAnnsplit,length))))
	{Noss_indel_grom1<-ifelse(!is.na(Noss_indel_gromAnnsplit2[,i]),strsplit(as.character(Noss_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Noss_indel_gromAnn1<-lapply(Noss_indel_grom1,"[",c(1:3,5))
	Noss_indel_gromAnndf1<-data.frame(matrix(unlist(Noss_indel_gromAnn1),nrow=length(Noss_indel_gromAnnsplit1),byrow=T))
	names(Noss_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Noss_indel_gromcomp1<-cbind(Noss_indel_grom[,1:7],Noss_indel_gromAnndf1)
	Noss_indel_gromdf<-rbind(Noss_indel_gromdf,na.exclude(Noss_indel_gromcomp1))
	}
Noss_indel_gromdf2<-Noss_indel_gromdf
Noss_indel_gromdf<-Noss_indel_gromdf2
Gene<-gsub("^.*-","",Noss_indel_gromdf2$Gene)
Noss_indel_gromdf$Gene<-Gene
Pais_indel_gromAnnsplit<-strsplit(as.character(Pais_indel_grom$ANN),",",fixed=T)
Pais_indel_gromAnnsplit1<-lapply(Pais_indel_gromAnnsplit,"[",1:max(unlist(lapply(Pais_indel_gromAnnsplit,length))))
Pais_indel_gromAnnsplit2<-data.frame(matrix(unlist(Pais_indel_gromAnnsplit1),nrow=length(Pais_indel_gromAnnsplit1),byrow=T))
Pais_indel_gromdf<-data.frame()
for (i in 1:max(unlist(lapply(Pais_indel_gromAnnsplit,length))))
	{Pais_indel_grom1<-ifelse(!is.na(Pais_indel_gromAnnsplit2[,i]),strsplit(as.character(Pais_indel_gromAnnsplit2[,i]),"|",fixed=T),NA)
	Pais_indel_gromAnn1<-lapply(Pais_indel_grom1,"[",c(1:3,5))
	Pais_indel_gromAnndf1<-data.frame(matrix(unlist(Pais_indel_gromAnn1),nrow=length(Pais_indel_gromAnnsplit1),byrow=T))
	names(Pais_indel_gromAnndf1)<-c("Base","Characterization","Effect","Gene")
	Pais_indel_gromcomp1<-cbind(Pais_indel_grom[,1:7],Pais_indel_gromAnndf1)
	Pais_indel_gromdf<-rbind(Pais_indel_gromdf,na.exclude(Pais_indel_gromcomp1))
	}
require(GenomicRanges)
Pais_indel_gromdf2<-Pais_indel_gromdf
Pais_indel_gromdf<-Pais_indel_gromdf2
Gene<-gsub("^.*-","",Pais_indel_gromdf2$Gene)
Pais_indel_gromdf$Gene<-Gene
NP_INDEL<-read.table("SNP_NP_Indels_diffabove04_hQD.txt",sep="\t",header=T,fill=T,quote="")

#subset to main scaffolds because otherwise NAs

require(GenomicRanges)
NossPais_Indels2<-rbind(Noss_indel_gromdf,Pais_indel_gromdf)
NossPais_Indels2<-NossPais_Indels2[NossPais_Indels2$Scaffold %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),]
NP_INDEL<-NP_INDEL[NP_INDEL$Chrom %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"),]
NP_INDEL_GR<-GRanges(seqnames=tolower(NP_INDEL$Chrom),ranges=IRanges(start=as.numeric(NP_INDEL$Pos),end=as.numeric(NP_INDEL$Pos)))
NossPais_Indels_GR<-GRanges(seqnames=tolower(NossPais_Indels2$Scaffold),ranges=IRanges(start=NossPais_Indels2$Pos,end=NossPais_Indels2$Pos))
values(NP_INDEL_GR)<-NP_INDEL[,c(1,4:64)]
values(NossPais_Indels_GR)<-NossPais_Indels2[,3:11]
NP_Indels_nearest<-nearest(NP_INDEL_GR,NossPais_Indels_GR,ignore.strand=T)
NP_nearest<-NossPais_Indels2[NP_Indels_nearest,1:11]
NP_Indel_table<-cbind(NP_INDEL,NP_nearest)
NP_Indel_table2<-NP_Indel_table[abs(NP_Indel_table[,3]-NP_Indel_table[,66])<=max(nchar(as.character(NP_Indel_table[,8])),nchar(as.character(NP_Indel_table[,9])),nchar(as.character(NP_Indel_table[,68])),nchar(as.character(NP_Indel_table[,69]))),]
SNPeffdf<-NP_Indel_table2

SNPeffdf2<-SNPeffdf[SNPeffdf[,12]=="HIGH"&SNPeffdf[,74]=="HIGH",]
SNPeffdf3<-SNPeffdf[SNPeffdf[,12]=="HIGH"|SNPeffdf[,74]=="HIGH",]
nrow(SNPeffdf2[!duplicated(SNPeffdf2),])
#50
nrow(SNPeffdf3[!duplicated(SNPeffdf3),])
#758
SNPeffdf2<-SNPeffdf2[!duplicated(SNPeffdf2),]
SNPeffdf3<-SNPeffdf3[!duplicated(SNPeffdf3),]

SNPeffdf2<-SNPeffdf2[,c(1:12,65:75,13:64)]
SNPeffdf3<-SNPeffdf3[,c(1:12,65:75,13:64)]

#remove duplicates where different high effect mutations occurred in both populations
dup_genes<-SNPeffdf2$Ah_ID[duplicated(SNPeffdf2$Ah_ID)]
unique_pos<-!(duplicated(SNPeffdf2$Pos))
Gene_tab<-tabulate(SNPeffdf2$Ah_ID)
sum_AF<-sum(SNPeffdf2$AF_Noss>SNPeffdf2$AF_Pais)

removal_df<-SNPeffdf2[SNPeffdf2$Ah_ID%in%dup_genes&unique_pos=="TRUE",]
removal<-c()
for(i in removal_df$Ah_ID)
	{for (j in unique(removal_df$Ath_ID[removal_df$Ah_ID==i]))
		{if (is.na(j)==TRUE)
			{removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i,])!=sum(SNPeffdf2$AF_Noss[SNPeffdf2$Ah_ID==i]>SNPeffdf2$AF_Pais[SNPeffdf2$Ah_ID==i])&nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i,])!=sum(SNPeffdf2$AF_Noss[SNPeffdf2$Ah_ID==i]<SNPeffdf2$AF_Pais[SNPeffdf2$Ah_ID==i]),paste0(i),"keep"))
			}
		else
			{
			removal<-c(removal,ifelse(nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j,])!=sum(SNPeffdf2$AF_Noss[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]>SNPeffdf2$AF_Pais[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j])&nrow(SNPeffdf2[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j,])!=sum(SNPeffdf2$AF_Noss[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]<SNPeffdf2$AF_Pais[SNPeffdf2$Ah_ID==i&SNPeffdf2$Ath_ID==j]),paste0(i,"_",j),"keep"))
			}
		}
	}


SNPeffdf2_clean<-SNPeffdf2[!(paste0(SNPeffdf2$Ah_ID,"_",SNPeffdf2$Ath_ID)%in%removal|(is.na(SNPeffdf2$Ath_ID)&SNPeffdf2$Ah_ID%in%removal)),]

write.table(SNPeffdf2_clean,"Indel_NP_HC_GROM.table",row.names=F,sep="\t")
write.xlsx(SNPeffdf2_clean,"Indel_NP_HC_GROM.xlsx",sheetName="Indel_NossPais_HC_GROM",row.names=F)


NP<-read.table("Indel_NP_HC_GROM.table",sep="\t",header=T)
KK<-read.table("Indel_KK_HC_GROM.table",sep="\t",header=T)
LB<-read.table("Indel_LB_HC_GROM.table",sep="\t",header=T)

write.xlsx2(KK,"Indel_HC_GROM.xlsx",sheetName="KromKosi",row.names=F)
write.xlsx2(NP,"Indel_HC_GROM.xlsx",sheetName="NossPais",row.names=F,append=T)
write.xlsx2(LB,"Indel_HC_GROM.xlsx",sheetName="LangBest",row.names=F,append=T)


