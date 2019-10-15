options(java.parameters = "-Xmx12000m")
require(xlsx)

FST<-read.table("Genes01percentFst_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")
DD<-read.table("Genes01percentDD_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")
NIELSEN<-read.table("Genes01percentNielsen_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")
DXY<-read.table("Genes01percentDxy_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")
FLK<-read.table("Genes01percentFlk_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")
VARLD<-read.table("Genes01percentVarLD_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")
AFDabs<-read.table("Genes01percentAFDabs_hallerigenome_TDmean_KromKosi.txt",sep="\t",header=T,fill=T,quote="")

TajimasD_Krom<-read.csv("TD_Krom.csv",header=T,sep="\t")
names(TajimasD_Krom)<-c("Scaffold","Start_pos","End_pos","TajimasDKrom")

TajimasD_Kosi<-read.csv("TD_Kosi.csv",header=T,sep="\t")
names(TajimasD_Kosi)<-c("Scaffold","Start_pos","End_pos","TajimasDKosi")

FST$Window_start<-as.numeric(as.character(FST$Window_start))
FST$Window_end<-as.numeric(as.character(FST$Window_end))
DD$Window_start<-as.numeric(as.character(DD$Window_start))
DD$Window_end<-as.numeric(as.character(DD$Window_end))
DXY$Window_start<-as.numeric(as.character(DXY$Window_start))
DXY$Window_end<-as.numeric(as.character(DXY$Window_end))
AFDabs$Window_start<-as.numeric(as.character(AFDabs$Window_start))
AFDabs$Window_end<-as.numeric(as.character(AFDabs$Window_end))
NIELSEN$Window_start<-as.numeric(as.character(NIELSEN$Window_start))
NIELSEN$Window_end<-as.numeric(as.character(NIELSEN$Window_end))
VARLD$Window_start<-as.numeric(as.character(VARLD$Window_start))
VARLD$Window_end<-as.numeric(as.character(VARLD$Window_end))
FLK$Window_start<-as.numeric(as.character(FLK$Window_start))
FLK$Window_end<-as.numeric(as.character(FLK$Window_end))

require(GenomicRanges)
FST_GRange<-GRanges(seqnames=tolower(FST$Chr),ranges=IRanges(start=FST$Window_start,end=FST$Window_end))
DXY_GRange<-GRanges(seqnames=tolower(DXY$Chr),ranges=IRanges(start=DXY$Window_start,end=DXY$Window_end))
AFDabs_GRange<-GRanges(seqnames=tolower(AFDabs$Chr),ranges=IRanges(start=AFDabs$Window_start,end=AFDabs$Window_end))
NIELSEN_GRange<-GRanges(seqnames=tolower(NIELSEN$Chr),ranges=IRanges(start=NIELSEN$Window_start,end=NIELSEN$Window_end))
VARLD_GRange<-GRanges(seqnames=tolower(VARLD$Chr),ranges=IRanges(start=VARLD$Window_start,end=VARLD$Window_end))
FLK_GRange<-GRanges(seqnames=tolower(FLK$Chr),ranges=IRanges(start=FLK$Window_start,end=FLK$Window_end))
DD_GRange<-GRanges(seqnames=tolower(DD$Chr),ranges=IRanges(start=DD$Window_start,end=DD$Window_end))

values(FST_GRange)<-FST[,c(1,12)]
values(DD_GRange)<-DD[,c(1,12)]
values(DXY_GRange)<-DXY[,c(1,12)]
values(AFDabs_GRange)<-AFDabs[,c(1,12)]
values(NIELSEN_GRange)<-NIELSEN[,c(1,12)]
values(VARLD_GRange)<-VARLD[,c(1,12)]
values(FLK_GRange)<-FLK[,c(1,12)]

TajimasD_Krom_GRange<-GRanges(seqnames=tolower(TajimasD_Krom$Scaffold),ranges=IRanges(start=TajimasD_Krom$Start_pos,end=TajimasD_Krom$End_pos))
TajimasD_Kosi_GRange<-GRanges(seqnames=tolower(TajimasD_Kosi$Scaffold),ranges=IRanges(start=TajimasD_Kosi$Start_pos,end=TajimasD_Kosi$End_pos))

values(TajimasD_Krom_GRange)<-TajimasD_Krom[,4]
values(TajimasD_Kosi_GRange)<-TajimasD_Kosi[,4]

TD_Krom_nearestFST<-nearest(FST_GRange,TajimasD_Krom_GRange,ignore.strand=T)
FST$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestFST,4]
TD_Krom_nearestDXY<-nearest(DXY_GRange,TajimasD_Krom_GRange,ignore.strand=T)
DXY$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestDXY,4]
TD_Krom_nearestAFDabs<-nearest(AFDabs_GRange,TajimasD_Krom_GRange,ignore.strand=T)
AFDabs$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestAFDabs,4]
TD_Krom_nearestNIELSEN<-nearest(NIELSEN_GRange,TajimasD_Krom_GRange,ignore.strand=T)
NIELSEN$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestNIELSEN,4]
TD_Krom_nearestVARLD<-nearest(VARLD_GRange,TajimasD_Krom_GRange,ignore.strand=T)
VARLD$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestVARLD,4]
TD_Krom_nearestFLK<-nearest(FLK_GRange,TajimasD_Krom_GRange,ignore.strand=T)
FLK$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestFLK,4]
TD_Krom_nearestDD<-nearest(DD_GRange,TajimasD_Krom_GRange,ignore.strand=T)
DD$TD_Krom_nearest<-TajimasD_Krom[TD_Krom_nearestDD,4]

TD_Kosi_nearestFST<-nearest(FST_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
FST$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestFST,4]
TD_Kosi_nearestDXY<-nearest(DXY_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
DXY$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestDXY,4]
TD_Kosi_nearestAFDabs<-nearest(AFDabs_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
AFDabs$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestAFDabs,4]
TD_Kosi_nearestNIELSEN<-nearest(NIELSEN_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
NIELSEN$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestNIELSEN,4]
TD_Kosi_nearestVARLD<-nearest(VARLD_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
VARLD$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestVARLD,4]
TD_Kosi_nearestFLK<-nearest(FLK_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
FLK$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestFLK,4]
TD_Kosi_nearestDD<-nearest(DD_GRange,TajimasD_Kosi_GRange,ignore.strand=T)
DD$TD_Kosi_nearest<-TajimasD_Kosi[TD_Kosi_nearestDD,4]

FST$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestFST,4]-TajimasD_Kosi[TD_Kosi_nearestFST,4]
NIELSEN$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestNIELSEN,4]-TajimasD_Kosi[TD_Kosi_nearestNIELSEN,4]
DXY$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestDXY,4]-TajimasD_Kosi[TD_Kosi_nearestDXY,4]
VARLD$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestVARLD,4]-TajimasD_Kosi[TD_Kosi_nearestVARLD,4]
FLK$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestFLK,4]-TajimasD_Kosi[TD_Kosi_nearestFLK,4]
AFDabs$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestAFDabs,4]-TajimasD_Kosi[TD_Kosi_nearestAFDabs,4]
DD$TD_diff_nearest_Krom_Kosi<-TajimasD_Krom[TD_Krom_nearestDD,4]-TajimasD_Kosi[TD_Kosi_nearestDD,4]



##############################
#TD and FWH standardized diff#
##############################
require(vegan)
TajimasD_Krom_stand<-decostand(TajimasD_Krom[,4],"standardize")
TajimasD_Kosi_stand<-decostand(TajimasD_Kosi[,4],"standardize")

FST<-cbind(FST,TajimasD_Krom_stand[TD_Krom_nearestFST],TajimasD_Kosi_stand[TD_Kosi_nearestFST])
NIELSEN<-cbind(NIELSEN,TajimasD_Krom_stand[TD_Krom_nearestNIELSEN],TajimasD_Kosi_stand[TD_Kosi_nearestNIELSEN])
DXY<-cbind(DXY,TajimasD_Krom_stand[TD_Krom_nearestDXY],TajimasD_Kosi_stand[TD_Kosi_nearestDXY])
FLK<-cbind(FLK,TajimasD_Krom_stand[TD_Krom_nearestFLK],TajimasD_Kosi_stand[TD_Kosi_nearestFLK])
AFDabs<-cbind(AFDabs,TajimasD_Krom_stand[TD_Krom_nearestAFDabs],TajimasD_Kosi_stand[TD_Kosi_nearestAFDabs])
VARLD<-cbind(VARLD,TajimasD_Krom_stand[TD_Krom_nearestVARLD],TajimasD_Kosi_stand[TD_Kosi_nearestVARLD])
DD<-cbind(DD,TajimasD_Krom_stand[TD_Krom_nearestDD],TajimasD_Kosi_stand[TD_Kosi_nearestDD])

TDdiffstandFST<-TajimasD_Krom_stand[TD_Krom_nearestFST]-TajimasD_Kosi_stand[TD_Kosi_nearestFST]
TDdiffstandDXY<-TajimasD_Krom_stand[TD_Krom_nearestDXY]-TajimasD_Kosi_stand[TD_Kosi_nearestDXY]
TDdiffstandAFDabs<-TajimasD_Krom_stand[TD_Krom_nearestAFDabs]-TajimasD_Kosi_stand[TD_Kosi_nearestAFDabs]
TDdiffstandNIELSEN<-TajimasD_Krom_stand[TD_Krom_nearestNIELSEN]-TajimasD_Kosi_stand[TD_Kosi_nearestNIELSEN]
TDdiffstandVARLD<-TajimasD_Krom_stand[TD_Krom_nearestVARLD]-TajimasD_Kosi_stand[TD_Kosi_nearestVARLD]
TDdiffstandFLK<-TajimasD_Krom_stand[TD_Krom_nearestFLK]-TajimasD_Kosi_stand[TD_Kosi_nearestFLK]
TDdiffstandDD<-TajimasD_Krom_stand[TD_Krom_nearestDD]-TajimasD_Kosi_stand[TD_Kosi_nearestDD]

FST$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandFST
DXY$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandDXY
AFDabs$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandAFDabs
NIELSEN$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandNIELSEN
VARLD$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandVARLD
FLK$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandFLK
DD$TD_diff_nearest_Krom_Kosi_stand<-TDdiffstandDD

names(FST)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")
names(DXY)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")
names(AFDabs)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")
names(NIELSEN)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")
names(VARLD)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")
names(FLK)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")
names(DD)[85:86]<-c("TajimasD_Krom_nearest_stand","TajimasD_Kosi_nearest_stand")

########################################
#replace SNPeff info to window specific#
########################################
#load SNPeff data
SNPeff<-read.csv("KromKosi.ann.table",header=T,sep="\t")

SNPeffAnnsplit<-strsplit(as.character(SNPeff$ANN),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",c(1:3,5))
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffcomp1<-cbind(SNPeff[,1:2],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,na.exclude(SNPeffcomp1))
	}
SNPeffdf2<-SNPeffdf
SPeffdf<-SNPeffdf2
Gene<-gsub("^.*-","",SNPeffdf2$Gene)
SNPeffdf$Gene<-Gene
write.table(SNPeffdf,"SNPeff_KromKosi.table",sep="\t",row.names=F)
SNPeffdf<-read.table("SNPeff_KromKosi.table",sep="\t",header=T)

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ah_ID")
gff<-gff[gff$Type=="gene",]

SNPeffdf3<-merge(SNPeffdf,gff,by.x="Gene",by.y="Ah_ID")
SNPeffdf<-SNPeffdf3
SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]
SNPeffdfMODIFIER_GRange<-GRanges(seqnames=SNPeffdfMODIFIER$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER$POS,end=SNPeffdfMODIFIER$POS))
values(SNPeffdfMODIFIER_GRange)<-SNPeffdfMODIFIER[,1:6]

FstallMODIFIER=mergeByOverlaps(FST_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
FstallMODIFIER_df=as.data.frame(FstallMODIFIER)
FstallMODIFIER_list<-FstallMODIFIER_df[as.character(FstallMODIFIER_df$Ah_ID)==as.character(FstallMODIFIER_df$Gene),]
FST$MODIFIER<-rep("NO",nrow(FST))
FST$MODIFIER<-ifelse(!is.na(match(paste(FST$Ah_ID,FST$Chr,FST$Window_start,FST$Window_end),paste(FstallMODIFIER_list$Gene,FstallMODIFIER_list$FST_GRange.seqnames,FstallMODIFIER_list$FST_GRange.start,FstallMODIFIER_list$FST_GRange.end))),"MODIFIER",FST$MODIFIER)

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
SNPeffdfMODERATE_GRange<-GRanges(seqnames=SNPeffdfMODERATE$CHROM,ranges=IRanges(start=SNPeffdfMODERATE$POS,end=SNPeffdfMODERATE$POS))
values(SNPeffdfMODERATE_GRange)<-SNPeffdfMODERATE[,1:6]

FstallMODERATE=mergeByOverlaps(FST_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
FstallMODERATE_df=as.data.frame(FstallMODERATE)
FstallMODERATE_list<-FstallMODERATE_df[as.character(FstallMODERATE_df$Ah_ID)==as.character(FstallMODERATE_df$Gene),]
FST$MODERATE<-rep("NO",nrow(FST))
FST$MODERATE<-ifelse(!is.na(match(paste(FST$Ah_ID,FST$Chr,FST$Window_start,FST$Window_end),paste(FstallMODERATE_list$Gene,FstallMODERATE_list$FST_GRange.seqnames,FstallMODERATE_list$FST_GRange.start,FstallMODERATE_list$FST_GRange.end))),"MODERATE",FST$MODERATE)

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH_GRange<-GRanges(seqnames=SNPeffdfHIGH$CHROM,ranges=IRanges(start=SNPeffdfHIGH$POS,end=SNPeffdfHIGH$POS))
values(SNPeffdfHIGH_GRange)<-SNPeffdfHIGH[,1:6]

FstallHIGH=mergeByOverlaps(FST_GRange,SNPeffdfHIGH_GRange,type=c("any"))
FstallHIGH_df=as.data.frame(FstallHIGH)
FstallHIGH_list<-FstallHIGH_df[as.character(FstallHIGH_df$Ah_ID)==as.character(FstallHIGH_df$Gene),]
FST$HIGH<-rep("NO",nrow(FST))
FST$HIGH<-ifelse(!is.na(match(paste(FST$Ah_ID,FST$Chr,FST$Window_start,FST$Window_end),paste(FstallHIGH_list$Gene,FstallHIGH_list$FST_GRange.seqnames,FstallHIGH_list$FST_GRange.start,FstallHIGH_list$FST_GRange.end))),"HIGH",FST$HIGH)

DDallMODIFIER=mergeByOverlaps(DD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
DDallMODIFIER_df=as.data.frame(DDallMODIFIER)
DDallMODIFIER_list<-DDallMODIFIER_df[as.character(DDallMODIFIER_df$Ah_ID)==as.character(DDallMODIFIER_df$Gene),]
DD$MODIFIER<-rep("NO",nrow(DD))
DD$MODIFIER<-ifelse(!is.na(match(paste(DD$Ah_ID,DD$Chr,DD$Window_start,DD$Window_end),paste(DDallMODIFIER_list$Gene,DDallMODIFIER_list$DD_GRange.seqnames,DDallMODIFIER_list$DD_GRange.start,DDallMODIFIER_list$DD_GRange.end))),"MODIFIER",DD$MODIFIER)

DDallMODERATE=mergeByOverlaps(DD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
DDallMODERATE_df=as.data.frame(DDallMODERATE)
DDallMODERATE_list<-DDallMODERATE_df[as.character(DDallMODERATE_df$Ah_ID)==as.character(DDallMODERATE_df$Gene),]
DD$MODERATE<-rep("NO",nrow(DD))
DD$MODERATE<-ifelse(!is.na(match(paste(DD$Ah_ID,DD$Chr,DD$Window_start,DD$Window_end),paste(DDallMODERATE_list$Gene,DDallMODERATE_list$DD_GRange.seqnames,DDallMODERATE_list$DD_GRange.start,DDallMODERATE_list$DD_GRange.end))),"MODERATE",DD$MODERATE)

DDallHIGH=mergeByOverlaps(DD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
DDallHIGH_df=as.data.frame(DDallHIGH)
DDallHIGH_list<-DDallHIGH_df[as.character(DDallHIGH_df$Ah_ID)==as.character(DDallHIGH_df$Gene),]
DD$HIGH<-rep("NO",nrow(DD))
DD$HIGH<-ifelse(!is.na(match(paste(DD$Ah_ID,DD$Chr,DD$Window_start,DD$Window_end),paste(DDallHIGH_list$Gene,DDallHIGH_list$DD_GRange.seqnames,DDallHIGH_list$DD_GRange.start,DDallHIGH_list$DD_GRange.end))),"HIGH",DD$HIGH)


NielsenallMODIFIER=mergeByOverlaps(NIELSEN_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
NielsenallMODIFIER_df=as.data.frame(NielsenallMODIFIER)
NielsenallMODIFIER_list<-NielsenallMODIFIER_df[as.character(NielsenallMODIFIER_df$Ah_ID)==as.character(NielsenallMODIFIER_df$Gene),]
NIELSEN$MODIFIER<-rep("NO",nrow(NIELSEN))
NIELSEN$MODIFIER<-ifelse(!is.na(match(paste(NIELSEN$Ah_ID,NIELSEN$Chr,NIELSEN$Window_start,NIELSEN$Window_end),paste(NielsenallMODIFIER_list$Gene,NielsenallMODIFIER_list$NIELSEN_GRange.seqnames,NielsenallMODIFIER_list$NIELSEN_GRange.start,NielsenallMODIFIER_list$NIELSEN_GRange.end))),"MODIFIER",NIELSEN$MODIFIER)

FlkallMODIFIER=mergeByOverlaps(FLK_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
FlkallMODIFIER_df=as.data.frame(FlkallMODIFIER)
FlkallMODIFIER_list<-FlkallMODIFIER_df[as.character(FlkallMODIFIER_df$Ah_ID)==as.character(FlkallMODIFIER_df$Gene),]
FLK$MODIFIER<-rep("NO",nrow(FLK))
FLK$MODIFIER<-ifelse(!is.na(match(paste(FLK$Ah_ID,FLK$Chr,FLK$Window_start,FLK$Window_end),paste(FlkallMODIFIER_list$Gene,FlkallMODIFIER_list$FLK_GRange.seqnames,FlkallMODIFIER_list$FLK_GRange.start,FlkallMODIFIER_list$FLK_GRange.end))),"MODIFIER",FLK$MODIFIER)

DxyallMODIFIER=mergeByOverlaps(DXY_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
DxyallMODIFIER_df=as.data.frame(DxyallMODIFIER)
DxyallMODIFIER_list<-DxyallMODIFIER_df[as.character(DxyallMODIFIER_df$Ah_ID)==as.character(DxyallMODIFIER_df$Gene),]
DXY$MODIFIER<-rep("NO",nrow(DXY))
DXY$MODIFIER<-ifelse(!is.na(match(paste(DXY$Ah_ID,DXY$Chr,DXY$Window_start,DXY$Window_end),paste(DxyallMODIFIER_list$Gene,DxyallMODIFIER_list$DXY_GRange.seqnames,DxyallMODIFIER_list$DXY_GRange.start,DxyallMODIFIER_list$DXY_GRange.end))),"MODIFIER",DXY$MODIFIER)

VARLDallMODIFIER=mergeByOverlaps(VARLD_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
VARLDallMODIFIER_df=as.data.frame(VARLDallMODIFIER)
VARLDallMODIFIER_list<-VARLDallMODIFIER_df[as.character(VARLDallMODIFIER_df$Ah_ID)==as.character(VARLDallMODIFIER_df$Gene),]
VARLD$MODIFIER<-rep("NO",nrow(VARLD))
VARLD$MODIFIER<-ifelse(!is.na(match(paste(VARLD$Ah_ID,VARLD$Chr,VARLD$Window_start,VARLD$Window_end),paste(VARLDallMODIFIER_list$Gene,VARLDallMODIFIER_list$VARLD_GRange.seqnames,VARLDallMODIFIER_list$VARLD_GRange.start,VARLDallMODIFIER_list$VARLD_GRange.end))),"MODIFIER",VARLD$MODIFIER)

AFDabsallMODIFIER=mergeByOverlaps(AFDabs_GRange,SNPeffdfMODIFIER_GRange,type=c("any"))
AFDabsallMODIFIER_df=as.data.frame(AFDabsallMODIFIER)
AFDabsallMODIFIER_list<-AFDabsallMODIFIER_df[as.character(AFDabsallMODIFIER_df$Ah_ID)==as.character(AFDabsallMODIFIER_df$Gene),]
AFDabs$MODIFIER<-rep("NO",nrow(AFDabs))
AFDabs$MODIFIER<-ifelse(!is.na(match(paste(AFDabs$Ah_ID,AFDabs$Chr,AFDabs$Window_start,AFDabs$Window_end),paste(AFDabsallMODIFIER_list$Gene,AFDabsallMODIFIER_list$AFDabs_GRange.seqnames,AFDabsallMODIFIER_list$AFDabs_GRange.start,AFDabsallMODIFIER_list$AFDabs_GRange.end))),"MODIFIER",AFDabs$MODIFIER)


NielsenallMODERATE=mergeByOverlaps(NIELSEN_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
NielsenallMODERATE_df=as.data.frame(NielsenallMODERATE)
NielsenallMODERATE_list<-NielsenallMODERATE_df[as.character(NielsenallMODERATE_df$Ah_ID)==as.character(NielsenallMODERATE_df$Gene),]
NIELSEN$MODERATE<-rep("NO",nrow(NIELSEN))
NIELSEN$MODERATE<-ifelse(match(paste(NIELSEN$Ah_ID,NIELSEN$Chr,NIELSEN$Window_start,NIELSEN$Window_end),paste(NielsenallMODERATE_list$Gene,NielsenallMODERATE_list$NIELSEN_GRange.seqnames,NielsenallMODERATE_list$NIELSEN_GRange.start,NielsenallMODERATE_list$NIELSEN_GRange.end)),"MODERATE",NIELSEN$MODERATE)

FlkallMODERATE=mergeByOverlaps(FLK_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
FlkallMODERATE_df=as.data.frame(FlkallMODERATE)
FlkallMODERATE_list<-FlkallMODERATE_df[as.character(FlkallMODERATE_df$Ah_ID)==as.character(FlkallMODERATE_df$Gene),]
FLK$MODERATE<-rep("NO",nrow(FLK))
FLK$MODERATE<-ifelse(!is.na(match(paste(FLK$Ah_ID,FLK$Chr,FLK$Window_start,FLK$Window_end),paste(FlkallMODERATE_list$Gene,FlkallMODERATE_list$FLK_GRange.seqnames,FlkallMODERATE_list$FLK_GRange.start,FlkallMODERATE_list$FLK_GRange.end))),"MODERATE",FLK$MODERATE)

DxyallMODERATE=mergeByOverlaps(DXY_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
DxyallMODERATE_df=as.data.frame(DxyallMODERATE)
DxyallMODERATE_list<-DxyallMODERATE_df[as.character(DxyallMODERATE_df$Ah_ID)==as.character(DxyallMODERATE_df$Gene),]
DXY$MODERATE<-rep("NO",nrow(DXY))
DXY$MODERATE<-ifelse(!is.na(match(paste(DXY$Ah_ID,DXY$Chr,DXY$Window_start,DXY$Window_end),paste(DxyallMODERATE_list$Gene,DxyallMODERATE_list$DXY_GRange.seqnames,DxyallMODERATE_list$DXY_GRange.start,DxyallMODERATE_list$DXY_GRange.end))),"MODERATE",DXY$MODERATE)

VARLDallMODERATE=mergeByOverlaps(VARLD_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
VARLDallMODERATE_df=as.data.frame(VARLDallMODERATE)
VARLDallMODERATE_list<-VARLDallMODERATE_df[as.character(VARLDallMODERATE_df$Ah_ID)==as.character(VARLDallMODERATE_df$Gene),]
VARLD$MODERATE<-rep("NO",nrow(VARLD))
VARLD$MODERATE<-ifelse(!is.na(match(paste(VARLD$Ah_ID,VARLD$Chr,VARLD$Window_start,VARLD$Window_end),paste(VARLDallMODERATE_list$Gene,VARLDallMODERATE_list$VARLD_GRange.seqnames,VARLDallMODERATE_list$VARLD_GRange.start,VARLDallMODERATE_list$VARLD_GRange.end))),"MODERATE",VARLD$MODERATE)

AFDabsallMODERATE=mergeByOverlaps(AFDabs_GRange,SNPeffdfMODERATE_GRange,type=c("any"))
AFDabsallMODERATE_df=as.data.frame(AFDabsallMODERATE)
AFDabsallMODERATE_list<-AFDabsallMODERATE_df[as.character(AFDabsallMODERATE_df$Ah_ID)==as.character(AFDabsallMODERATE_df$Gene),]
AFDabs$MODERATE<-rep("NO",nrow(AFDabs))
AFDabs$MODERATE<-ifelse(!is.na(match(paste(AFDabs$Ah_ID,AFDabs$Chr,AFDabs$Window_start,AFDabs$Window_end),paste(AFDabsallMODERATE_list$Gene,AFDabsallMODERATE_list$AFDabs_GRange.seqnames,AFDabsallMODERATE_list$AFDabs_GRange.start,AFDabsallMODERATE_list$AFDabs_GRange.end))),"MODERATE",AFDabs$MODERATE)


NielsenallHIGH=mergeByOverlaps(NIELSEN_GRange,SNPeffdfHIGH_GRange,type=c("any"))
NielsenallHIGH_df=as.data.frame(NielsenallHIGH)
NielsenallHIGH_list<-NielsenallHIGH_df[as.character(NielsenallHIGH_df$Ah_ID)==as.character(NielsenallHIGH_df$Gene),]
NIELSEN$HIGH<-rep("NO",nrow(NIELSEN))
NIELSEN$HIGH<-ifelse(!is.na(match(paste(NIELSEN$Ah_ID,NIELSEN$Chr,NIELSEN$Window_start,NIELSEN$Window_end),paste(NielsenallHIGH_list$Gene,NielsenallHIGH_list$NIELSEN_GRange.seqnames,NielsenallHIGH_list$NIELSEN_GRange.start,NielsenallHIGH_list$NIELSEN_GRange.end))),"HIGH",NIELSEN$HIGH)

FlkallHIGH=mergeByOverlaps(FLK_GRange,SNPeffdfHIGH_GRange,type=c("any"))
FlkallHIGH_df=as.data.frame(FlkallHIGH)
FlkallHIGH_list<-FlkallHIGH_df[as.character(FlkallHIGH_df$Ah_ID)==as.character(FlkallHIGH_df$Gene),]
FLK$HIGH<-rep("NO",nrow(FLK))
FLK$HIGH<-ifelse(!is.na(match(paste(FLK$Ah_ID,FLK$Chr,FLK$Window_start,FLK$Window_end),paste(FlkallHIGH_list$Gene,FlkallHIGH_list$FLK_GRange.seqnames,FlkallHIGH_list$FLK_GRange.start,FlkallHIGH_list$FLK_GRange.end))),"HIGH",FLK$HIGH)

DxyallHIGH=mergeByOverlaps(DXY_GRange,SNPeffdfHIGH_GRange,type=c("any"))
DxyallHIGH_df=as.data.frame(DxyallHIGH)
DxyallHIGH_list<-DxyallHIGH_df[as.character(DxyallHIGH_df$Ah_ID)==as.character(DxyallHIGH_df$Gene),]
DXY$HIGH<-rep("NO",nrow(DXY))
DXY$HIGH<-ifelse(!is.na(match(paste(DXY$Ah_ID,DXY$Chr,DXY$Window_start,DXY$Window_end),paste(DxyallHIGH_list$Gene,DxyallHIGH_list$DXY_GRange.seqnames,DxyallHIGH_list$DXY_GRange.start,DxyallHIGH_list$DXY_GRange.end))),"HIGH",DXY$HIGH)

VARLDallHIGH=mergeByOverlaps(VARLD_GRange,SNPeffdfHIGH_GRange,type=c("any"))
VARLDallHIGH_df=as.data.frame(VARLDallHIGH)
VARLDallHIGH_list<-VARLDallHIGH_df[as.character(VARLDallHIGH_df$Ah_ID)==as.character(VARLDallHIGH_df$Gene),]
VARLD$HIGH<-rep("NO",nrow(VARLD))
VARLD$HIGH<-ifelse(!is.na(match(paste(VARLD$Ah_ID,VARLD$Chr,VARLD$Window_start,VARLD$Window_end),paste(VARLDallHIGH_list$Gene,VARLDallHIGH_list$VARLD_GRange.seqnames,VARLDallHIGH_list$VARLD_GRange.start,VARLDallHIGH_list$VARLD_GRange.end))),"HIGH",VARLD$HIGH)

AFDabsallHIGH=mergeByOverlaps(AFDabs_GRange,SNPeffdfHIGH_GRange,type=c("any"))
AFDabsallHIGH_df=as.data.frame(AFDabsallHIGH)
AFDabsallHIGH_list<-AFDabsallHIGH_df[as.character(AFDabsallHIGH_df$Ah_ID)==as.character(AFDabsallHIGH_df$Gene),]
AFDabs$HIGH<-rep("NO",nrow(AFDabs))
AFDabs$HIGH<-ifelse(!is.na(match(paste(AFDabs$Ah_ID,AFDabs$Chr,AFDabs$Window_start,AFDabs$Window_end),paste(AFDabsallHIGH_list$Gene,AFDabsallHIGH_list$AFDabs_GRange.seqnames,AFDabsallHIGH_list$AFDabs_GRange.start,AFDabsallHIGH_list$AFDabs_GRange.end))),"HIGH",AFDabs$HIGH)


#############################################################################
#add SNPeff info with info for whole gene#
#############################################################################
data1<-read.table("Kromnew.table",header=TRUE,sep="\t")
data2<-read.table("Kosinew.table",header=TRUE,sep="\t")

AFdata<-read.table("HalleriKromKosiGS_new.csv",header=TRUE,sep="\t")
test<-as.character(AFdata$CHROM)
AFdata[,1]<-test
AFdata$AF<-AFdata$AC/AFdata$AN
AFdata$AF.1<-AFdata$AC.1/AFdata$AN.1
AFdata2<-data.frame(AFdata[,c(1,2,9)],AFdata[,10])
names(AFdata2)<-c("CHROM","POS","AF_Krom","AF_Kosi")

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH2<-merge(SNPeffdfHIGH,AFdata2)

FSTHighKrom<-vector(mode="character",length=nrow(FST))
for (i in 1:nrow(FST)) {
	if(FST$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		FSTHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(FST$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FSTHighKrom[i]<-"NO"
		}
	}
FSTHighKosi<-vector(mode="character",length=nrow(FST))
for (i in 1:nrow(FST)) {
	if(FST$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		FSTHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(FST$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FSTHighKosi[i]<-"NO"
		}
	}

test<-FST
FST<-data.frame(append(test,list(AF_HIGH_Krom_gene=FSTHighKrom),after=match("HIGH",names(test))))
FST<-data.frame(append(FST,list(AF_HIGH_Kosi_gene=FSTHighKosi),after=match("AF_HIGH_Krom_gene",names(FST))))

DDHighKrom<-vector(mode="character",length=nrow(DD))
for (i in 1:nrow(DD)) {
	if(DD$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		DDHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(DD$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DDHighKrom[i]<-"NO"
		}
	}
DDHighKosi<-vector(mode="character",length=nrow(DD))
for (i in 1:nrow(DD)) {
	if(DD$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		DDHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(DD$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DDHighKosi[i]<-"NO"
		}
	}

test<-DD
DD<-data.frame(append(test,list(AF_HIGH_Krom_gene=DDHighKrom),after=match("HIGH",names(test))))
DD<-data.frame(append(DD,list(AF_HIGH_Kosi_gene=DDHighKosi),after=match("AF_HIGH_Krom_gene",names(DD))))


NIELSENHighKrom<-vector(mode="character",length=nrow(NIELSEN))
for (i in 1:nrow(NIELSEN)) {
	if(NIELSEN$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		NIELSENHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(NIELSEN$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		NIELSENHighKrom[i]<-"NO"
		}
	}
NIELSENHighKosi<-vector(mode="character",length=nrow(NIELSEN))
for (i in 1:nrow(NIELSEN)) {
	if(NIELSEN$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		NIELSENHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(NIELSEN$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		NIELSENHighKosi[i]<-"NO"
		}
	}

test<-NIELSEN
NIELSEN<-data.frame(append(test,list(AF_HIGH_Krom_gene=NIELSENHighKrom),after=match("HIGH",names(test))))
NIELSEN<-data.frame(append(NIELSEN,list(AF_HIGH_Kosi_gene=NIELSENHighKosi),after=match("AF_HIGH_Krom_gene",names(NIELSEN))))

DXYHighKrom<-vector(mode="character",length=nrow(DXY))
for (i in 1:nrow(DXY)) {
	if(DXY$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		DXYHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(DXY$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DXYHighKrom[i]<-"NO"
		}
	}
DXYHighKosi<-vector(mode="character",length=nrow(DXY))
for (i in 1:nrow(DXY)) {
	if(DXY$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		DXYHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(DXY$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		DXYHighKosi[i]<-"NO"
		}
	}

test<-DXY
DXY<-data.frame(append(test,list(AF_HIGH_Krom_gene=DXYHighKrom),after=match("HIGH",names(test))))
DXY<-data.frame(append(DXY,list(AF_HIGH_Kosi_gene=DXYHighKosi),after=match("AF_HIGH_Krom_gene",names(DXY))))

AFDabsHighKrom<-vector(mode="character",length=nrow(AFDabs))
for (i in 1:nrow(AFDabs)) {
	if(AFDabs$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		AFDabsHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(AFDabs$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		AFDabsHighKrom[i]<-"NO"
		}
	}
AFDabsHighKosi<-vector(mode="character",length=nrow(AFDabs))
for (i in 1:nrow(AFDabs)) {
	if(AFDabs$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		AFDabsHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(AFDabs$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		AFDabsHighKosi[i]<-"NO"
		}
	}

test<-AFDabs
AFDabs<-data.frame(append(test,list(AF_HIGH_Krom_gene=AFDabsHighKrom),after=match("HIGH",names(test))))
AFDabs<-data.frame(append(AFDabs,list(AF_HIGH_Kosi_gene=AFDabsHighKosi),after=match("AF_HIGH_Krom_gene",names(AFDabs))))

FLKHighKrom<-vector(mode="character",length=nrow(FLK))
for (i in 1:nrow(FLK)) {
	if(FLK$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		FLKHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(FLK$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FLKHighKrom[i]<-"NO"
		}
	}
FLKHighKosi<-vector(mode="character",length=nrow(FLK))
for (i in 1:nrow(FLK)) {
	if(FLK$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		FLKHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(FLK$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		FLKHighKosi[i]<-"NO"
		}
	}

test<-FLK
FLK<-data.frame(append(test,list(AF_HIGH_Krom_gene=FLKHighKrom),after=match("HIGH",names(test))))
FLK<-data.frame(append(FLK,list(AF_HIGH_Kosi_gene=FLKHighKosi),after=match("AF_HIGH_Krom_gene",names(FLK))))


VARLDHighKrom<-vector(mode="character",length=nrow(VARLD))
for (i in 1:nrow(VARLD)) {
	if(VARLD$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		VARLDHighKrom[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Krom[match(VARLD$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		VARLDHighKrom[i]<-"NO"
		}
	}
VARLDHighKosi<-vector(mode="character",length=nrow(VARLD))
for (i in 1:nrow(VARLD)) {
	if(VARLD$Ah_ID[i]%in%SNPeffdfHIGH2$Gene) {
		VARLDHighKosi[i]<-max(na.exclude(SNPeffdfHIGH2$AF_Kosi[match(VARLD$Ah_ID[i],SNPeffdfHIGH2$Gene)]))
		} else {
		VARLDHighKosi[i]<-"NO"
		}
	}

test<-VARLD
VARLD<-data.frame(append(test,list(AF_HIGH_Krom_gene=VARLDHighKrom),after=match("HIGH",names(test))))
VARLD<-data.frame(append(VARLD,list(AF_HIGH_Kosi_gene=VARLDHighKosi),after=match("AF_HIGH_Krom_gene",names(VARLD))))


#############################################
#add Rank info#
#############################################
FST$Rank_Fst<-rank(-as.numeric(as.character(FST$Fst)),ties.method="min")
DD$Rank_DD<-rank(as.numeric(as.character(DD$DD)),ties.method="min")
NIELSEN$Rank_Nielsen<-rank(-as.numeric(as.character(NIELSEN$Nielsen)),ties.method="min")
DXY$Rank_Dxy<-rank(-as.numeric(as.character(DXY$Dxy)),ties.method="min")
AFDabs$Rank_AFDabs<-rank(-as.numeric(as.character(AFDabs$AFDabs)),ties.method="min")
FLK$Rank_Flk<-rank(-as.numeric(as.character(FLK$Flk)),ties.method="min")
VARLD$Rank_Varld<-rank(-as.numeric(as.character(VARLD$VarLD)),ties.method="min")

######################################################
#add Sweedinfo#
######################################################
filelist = list.files(pattern = "SweeD_Report.Kosi_SweeD_allSNPs*") 
require(gtools)
filelist <- mixedsort(filelist)
scaff1<-read.table("SweeD_Report.Kosi_SweeD_allSNPs.chr1",sep="\t",header=T)
datafr = do.call(rbind,lapply(filelist[-1],function(fn)read.table(fn,header=F,sep="\t",colClasses =c("numeric","numeric","numeric") ))) #merge all count files
chr<-strsplit(filelist,"\\.")
chr1<-lapply(chr, tail, n = 1L)
chr2<-data.frame(matrix(unlist(chr1),nrow=length(chr1),byrow=T))[,1]
#datanrow = do.call(rbind,lapply(filelist,function(fn)nrow(read.table(fn,header=FALSE, sep="\t"))))
#all 37500
Scaffold<-rep(as.character(chr2),each=37500)
colnames(datafr)<-colnames(scaff1)
datafr2<-rbind(scaff1,datafr)
dataSweed_Kosi<-data.frame(Scaffold,datafr2)
names(dataSweed_Kosi)<-c("Scaffold","Position","Likelihood","Alpha")

filelist = list.files(pattern = "SweeD_Report.Krom_SweeD_allSNPs*") 
require(gtools)
filelist <- mixedsort(filelist)
scaff1<-read.table("SweeD_Report.Krom_SweeD_allSNPs.chr1",sep="\t",header=T)
datafr = do.call(rbind,lapply(filelist[-1],function(fn)read.table(fn,header=F,sep="\t",colClasses =c("numeric","numeric","numeric") ))) #merge all count files
chr<-strsplit(filelist,"\\.")
chr1<-lapply(chr, tail, n = 1L)
chr2<-data.frame(matrix(unlist(chr1),nrow=length(chr1),byrow=T))[,1]
#datanrow = do.call(rbind,lapply(filelist,function(fn)nrow(read.table(fn,header=FALSE, sep="\t"))))
#all 37500
Scaffold<-rep(as.character(chr2),each=37500)
colnames(datafr)<-colnames(scaff1)
datafr2<-rbind(scaff1,datafr)
dataSweed_Krom<-data.frame(Scaffold,datafr2)
names(dataSweed_Krom)<-c("Scaffold","Position","Likelihood","Alpha")

require(GenomicRanges)
Sweed_Krom_GR<-GRanges(seqnames=tolower(dataSweed_Krom$Scaffold),ranges=IRanges(start=dataSweed_Krom$Position,end=dataSweed_Krom$Position))
values(Sweed_Krom_GR)<-dataSweed_Krom[,3:4]


Sweed_Krom_nearestFST<-nearest(FST_GRange,Sweed_Krom_GR,ignore.strand=T)
FST$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestFST,3]
FST$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestFST,4]
Sweed_Krom_nearestDXY<-nearest(DXY_GRange,Sweed_Krom_GR,ignore.strand=T)
DXY$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestDXY,3]
DXY$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestDXY,4]
Sweed_Krom_nearestNIELSEN<-nearest(NIELSEN_GRange,Sweed_Krom_GR,ignore.strand=T)
NIELSEN$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestNIELSEN,3]
NIELSEN$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestNIELSEN,4]
Sweed_Krom_nearestFLK<-nearest(FLK_GRange,Sweed_Krom_GR,ignore.strand=T)
FLK$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestFLK,3]
FLK$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestFLK,4]
Sweed_Krom_nearestVARLD<-nearest(VARLD_GRange,Sweed_Krom_GR,ignore.strand=T)
VARLD$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestVARLD,3]
VARLD$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestVARLD,4]
Sweed_Krom_nearestAFDabs<-nearest(AFDabs_GRange,Sweed_Krom_GR,ignore.strand=T)
AFDabs$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestAFDabs,3]
AFDabs$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestAFDabs,4]
Sweed_Krom_nearestDD<-nearest(DD_GRange,Sweed_Krom_GR,ignore.strand=T)
DD$Sweed_Krom_nearest_Likelihood<-dataSweed_Krom[Sweed_Krom_nearestDD,3]
DD$Sweed_Krom_nearest_alpha<-dataSweed_Krom[Sweed_Krom_nearestDD,4]

Sweed_Kosi_GR<-GRanges(seqnames=tolower(dataSweed_Kosi$Scaffold),ranges=IRanges(start=dataSweed_Kosi$Position,end=dataSweed_Kosi$Position))
values(Sweed_Kosi_GR)<-dataSweed_Kosi[,3:4]

Sweed_Kosi_nearestFST<-nearest(FST_GRange,Sweed_Kosi_GR,ignore.strand=T)
FST$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestFST,3]
FST$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestFST,4]
Sweed_Kosi_nearestDXY<-nearest(DXY_GRange,Sweed_Kosi_GR,ignore.strand=T)
DXY$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestDXY,3]
DXY$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestDXY,4]
Sweed_Kosi_nearestNIELSEN<-nearest(NIELSEN_GRange,Sweed_Kosi_GR,ignore.strand=T)
NIELSEN$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestNIELSEN,3]
NIELSEN$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestNIELSEN,4]
Sweed_Kosi_nearestFLK<-nearest(FLK_GRange,Sweed_Kosi_GR,ignore.strand=T)
FLK$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestFLK,3]
FLK$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestFLK,4]
Sweed_Kosi_nearestVARLD<-nearest(VARLD_GRange,Sweed_Kosi_GR,ignore.strand=T)
VARLD$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestVARLD,3]
VARLD$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestVARLD,4]
Sweed_Kosi_nearestAFDabs<-nearest(AFDabs_GRange,Sweed_Kosi_GR,ignore.strand=T)
AFDabs$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestAFDabs,3]
AFDabs$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestAFDabs,4]
Sweed_Kosi_nearestDD<-nearest(DD_GRange,Sweed_Kosi_GR,ignore.strand=T)
DD$Sweed_Kosi_nearest_Likelihood<-dataSweed_Kosi[Sweed_Kosi_nearestDD,3]
DD$Sweed_Kosi_nearest_alpha<-dataSweed_Kosi[Sweed_Kosi_nearestDD,4]


#################################################
#add info in how many metrics among 0.1% outlier#
#################################################
Quantiles<-read.xlsx2("Quantiles_KromKosi.xlsx",1,header=F)
QuanileFst<-as.numeric(as.character(Quantiles$X11[2]))
QuanileAFDabs<-as.numeric(as.character(Quantiles$X11[4]))
QuanileNielsen<-as.numeric(as.character(Quantiles$X11[6]))
QuanileDxy<-as.numeric(as.character(Quantiles$X11[8]))
QuanileVarLD<-as.numeric(as.character(Quantiles$X11[10]))
QuanileFlk<-as.numeric(as.character(Quantiles$X11[12]))
QuanileDD<-as.numeric(as.character(Quantiles$X3[14]))

countFST<-vector(mode="numeric",length=nrow(FST))
for (i in 1:nrow(FST))
	{countFST[i]=0
	if (as.numeric(as.character(FST$Nielsen[i]))>=QuanileNielsen)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$AFDabs[i]))>=QuanileAFDabs)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$Dxy[i]))>=QuanileDxy)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$VarLD[i]))>=QuanileVarLD)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$Flk[i]))>=QuanileFlk&as.numeric(as.character(FST$Flk_pvalue[i]))<=0.05)
		{countFST[i]=countFST[i]+1
		}
	if (as.numeric(as.character(FST$DD[i]))<=QuanileDD)
		{countFST[i]=countFST[i]+1
		}
	}
FST$Number_of_metrics_0.1percent<-countFST


countNIELSEN<-vector(mode="numeric",length=nrow(NIELSEN))
for (i in 1:nrow(NIELSEN))
	{countNIELSEN[i]=0
	if (as.numeric(as.character(NIELSEN$Fst[i]))>=QuanileFst)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$AFDabs[i]))>=QuanileAFDabs)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$Dxy[i]))>=QuanileDxy)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$VarLD[i]))>=QuanileVarLD)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$Flk[i]))>=QuanileFlk&as.numeric(as.character(NIELSEN$Flk_pvalue[i]))<=0.05)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	if (as.numeric(as.character(NIELSEN$DD[i]))<=QuanileDD)
		{countNIELSEN[i]=countNIELSEN[i]+1
		}
	}
NIELSEN$Number_of_metrics_0.1percent<-countNIELSEN

countDXY<-vector(mode="numeric",length=nrow(DXY))
for (i in 1:nrow(DXY))
	{countDXY[i]=0
	if (as.numeric(as.character(DXY$Fst[i]))>=QuanileFst)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$AFDabs[i]))>=QuanileAFDabs)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$Nielsen[i]))>=QuanileNielsen)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$VarLD[i]))>=QuanileVarLD)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$Flk[i]))>=QuanileFlk&as.numeric(as.character(DXY$Flk_pvalue[i]))<=0.05)
		{countDXY[i]=countDXY[i]+1
		}
	if (as.numeric(as.character(DXY$DD[i]))<=QuanileDD)
		{countDXY[i]=countDXY[i]+1
		}
	}
DXY$Number_of_metrics_0.1percent<-countDXY

countVARLD<-vector(mode="numeric",length=nrow(VARLD))
for (i in 1:nrow(VARLD))
	{countVARLD[i]=0
	if (as.numeric(as.character(VARLD$Fst[i]))>=QuanileFst)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$AFDabs[i]))>=QuanileAFDabs)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$Nielsen[i]))>=QuanileNielsen)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$Dxy[i]))>=QuanileDxy)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$Flk[i]))>=QuanileFlk&as.numeric(as.character(VARLD$Flk_pvalue[i]))<=0.05)
		{countVARLD[i]=countVARLD[i]+1
		}
	if (as.numeric(as.character(VARLD$DD[i]))<=QuanileDD)
		{countVARLD[i]=countVARLD[i]+1
		}
	}
VARLD$Number_of_metrics_0.1percent<-countVARLD

countFLK<-vector(mode="numeric",length=nrow(FLK))
for (i in 1:nrow(FLK))
	{countFLK[i]=0
	if (as.numeric(as.character(FLK$Fst[i]))>=QuanileFst)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$AFDabs[i]))>=QuanileAFDabs)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$Nielsen[i]))>=QuanileNielsen)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$Dxy[i]))>=QuanileDxy)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$VarLD[i]))>=QuanileVarLD)
		{countFLK[i]=countFLK[i]+1
		}
	if (as.numeric(as.character(FLK$DD[i]))<=QuanileDD)
		{countFLK[i]=countFLK[i]+1
		}
	}
FLK$Number_of_metrics_0.1percent<-countFLK

countAFDabs<-vector(mode="numeric",length=nrow(AFDabs))
for (i in 1:nrow(AFDabs))
	{countAFDabs[i]=0
	if (as.numeric(as.character(AFDabs$Fst[i]))>=QuanileFst)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$Flk[i]))>=QuanileFlk&as.numeric(as.character(AFDabs$Flk_pvalue[i]))<=0.05)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$Nielsen[i]))>=QuanileNielsen)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$Dxy[i]))>=QuanileDxy)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$VarLD[i]))>=QuanileVarLD)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	if (as.numeric(as.character(AFDabs$DD[i]))<=QuanileDD)
		{countAFDabs[i]=countAFDabs[i]+1
		}
	}
AFDabs$Number_of_metrics_0.1percent<-countAFDabs

countDD<-vector(mode="numeric",length=nrow(DD))
for (i in 1:nrow(DD))
	{countDD[i]=0
	if (as.numeric(as.character(DD$Fst[i]))>=QuanileFst)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$Flk[i]))>=QuanileFlk&as.numeric(as.character(DD$Flk_pvalue[i]))<=0.05)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$Nielsen[i]))>=QuanileNielsen)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$Dxy[i]))>=QuanileDxy)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$VarLD[i]))>=QuanileVarLD)
		{countDD[i]=countDD[i]+1
		}
	if (as.numeric(as.character(DD$AFDabs[i]))>=QuanileAFDabs)
		{countDD[i]=countDD[i]+1
		}
	}
DD$Number_of_metrics_0.1percent<-countDD


#################################################################
#add genes with high impact SNPs and allele frequency difference#
#################################################################
#SNPeffdfHIGH2<-merge(SNPeffdfHIGH,AFdata2)

HIGHlist<-SNPeffdfHIGH2[abs(SNPeffdfHIGH2$AF_Krom-SNPeffdfHIGH2$AF_Kosi)>0.4,]
Krom_UG<-read.table("KromUG.table",header=T,sep="\t")
Kosi_UG<-read.table("KosiUG.table",header=T,sep="\t")
Krom_UG$AF<-Krom_UG$AC/Krom_UG$AN
Kosi_UG$AF<-Kosi_UG$AC/Kosi_UG$AN

KK_UG<-cbind(Krom_UG[,c(1,2,5)],Kosi_UG[,c(5)])
names(KK_UG)[4]<-"AF.1"

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

#add Mapman categories
Mapman<-read.table("Mapman_for_merge.table",sep="\t",header=T,fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
Mapman<-Mapman[Mapman$TYPE=="T",]
Mapman_df<-unique(Mapman[,c(3,9)])


#add Lists
All_lists<-read.table("All_lists_Cu_project.table",sep="\t",header=T,fill=T,quote="")

## Load GFF files and add promotor region to genes

### Step 1: import genes from GFF files

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ah_ID")
#gff2<-gff
gff<-gff[gff$Type=="gene",]
#gff$Start[gff$Strand=="+"]<-gff$Start[gff$Strand=="+"]-2000
#gff$End[gff$Strand=="-"]<-gff$End[gff$Strand=="-"]+2000
gff_GRange<-GRanges(seqnames=tolower(gff$Scaffold),ranges=IRanges(start=gff$Start,end=gff$End))
AhalGenes_plusPromot <- gff_GRange
values(AhalGenes_plusPromot)<-gff$Ah_ID
### Create a file with the list of Thaliana genes and descriptions

HIGHlist_hallerigenes_df_w_orthogroup= merge(HIGHlist,OG_Ahalleri, by.x="Gene",by.y="Ah_ID", all.x=TRUE)

HIGHlist_hallerigenes_df_w_MM<-merge(HIGHlist_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)

colnames(HIGHlist_hallerigenes_df_w_MM)[26]<-"Mapman_category"

HIGHlist_hallerigenes_df_w_MM_lists<-merge(HIGHlist_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")

#change order

HIGHlist_hallerigenes_df_w_MM_lists<-HIGHlist_hallerigenes_df_w_MM_lists[,c(2:18,1,19:69)]

HighSNPeff<-HIGHlist_hallerigenes_df_w_MM_lists
diff<-ifelse(paste(HighSNPeff$CHROM,HighSNPeff$POS,sep="_")%in%paste(KK_UG$CHROM,KK_UG$POS,sep="_"),"YES","NO")

summarydiff<-ifelse(diff=="YES",(KK_UG$AF[match(paste(HighSNPeff$CHROM,HighSNPeff$POS,sep="_"),paste(KK_UG$CHROM,KK_UG$POS,sep="_"))])-(KK_UG$AF.1[match(paste(HighSNPeff$CHROM,HighSNPeff$POS,sep="_"),paste(KK_UG$CHROM,KK_UG$POS,sep="_"))]),"NO")
HighSNPeffdiff<-cbind(HighSNPeff,HighSNPeff$AF_Krom-HighSNPeff$AF_Kosi,summarydiff)
names(HighSNPeffdiff)[70:71]<-c("AF_Krom-AF_Kosi","UG_AF_diff")
HighSNPeffdiff<-HighSNPeffdiff[,c(1:17,70:71,18:69)]

write.table(HighSNPeffdiff,"GenesHIGHeffect_above40percentAFdiff_withUG.txt", sep="\t", row.names=F,quote=F)
write.xlsx(HighSNPeffdiff,"GenesHIGHeffect_above40percentAFdiff_withUG.xlsx",row.names=F,sheetName="KromKosi")


write.xlsx2(FST,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="Fst",col.names=TRUE,row.names=FALSE)
write.xlsx2(DXY,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="Dxy",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(AFDabs,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="AFDabs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(NIELSEN,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="Nielsen",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(VARLD,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="VarLD",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(FLK,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="Flk",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(DD,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="DD",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx2(HighSNPeffdiff,"Genes_01percent_KromKosi_refhalleri.xlsx",sheetName="Highimpact",col.names=TRUE,row.names=FALSE,append=T)

write.csv(FST,"FST_KromKosi_01percent.csv",row.names=F)
write.csv(DXY,"DXY_KromKosi_01percent.csv",row.names=F)
write.csv(AFDabs,"AFDabs_KromKosi_01percent.csv",row.names=F)
write.csv(NIELSEN,"NIELSEN_KromKosi_01percent.csv",row.names=F)
write.csv(VARLD,"VARLD_KromKosi_01percent.csv",row.names=F)
write.csv(FLK,"FLK_KromKosi_01percent.csv",row.names=F)
write.csv(DD,"DD_KromKosi_01percent.csv",row.names=F)




