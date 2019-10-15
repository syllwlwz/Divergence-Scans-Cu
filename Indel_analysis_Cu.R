data1<-read.table("KromINDELnew.table",header=TRUE,sep="\t")
data2<-read.table("KosiINDELnew.table",header=TRUE,sep="\t")
data<-data.frame(data1,data2)
str(data)
KK_Indels<-data[!((data$AC/data$AN)==(data$AC.1/data$AN.1)),]

SNPeff<-read.csv("Krom.Indels.ann.table",header=T,sep="\t")
KK_HC<-read.table("KromINDELnew3.vcf",sep="\t",header=F)

SNPeff2<-data.frame(SNPeff,KK_HC[,4:5])
names(SNPeff2)[4:5]<-c("REF","ALT")

SNPeffAnnsplit<-strsplit(as.character(SNPeff2$ANN),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",c(1:3,5))
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffcomp1<-cbind(SNPeff2[,c(1:2,4:5)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,na.exclude(SNPeffcomp1))
	}
require(GenomicRanges)
SNPeffdf2<-SNPeffdf
SNPeffdf<-SNPeffdf2
Gene<-gsub("^.*-","",SNPeffdf2$Gene)
SNPeffdf$Gene<-Gene

write.table(SNPeffdf,"SNPeff_KromKosi_Indels.table",sep="\t",row.names=F)
SNPeffdf<-read.table("SNPeff_KromKosi_Indels.table",sep="\t",header=T)

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ahal_Gene")
#gff2<-gff
gff<-gff[gff$Type=="gene",]
gff$Start<-ifelse(gff$Strand=="+",gff$Start+2000,gff$Start)
gff$End<-ifelse(gff$Strand=="-",gff$End-2000,gff$End)

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]

SNPeffdfMODIFIER_halgenes_df=merge(gff,SNPeffdfMODIFIER,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfMODIFIER_halgenes_GRange<-GRanges(seqnames=SNPeffdfMODIFIER_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER_halgenes_df$POS,end=SNPeffdfMODIFIER_halgenes_df$POS))
values(SNPeffdfMODIFIER_halgenes_GRange)<-SNPeffdfMODIFIER_halgenes_df[,c(1:10,13:17)]

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
SNPeffdfMODERATE_halgenes_df=merge(gff,SNPeffdfMODERATE,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfMODERATE_halgenes_GRange<-GRanges(seqnames=SNPeffdfMODERATE_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfMODERATE_halgenes_df$POS,end=SNPeffdfMODERATE_halgenes_df$POS))
values(SNPeffdfMODERATE_halgenes_GRange)<-SNPeffdfMODERATE_halgenes_df[,c(1:10,13:17)]

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH_halgenes_df=merge(gff,SNPeffdfHIGH,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfHIGH_halgenes_GRange<-GRanges(seqnames=SNPeffdfHIGH_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfHIGH_halgenes_df$POS,end=SNPeffdfHIGH_halgenes_df$POS))
values(SNPeffdfHIGH_halgenes_GRange)<-SNPeffdfHIGH_halgenes_df[,c(1:10,13:17)]

SNPeffdf_GRange<-GRanges(seqnames=SNPeffdf$CHROM,ranges=IRanges(start=SNPeffdf$POS,end=SNPeffdf$POS))
values(SNPeffdf_GRange)<-SNPeffdf[,3:8]

KK_Indels_GRange<-GRanges(seqnames=KK_Indels$CHROM,ranges=IRanges(start=KK_Indels$POS,end=KK_Indels$POS))
values(KK_Indels_GRange)<-KK_Indels[,c("AC","AC.1","AN","AN.1")]
SNP_KK_Indels=mergeByOverlaps(KK_Indels_GRange,SNPeffdf_GRange,type=c("any"))
SNP_KK_Indels_df=as.data.frame(SNP_KK_Indels)
SNP_KK_Indels_df2<-cbind(SNP_KK_Indels_df[,1:2],SNP_KK_Indels_df[,6:9],SNP_KK_Indels_df[,25:30])
names(SNP_KK_Indels_df2)<-c("Chrom","Pos","AC_Krom","AC_Kosi","AN_Krom","AN_Kosi","REF","ALT","Base","Characterization","Effect","Gene")

SNP_KK_Indels_diff<-SNP_KK_Indels_df2[abs((SNP_KK_Indels_df2$AC_Krom/SNP_KK_Indels_df2$AN_Krom)-(SNP_KK_Indels_df2$AC_Kosi/SNP_KK_Indels_df2$AN_Kosi))>0.4&(((SNP_KK_Indels_df2$AC_Krom/SNP_KK_Indels_df2$AN_Krom)<=0.1)|((SNP_KK_Indels_df2$AC_Krom/SNP_KK_Indels_df2$AN_Krom)>=0.9)),]

options(java.parameters = "-Xmx12000m")
require("xlsx")

#add A.thaliana orthologues
OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

#add Mapman categories
Mapman<-read.table("Mapman_for_merge.table",sep="\t",header=T,fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
Mapman<-Mapman[Mapman$TYPE=="T",]
Mapman_df<-unique(Mapman[,c(3,9)])

#add gene lists
All_lists<-read.table("All_lists_Cu_project.table",sep="\t",header=T,fill=T,quote="")

SNP_KK_Indels_diff_w_orthogroup= merge(OG_Ahalleri,SNP_KK_Indels_diff, by.y="Gene",by.x="Ah_ID", all.y=TRUE)
SNP_KK_Indels_diff_w_MM<-merge(SNP_KK_Indels_diff_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
SNP_KK_Indels_diff_w_MM_lists=merge(SNP_KK_Indels_diff_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
SNP_KK_Indels_diff_w_MM_lists<-SNP_KK_Indels_diff_w_MM_lists[,c(2,10:20,1,3:9,21:64)]
colnames(SNP_KK_Indels_diff_w_MM_lists)[21]<-"Mapman_category"

Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist<-SNP_KK_Indels_diff_w_MM_lists
Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Krom<-Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Krom/Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AN_Krom
Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Kosi<-Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AC_Kosi/Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist$AN_Kosi
names(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist)[4]<-"AF_Krom"
names(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist)[5]<-"AF_Kosi"
write.table(Thalgenes_described_SNP_KK_Indels_diff_w_orthogroupmetlist,"SNP_KK_Indels_diffabove04_hQD.txt", sep="\t", row.names=F,quote=F)


#############################################
#NossPais#
#############################################
data1<-read.table("NossINDELnew.table",header=TRUE,sep="\t")
data2<-read.table("PaisINDELnew.table",header=TRUE,sep="\t")
data<-data.frame(data1,data2)
str(data)
NP_Indels<-data[!((data$AC/data$AN)==(data$AC.1/data$AN.1)),]

SNPeff<-read.csv("Noss.Indels.ann.table",header=T,sep="\t")
NP_HC<-read.table("NossINDELnew3.vcf",sep="\t",header=F)

SNPeff2<-data.frame(SNPeff,NP_HC[,4:5])
names(SNPeff2)[4:5]<-c("REF","ALT")

SNPeffAnnsplit<-strsplit(as.character(SNPeff2$ANN),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",c(1:3,5))
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffcomp1<-cbind(SNPeff2[,c(1:2,4:5)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,na.exclude(SNPeffcomp1))
	}
require(GenomicRanges)
SNPeffdf2<-SNPeffdf
SNPeffdf<-SNPeffdf2
Gene<-gsub("^.*-","",SNPeffdf2$Gene)
SNPeffdf$Gene<-Gene

write.table(SNPeffdf,"SNPeff_NossPais_Indels.table",sep="\t",row.names=F)
SNPeffdf<-read.table("SNPeff_NossPais_Indels.table",sep="\t",header=T)

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ahal_Gene")
#gff2<-gff
gff<-gff[gff$Type=="gene",]
gff$Start<-ifelse(gff$Strand=="+",gff$Start+2000,gff$Start)
gff$End<-ifelse(gff$Strand=="-",gff$End-2000,gff$End)

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]

SNPeffdfMODIFIER_halgenes_df=merge(gff,SNPeffdfMODIFIER,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfMODIFIER_halgenes_GRange<-GRanges(seqnames=SNPeffdfMODIFIER_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER_halgenes_df$POS,end=SNPeffdfMODIFIER_halgenes_df$POS))
values(SNPeffdfMODIFIER_halgenes_GRange)<-SNPeffdfMODIFIER_halgenes_df[,c(1:10,13:17)]

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
SNPeffdfMODERATE_halgenes_df=merge(gff,SNPeffdfMODERATE,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfMODERATE_halgenes_GRange<-GRanges(seqnames=SNPeffdfMODERATE_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfMODERATE_halgenes_df$POS,end=SNPeffdfMODERATE_halgenes_df$POS))
values(SNPeffdfMODERATE_halgenes_GRange)<-SNPeffdfMODERATE_halgenes_df[,c(1:10,13:17)]

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH_halgenes_df=merge(gff,SNPeffdfHIGH,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfHIGH_halgenes_GRange<-GRanges(seqnames=SNPeffdfHIGH_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfHIGH_halgenes_df$POS,end=SNPeffdfHIGH_halgenes_df$POS))
values(SNPeffdfHIGH_halgenes_GRange)<-SNPeffdfHIGH_halgenes_df[,c(1:10,13:17)]

SNPeffdf_GRange<-GRanges(seqnames=SNPeffdf$CHROM,ranges=IRanges(start=SNPeffdf$POS,end=SNPeffdf$POS))
values(SNPeffdf_GRange)<-SNPeffdf[,3:8]

NP_Indels_GRange<-GRanges(seqnames=NP_Indels$CHROM,ranges=IRanges(start=NP_Indels$POS,end=NP_Indels$POS))
values(NP_Indels_GRange)<-NP_Indels[,c("AC","AC.1","AN","AN.1")]
SNP_NP_Indels=mergeByOverlaps(NP_Indels_GRange,SNPeffdf_GRange,type=c("any"))
SNP_NP_Indels_df=as.data.frame(SNP_NP_Indels)
SNP_NP_Indels_df2<-cbind(SNP_NP_Indels_df[,1:2],SNP_NP_Indels_df[,6:9],SNP_NP_Indels_df[,25:30])
names(SNP_NP_Indels_df2)<-c("Chrom","Pos","AC_Noss","AC_Pais","AN_Noss","AN_Pais","REF","ALT","Base","Characterization","Effect","Gene")

SNP_NP_Indels_diff<-SNP_NP_Indels_df2[abs((SNP_NP_Indels_df2$AC_Noss/SNP_NP_Indels_df2$AN_Noss)-(SNP_NP_Indels_df2$AC_Pais/SNP_NP_Indels_df2$AN_Pais))>0.4&(((SNP_NP_Indels_df2$AC_Noss/SNP_NP_Indels_df2$AN_Noss)<=0.1)|((SNP_NP_Indels_df2$AC_Noss/SNP_NP_Indels_df2$AN_Noss)>=0.9)),]

options(java.parameters = "-Xmx12000m")
require("xlsx")

#add A.thaliana orthologues
OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

#add Mapman categories
Mapman<-read.table("Mapman_for_merge.table",sep="\t",header=T,fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
Mapman<-Mapman[Mapman$TYPE=="T",]
Mapman_df<-unique(Mapman[,c(3,9)])

#add gene lists
All_lists<-read.table("All_lists_Cu_project.table",sep="\t",header=T,fill=T,quote="")

SNP_NP_Indels_diff_w_orthogroup= merge(OG_Ahalleri,SNP_NP_Indels_diff, by.y="Gene",by.x="Ah_ID", all.y=TRUE)
SNP_NP_Indels_diff_w_MM<-merge(SNP_NP_Indels_diff_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
SNP_NP_Indels_diff_w_MM_lists=merge(SNP_NP_Indels_diff_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
SNP_NP_Indels_diff_w_MM_lists<-SNP_NP_Indels_diff_w_MM_lists[,c(2,10:20,1,3:9,21:64)]
colnames(SNP_NP_Indels_diff_w_MM_lists)[21]<-"Mapman_category"

Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist<-SNP_NP_Indels_diff_w_MM_lists
Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist$AC_Noss<-Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist$AC_Noss/Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist$AN_Noss
Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist$AC_Pais<-Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist$AC_Pais/Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist$AN_Pais
names(Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist)[4]<-"AF_Noss"
names(Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist)[5]<-"AF_Pais"
write.table(Thalgenes_described_SNP_NP_Indels_diff_w_orthogroupmetlist,"SNP_NP_Indels_diffabove04_hQD.txt", sep="\t", row.names=F,quote=F)


#############################################
#LangBest#
#############################################
data1<-read.table("LangINDELnew.table",header=TRUE,sep="\t")
data2<-read.table("BestINDELnew.table",header=TRUE,sep="\t")
data<-data.frame(data1,data2)
str(data)
LB_Indels<-data[!((data$AC/data$AN)==(data$AC.1/data$AN.1)),]

SNPeff<-read.csv("Lang.Indels.ann.table",header=T,sep="\t")
LB_HC<-read.table("LangINDELnew3.vcf",sep="\t",header=F)

SNPeff2<-data.frame(SNPeff,LB_HC[,4:5])
names(SNPeff2)[4:5]<-c("REF","ALT")

SNPeffAnnsplit<-strsplit(as.character(SNPeff2$ANN),",",fixed=T)
SNPeffAnnsplit1<-lapply(SNPeffAnnsplit,"[",1:max(unlist(lapply(SNPeffAnnsplit,length))))
SNPeffAnnsplit2<-data.frame(matrix(unlist(SNPeffAnnsplit1),nrow=length(SNPeffAnnsplit1),byrow=T))
SNPeffdf<-data.frame()
for (i in 1:max(unlist(lapply(SNPeffAnnsplit,length))))
	{SNPeff1<-ifelse(!is.na(SNPeffAnnsplit2[,i]),strsplit(as.character(SNPeffAnnsplit2[,i]),"|",fixed=T),NA)
	SNPeffAnn1<-lapply(SNPeff1,"[",c(1:3,5))
	SNPeffAnndf1<-data.frame(matrix(unlist(SNPeffAnn1),nrow=length(SNPeffAnnsplit1),byrow=T))
	names(SNPeffAnndf1)<-c("Base","Characterization","Effect","Gene")
	SNPeffcomp1<-cbind(SNPeff2[,c(1:2,4:5)],SNPeffAnndf1)
	SNPeffdf<-rbind(SNPeffdf,na.exclude(SNPeffcomp1))
	}
require(GenomicRanges)
SNPeffdf2<-SNPeffdf
SNPeffdf<-SNPeffdf2
Gene<-gsub("^.*-","",SNPeffdf2$Gene)
SNPeffdf$Gene<-Gene

write.table(SNPeffdf,"SNPeff_LangBest_Indels.table",sep="\t",row.names=F)
SNPeffdf<-read.table("SNPeff_LangBest_Indels.table",sep="\t",header=T)

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ahal_Gene")
#gff2<-gff
gff<-gff[gff$Type=="gene",]
gff$Start<-ifelse(gff$Strand=="+",gff$Start+2000,gff$Start)
gff$End<-ifelse(gff$Strand=="-",gff$End-2000,gff$End)

SNPeffdfMODIFIER<-SNPeffdf[SNPeffdf$Effect=="MODIFIER",]

SNPeffdfMODIFIER_halgenes_df=merge(gff,SNPeffdfMODIFIER,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfMODIFIER_halgenes_GRange<-GRanges(seqnames=SNPeffdfMODIFIER_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfMODIFIER_halgenes_df$POS,end=SNPeffdfMODIFIER_halgenes_df$POS))
values(SNPeffdfMODIFIER_halgenes_GRange)<-SNPeffdfMODIFIER_halgenes_df[,c(1:10,13:17)]

SNPeffdfMODERATE<-SNPeffdf[SNPeffdf$Effect=="MODERATE",]
SNPeffdfMODERATE_halgenes_df=merge(gff,SNPeffdfMODERATE,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfMODERATE_halgenes_GRange<-GRanges(seqnames=SNPeffdfMODERATE_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfMODERATE_halgenes_df$POS,end=SNPeffdfMODERATE_halgenes_df$POS))
values(SNPeffdfMODERATE_halgenes_GRange)<-SNPeffdfMODERATE_halgenes_df[,c(1:10,13:17)]

SNPeffdfHIGH<-SNPeffdf[SNPeffdf$Effect=="HIGH",]
SNPeffdfHIGH_halgenes_df=merge(gff,SNPeffdfHIGH,by.y="Gene",by.x="Ahal_Gene")
SNPeffdfHIGH_halgenes_GRange<-GRanges(seqnames=SNPeffdfHIGH_halgenes_df$CHROM,ranges=IRanges(start=SNPeffdfHIGH_halgenes_df$POS,end=SNPeffdfHIGH_halgenes_df$POS))
values(SNPeffdfHIGH_halgenes_GRange)<-SNPeffdfHIGH_halgenes_df[,c(1:10,13:17)]

SNPeffdf_GRange<-GRanges(seqnames=SNPeffdf$CHROM,ranges=IRanges(start=SNPeffdf$POS,end=SNPeffdf$POS))
values(SNPeffdf_GRange)<-SNPeffdf[,3:8]

LB_Indels_GRange<-GRanges(seqnames=LB_Indels$CHROM,ranges=IRanges(start=LB_Indels$POS,end=LB_Indels$POS))
values(LB_Indels_GRange)<-LB_Indels[,c("AC","AC.1","AN","AN.1")]
SNP_LB_Indels=mergeByOverlaps(LB_Indels_GRange,SNPeffdf_GRange,type=c("any"))
SNP_LB_Indels_df=as.data.frame(SNP_LB_Indels)
SNP_LB_Indels_df2<-cbind(SNP_LB_Indels_df[,1:2],SNP_LB_Indels_df[,6:9],SNP_LB_Indels_df[,25:30])
names(SNP_LB_Indels_df2)<-c("Chrom","Pos","AC_Lang","AC_Best","AN_Lang","AN_Best","REF","ALT","Base","Characterization","Effect","Gene")

SNP_LB_Indels_diff<-SNP_LB_Indels_df2[abs((SNP_LB_Indels_df2$AC_Lang/SNP_LB_Indels_df2$AN_Lang)-(SNP_LB_Indels_df2$AC_Best/SNP_LB_Indels_df2$AN_Best))>0.4&(((SNP_LB_Indels_df2$AC_Lang/SNP_LB_Indels_df2$AN_Lang)<=0.1)|((SNP_LB_Indels_df2$AC_Lang/SNP_LB_Indels_df2$AN_Lang)>=0.9)),]

options(java.parameters = "-Xmx12000m")
require("xlsx")

#add A.thaliana orthologues
OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

#add Mapman categories
Mapman<-read.table("Mapman_for_merge.table",sep="\t",header=T,fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
Mapman<-Mapman[Mapman$TYPE=="T",]
Mapman_df<-unique(Mapman[,c(3,9)])

#add gene lists
All_lists<-read.table("All_lists_Cu_project.table",sep="\t",header=T,fill=T,quote="")

SNP_LB_Indels_diff_w_orthogroup= merge(OG_Ahalleri,SNP_LB_Indels_diff, by.y="Gene",by.x="Ah_ID", all.y=TRUE)
SNP_LB_Indels_diff_w_MM<-merge(SNP_LB_Indels_diff_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
SNP_LB_Indels_diff_w_MM_lists=merge(SNP_LB_Indels_diff_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
SNP_LB_Indels_diff_w_MM_lists<-SNP_LB_Indels_diff_w_MM_lists[,c(2,10:20,1,3:9,21:64)]
colnames(SNP_LB_Indels_diff_w_MM_lists)[21]<-"Mapman_category"

Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist<-SNP_LB_Indels_diff_w_MM_lists
Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist$AC_Lang<-Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist$AC_Lang/Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist$AN_Lang
Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist$AC_Best<-Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist$AC_Best/Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist$AN_Best
names(Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist)[4]<-"AF_Lang"
names(Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist)[5]<-"AF_Best"
write.table(Thalgenes_described_SNP_LB_Indels_diff_w_orthogroupmetlist,"SNP_LB_Indels_diffabove04_hQD.txt", sep="\t", row.names=F,quote=F)


