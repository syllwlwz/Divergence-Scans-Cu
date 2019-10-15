######################################################################
#Analysis with potential promoter region before gene of 2 kb included#
######################################################################

#################################################
#load metrics and assign headers where necessary#
#################################################
require(gtools)

fst<-read.table("FstKromKosiHalleri_10SNPwin_new.csv",header=FALSE,sep="\t")
names(fst)<-c("Scaffold","Start_pos","End_pos","Fst")
fst<-subset(fst,!is.na(fst$Start_pos))
fst$Scaffold<-factor(fst$Scaffold,levels=mixedsort(levels(fst$Scaffold)))
fst<-fst[order(fst$Scaffold),]

DDdata<-read.table("DDKromKosiHalleri_10SNPwin_new.csv",sep="\t",header=F)
names(DDdata)<-c("Scaffold","Start_pos","End_pos","Mean_pos","Mean_allele_frequency_cohort1","Mean_allele_frequency_cohort2","Mean_abs_diff","Mean_pi_cohort1","Mean_pi_cohort2","Mean_raw_diff","DD")
DD<-data.frame(DDdata[,1:3],DDdata$DD)
names(DD)<-c("Scaffold","Start_pos","End_pos","DD")
DD$Scaffold<-factor(DD$Scaffold,levels=mixedsort(levels(DD$Scaffold)))
DD<-DD[order(DD$Scaffold),]

AFD<-data.frame(DDdata[,1:3],DDdata$Mean_raw_diff)
names(AFD)<-c("Scaffold","Start_pos","End_pos","AFD")
AFD$Scaffold<-factor(AFD$Scaffold,levels=mixedsort(levels(AFD$Scaffold)))
AFD<-AFD[order(AFD$Scaffold),]

AFDabs<-data.frame(DDdata[,1:3],DDdata$Mean_abs_diff)
names(AFDabs)<-c("Scaffold","Start_pos","End_pos","AFDabs")
AFDabs$Scaffold<-factor(AFDabs$Scaffold,levels=mixedsort(levels(AFDabs$Scaffold)))
AFDabs<-AFDabs[order(AFDabs$Scaffold),]

Nielsen<-read.table("NielsenKromKosi_folded_parallel2.csv",header=FALSE)
Nielsen<-Nielsen[-c(2,3)]
names(Nielsen)<-c("Scaffold","Start_pos","End_pos","Nielsen_diversity")
Nielsen$Scaffold<-factor(Nielsen$Scaffold,levels=mixedsort(levels(Nielsen$Scaffold)))
Nielsen<-Nielsen[order(Nielsen$Scaffold),]

Dxy<-read.table("DxyKromKosi_new.csv",header=FALSE,sep="\t")
Dxy<-subset(Dxy,!is.na(Dxy[,2]))
names(Dxy)<-c("Scaffold","Start_pos","End_pos","Dxy")
Dxy$Scaffold<-factor(Dxy$Scaffold,levels=mixedsort(levels(Dxy$Scaffold)))
Dxy<-Dxy[order(Dxy$Scaffold),]

Flk<-read.table("FlkKromKosi_min.csv",header=FALSE,sep="\t")
names(Flk)<-c("Scaffold","Start_pos","End_pos","pvalue","Flk")
Flk2<-Flk
library(gtools)
Flk2$Scaffold<-factor(Flk2$Scaffold,levels=mixedsort(levels(Flk2$Scaffold)))
Flk2<-Flk2[order(Flk2$Scaffold),]
plot(Flk2$Flk~fst$Fst)
Flk<-Flk2

VarLD<-read.table("KromKosi_VarLD_10SNP_windows",header=T,sep=",")
VarLD<-VarLD[-c(2,3)]
names(VarLD)<-c("Scaffold","Position","VarLD")
VarLD$Scaffold<-factor(VarLD$Scaffold,levels=mixedsort(levels(VarLD$Scaffold)))
VarLD<-VarLD[order(VarLD$Scaffold),]

TajimasD_Krom<-read.table("TajimasDKromHalleri_allSNPs.csv",header=FALSE,sep="\t")
names(TajimasD_Krom)<-c("Scaffold","Start_pos","End_pos","TajimasDKrom")
TajimasD_Krom$Scaffold<-factor(TajimasD_Krom$Scaffold,levels=mixedsort(levels(TajimasD_Krom$Scaffold)))
TajimasD_Krom<-TajimasD_Krom[order(TajimasD_Krom$Scaffold),]

TajimasD_Kosi<-read.table("TajimasDKosiHalleri_allSNPs.csv",header=FALSE,sep="\t")
names(TajimasD_Kosi)<-c("Scaffold","Start_pos","End_pos","TajimasDKosi")
TajimasD_Kosi$Scaffold<-factor(TajimasD_Kosi$Scaffold,levels=mixedsort(levels(TajimasD_Kosi$Scaffold)))
TajimasD_Kosi<-TajimasD_Kosi[order(TajimasD_Kosi$Scaffold),]

####################################
#combine all metrics into one table#
####################################

wholedata<-data.frame(fst,DD$DD,Nielsen$Nielsen_diversity,Dxy$Dxy,Flk$Flk,Flk$pvalue,VarLD$VarLD,AFD$AFD,AFDabs$AFD,DDdata$Mean_pi_cohort1,DDdata$Mean_pi_cohort2,DDdata$Mean_allele_frequency_cohort1,DDdata$Mean_allele_frequency_cohort2)
names(wholedata)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AFKrom","AFKosi")
#write.table(wholedata,"AllpostUGtestsprewikill.csv",col.names=TRUE,sep="\t",row.names=FALSE)

#################################
#kill windows bigger than 100 kb#
#################################
wholedata<-wholedata[wholedata$End_pos-wholedata$Start_pos<100000,]
write.table(wholedata,"AllpostUGtests_KromKosi.csv",col.names=TRUE,sep="\t",row.names=FALSE)
#2 windows excluded
#282,145 windows left

#############################################
#calculate average and median of window size#
#############################################
mean(wholedata$End_pos-wholedata$Start_pos)
#584.008
median(wholedata$End_pos-wholedata$Start_pos)
#295

#quantiles
quantile(wholedata$End_pos-wholedata$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#1%  10%  25%  75%  90%  99% 
#  47   99  163  579 1077 5269

#########################################
#find overlaps between windows and genes#
#########################################

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

## Load GFF files and add promotor region to genes

### Step 1: import genes from GFF files

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ahal_Gene")
#gff2<-gff
gff<-gff[gff$Type=="gene",]
gff$Start[gff$Strand=="+"]<-gff$Start[gff$Strand=="+"]-2000
gff$End[gff$Strand=="-"]<-gff$End[gff$Strand=="-"]+2000
gff_GRange<-GRanges(seqnames=tolower(gff$Scaffold),ranges=IRanges(start=gff$Start,end=gff$End),strand=gff$Strand)
AhalGenes_plusPromot <- gff_GRange
values(AhalGenes_plusPromot)<-gff$Ahal_Gene

###########################
#select 1% outlier windows#
###########################

Perc1Fst<-wholedata[wholedata$Fst>=quantile(wholedata$Fst,0.999)[[1]],]
Perc1DD<-wholedata[wholedata$DD<=quantile(wholedata$DD,0.001)[[1]],]
Perc1Nielsen<-wholedata[wholedata$Nielsen>=quantile(wholedata$Nielsen,0.999)[[1]],]
Perc1Dxy<-wholedata[wholedata$Dxy>=quantile(wholedata$Dxy,0.999)[[1]],]
Perc1Flk<-wholedata[(wholedata$Flk>=quantile(wholedata$Flk,0.999)[[1]])&wholedata$Flk_pvalue<0.05,]
Perc1VarLD<-wholedata[wholedata$VarLD>=quantile(wholedata$VarLD,0.999)[[1]],]
Perc1AFD<-wholedata[wholedata$AFD>=quantile(wholedata$AFD,0.999)[[1]],]
Perc1AFDabs<-wholedata[wholedata$AFDabs>=quantile(wholedata$AFDabs,0.999)[[1]],]

### Convert list of significant markers to a "Genomic Ranges" table. 
Perc1Fst_GRange<-GRanges(seqnames=tolower(Perc1Fst$Scaffold),ranges=IRanges(start=Perc1Fst$Start_pos,end=Perc1Fst$End_pos))
Perc1DD_GRange<-GRanges(seqnames=tolower(Perc1DD$Scaffold),ranges=IRanges(start=Perc1DD$Start_pos,end=Perc1DD$End_pos))
Perc1Nielsen_GRange<-GRanges(seqnames=tolower(Perc1Nielsen$Scaffold),ranges=IRanges(start=Perc1Nielsen$Start_pos,end=Perc1Nielsen$End_pos))
Perc1DD_GRange<-GRanges(seqnames=tolower(Perc1DD$Scaffold),ranges=IRanges(start=Perc1DD$Start_pos,end=Perc1DD$End_pos))
Perc1Dxy_GRange<-GRanges(seqnames=tolower(Perc1Dxy$Scaffold),ranges=IRanges(start=Perc1Dxy$Start_pos,end=Perc1Dxy$End_pos))
Perc1Flk_GRange<-GRanges(seqnames=tolower(Perc1Flk$Scaffold),ranges=IRanges(start=Perc1Flk$Start_pos,end=Perc1Flk$End_pos))
Perc1VarLD_GRange<-GRanges(seqnames=tolower(Perc1VarLD$Scaffold),ranges=IRanges(start=Perc1VarLD$Start_pos,end=Perc1VarLD$End_pos))
Perc1AFD_GRange<-GRanges(seqnames=tolower(Perc1AFD$Scaffold),ranges=IRanges(start=Perc1AFD$Start_pos,end=Perc1AFD$End_pos))
Perc1AFDabs_GRange<-GRanges(seqnames=tolower(Perc1AFDabs$Scaffold),ranges=IRanges(start=Perc1AFDabs$Start_pos,end=Perc1AFDabs$End_pos))
values(Perc1Fst_GRange)<-Perc1Fst[,4:16]
values(Perc1DD_GRange)<-Perc1DD[,4:16]
values(Perc1Nielsen_GRange)<-Perc1Nielsen[,4:16]
values(Perc1DD_GRange)<-Perc1DD[,4:16]
values(Perc1Dxy_GRange)<-Perc1Dxy[,4:16]
values(Perc1Flk_GRange)<-Perc1Flk[,4:16]
values(Perc1VarLD_GRange)<-Perc1VarLD[,4:16]
values(Perc1AFD_GRange)<-Perc1AFD[,4:16]
values(Perc1AFDabs_GRange)<-Perc1AFDabs[,4:16]

## Merge Significant Regions with Gene list in the halleri Genome

Perc1Fst_hallerigenes=mergeByOverlaps(Perc1Fst_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1Fst_hallerigenes_df=data.frame(as.data.frame(Perc1Fst_hallerigenes$Perc1Fst_GRange),as.data.frame(Perc1Fst_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1Fst_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Fst_hallerigenes_df$gene_start<-ifelse(Perc1Fst_hallerigenes_df$gene_strand=="+",Perc1Fst_hallerigenes_df$gene_start+2000,Perc1Fst_hallerigenes_df$gene_start)
Perc1Fst_hallerigenes_df$gene_end<-ifelse(Perc1Fst_hallerigenes_df$gene_strand=="-",Perc1Fst_hallerigenes_df$gene_end-2000,Perc1Fst_hallerigenes_df$gene_end)
Perc1Fst_hallerigenes_df$gene_size<-Perc1Fst_hallerigenes_df$gene_size-2000

Perc1DD_hallerigenes=mergeByOverlaps(Perc1DD_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1DD_hallerigenes_df=data.frame(as.data.frame(Perc1DD_hallerigenes$Perc1DD_GRange),as.data.frame(Perc1DD_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1DD_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1DD_hallerigenes_df$gene_start<-ifelse(Perc1DD_hallerigenes_df$gene_strand=="+",Perc1DD_hallerigenes_df$gene_start+2000,Perc1DD_hallerigenes_df$gene_start)
Perc1DD_hallerigenes_df$gene_end<-ifelse(Perc1DD_hallerigenes_df$gene_strand=="-",Perc1DD_hallerigenes_df$gene_end-2000,Perc1DD_hallerigenes_df$gene_end)
Perc1DD_hallerigenes_df$gene_size<-Perc1DD_hallerigenes_df$gene_size-2000

Perc1Nielsen_hallerigenes=mergeByOverlaps(Perc1Nielsen_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1Nielsen_hallerigenes_df=data.frame(as.data.frame(Perc1Nielsen_hallerigenes$Perc1Nielsen_GRange),as.data.frame(Perc1Nielsen_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1Nielsen_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Nielsen_hallerigenes_df$gene_start<-ifelse(Perc1Nielsen_hallerigenes_df$gene_strand=="+",Perc1Nielsen_hallerigenes_df$gene_start+2000,Perc1Nielsen_hallerigenes_df$gene_start)
Perc1Nielsen_hallerigenes_df$gene_end<-ifelse(Perc1Nielsen_hallerigenes_df$gene_strand=="-",Perc1Nielsen_hallerigenes_df$gene_end-2000,Perc1Nielsen_hallerigenes_df$gene_end)
Perc1Nielsen_hallerigenes_df$gene_size<-Perc1Nielsen_hallerigenes_df$gene_size-2000

Perc1Dxy_hallerigenes=mergeByOverlaps(Perc1Dxy_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1Dxy_hallerigenes_df=data.frame(as.data.frame(Perc1Dxy_hallerigenes$Perc1Dxy_GRange),as.data.frame(Perc1Dxy_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1Dxy_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi", "Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Dxy_hallerigenes_df$gene_start<-ifelse(Perc1Dxy_hallerigenes_df$gene_strand=="+",Perc1Dxy_hallerigenes_df$gene_start+2000,Perc1Dxy_hallerigenes_df$gene_start)
Perc1Dxy_hallerigenes_df$gene_end<-ifelse(Perc1Dxy_hallerigenes_df$gene_strand=="-",Perc1Dxy_hallerigenes_df$gene_end-2000,Perc1Dxy_hallerigenes_df$gene_end)
Perc1Dxy_hallerigenes_df$gene_size<-Perc1Dxy_hallerigenes_df$gene_size-2000

Perc1Flk_hallerigenes=mergeByOverlaps(Perc1Flk_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1Flk_hallerigenes_df=data.frame(as.data.frame(Perc1Flk_hallerigenes$Perc1Flk_GRange),as.data.frame(Perc1Flk_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1Flk_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1Flk_hallerigenes_df$gene_start<-ifelse(Perc1Flk_hallerigenes_df$gene_strand=="+",Perc1Flk_hallerigenes_df$gene_start+2000,Perc1Flk_hallerigenes_df$gene_start)
Perc1Flk_hallerigenes_df$gene_end<-ifelse(Perc1Flk_hallerigenes_df$gene_strand=="-",Perc1Flk_hallerigenes_df$gene_end-2000,Perc1Flk_hallerigenes_df$gene_end)
Perc1Flk_hallerigenes_df$gene_size<-Perc1Flk_hallerigenes_df$gene_size-2000

Perc1VarLD_hallerigenes=mergeByOverlaps(Perc1VarLD_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1VarLD_hallerigenes_df=data.frame(as.data.frame(Perc1VarLD_hallerigenes$Perc1VarLD_GRange),as.data.frame(Perc1VarLD_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1VarLD_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1VarLD_hallerigenes_df$gene_start<-ifelse(Perc1VarLD_hallerigenes_df$gene_strand=="+",Perc1VarLD_hallerigenes_df$gene_start+2000,Perc1VarLD_hallerigenes_df$gene_start)
Perc1VarLD_hallerigenes_df$gene_end<-ifelse(Perc1VarLD_hallerigenes_df$gene_strand=="-",Perc1VarLD_hallerigenes_df$gene_end-2000,Perc1VarLD_hallerigenes_df$gene_end)
Perc1VarLD_hallerigenes_df$gene_size<-Perc1VarLD_hallerigenes_df$gene_size-2000

Perc1AFD_hallerigenes=mergeByOverlaps(Perc1AFD_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1AFD_hallerigenes_df=data.frame(as.data.frame(Perc1AFD_hallerigenes$Perc1AFD_GRange),as.data.frame(Perc1AFD_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1AFD_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1AFD_hallerigenes_df$gene_start<-ifelse(Perc1AFD_hallerigenes_df$gene_strand=="+",Perc1AFD_hallerigenes_df$gene_start+2000,Perc1AFD_hallerigenes_df$gene_start)
Perc1AFD_hallerigenes_df$gene_end<-ifelse(Perc1AFD_hallerigenes_df$gene_strand=="-",Perc1AFD_hallerigenes_df$gene_end-2000,Perc1AFD_hallerigenes_df$gene_end)
Perc1AFD_hallerigenes_df$gene_size<-Perc1AFD_hallerigenes_df$gene_size-2000

Perc1AFDabs_hallerigenes=mergeByOverlaps(Perc1AFDabs_GRange,AhalGenes_plusPromot,type=c("any"))
Perc1AFDabs_hallerigenes_df=data.frame(as.data.frame(Perc1AFDabs_hallerigenes$Perc1AFDabs_GRange),as.data.frame(Perc1AFDabs_hallerigenes$AhalGenes_plusPromot))
colnames(Perc1AFDabs_hallerigenes_df)=c("Chr", "Window_start", "Window_end", "Window_width", "Window_strand", "Fst","DD","Nielsen","Dxy","Flk","Flk_pvalue","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AF_Krom","AF_Kosi","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene" )
Perc1AFDabs_hallerigenes_df$gene_start<-ifelse(Perc1AFDabs_hallerigenes_df$gene_strand=="+",Perc1AFDabs_hallerigenes_df$gene_start+2000,Perc1AFDabs_hallerigenes_df$gene_start)
Perc1AFDabs_hallerigenes_df$gene_end<-ifelse(Perc1AFDabs_hallerigenes_df$gene_strand=="-",Perc1AFDabs_hallerigenes_df$gene_end-2000,Perc1AFDabs_hallerigenes_df$gene_end)
Perc1AFDabs_hallerigenes_df$gene_size<-Perc1AFDabs_hallerigenes_df$gene_size-2000


######                                                                     #####
###### REGIONS NOT OVERLAPPING GENES WILL BE LOST FROM LIST AT THIS POINT  #####
######                                                                     #####

### Create a file with the list of Thaliana genes and descriptions

OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

#halleri genes to thaliana orthologues
Perc1Fst_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1Fst_hallerigenes_df,  by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1DD_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1DD_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1Nielsen_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1Nielsen_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1Dxy_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri, Perc1Dxy_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1Flk_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1Flk_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1VarLD_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1VarLD_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1AFD_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1AFD_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)
Perc1AFDabs_hallerigenes_df_w_orthogroup= merge(OG_Ahalleri,Perc1AFDabs_hallerigenes_df,by.y="Gene",by.x="Ah_ID", all.y=TRUE)

#add Mapman categories
Mapman<-read.table("Mapman_for_merge.table",sep="\t",header=T,fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
Mapman<-Mapman[Mapman$TYPE=="T",]
Mapman_df<-unique(Mapman[,c(3,9)])

Perc1Fst_hallerigenes_df_w_MM<-merge(Perc1Fst_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1DD_hallerigenes_df_w_MM<-merge(Perc1DD_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1Nielsen_hallerigenes_df_w_MM<-merge(Perc1Nielsen_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1Dxy_hallerigenes_df_w_MM<-merge(Perc1Dxy_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1Flk_hallerigenes_df_w_MM<-merge(Perc1Flk_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1VarLD_hallerigenes_df_w_MM<-merge(Perc1VarLD_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1AFD_hallerigenes_df_w_MM<-merge(Perc1AFD_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Perc1AFDabs_hallerigenes_df_w_MM<-merge(Perc1AFDabs_hallerigenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)

colnames(Perc1Fst_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1DD_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1Nielsen_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1Dxy_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1Flk_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1VarLD_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1AFD_hallerigenes_df_w_MM)[33]<-"Mapman_category"
colnames(Perc1AFDabs_hallerigenes_df_w_MM)[33]<-"Mapman_category"

#add Lists
All_lists<-read.table("All_lists_Cu_project.table",sep="\t",header=T,fill=T,quote="")

Perc1Fst_hallerigenes_df_w_MM_lists<-merge(Perc1Fst_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1DD_hallerigenes_df_w_MM_lists<-merge(Perc1DD_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1Nielsen_hallerigenes_df_w_MM_lists<-merge(Perc1Nielsen_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1Dxy_hallerigenes_df_w_MM_lists<-merge(Perc1Dxy_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1Flk_hallerigenes_df_w_MM_lists<-merge(Perc1Flk_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1VarLD_hallerigenes_df_w_MM_lists<-merge(Perc1VarLD_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1AFD_hallerigenes_df_w_MM_lists<-merge(Perc1AFD_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Perc1AFDabs_hallerigenes_df_w_MM_lists<-merge(Perc1AFDabs_hallerigenes_df_w_MM,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")

#change order

Perc1Fst_hallerigenes_df_w_MM_lists<-Perc1Fst_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1DD_hallerigenes_df_w_MM_lists<-Perc1DD_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1Nielsen_hallerigenes_df_w_MM_lists<-Perc1Nielsen_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1Dxy_hallerigenes_df_w_MM_lists<-Perc1Dxy_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1Flk_hallerigenes_df_w_MM_lists<-Perc1Flk_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1VarLD_hallerigenes_df_w_MM_lists<-Perc1VarLD_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1AFD_hallerigenes_df_w_MM_lists<-Perc1AFD_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]
Perc1AFDabs_hallerigenes_df_w_MM_lists<-Perc1AFDabs_hallerigenes_df_w_MM_lists[,c(2,28:32,10:14,1,3:9,15:76)]

#Generate file
write.table(Perc1Fst_hallerigenes_df_w_MM_lists,"Genes01percentFst_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1DD_hallerigenes_df_w_MM_lists,"Genes01percentDD_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1Nielsen_hallerigenes_df_w_MM_lists,"Genes01percentNielsen_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1Dxy_hallerigenes_df_w_MM_lists,"Genes01percentDxy_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1Flk_hallerigenes_df_w_MM_lists,"Genes01percentFlk_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1VarLD_hallerigenes_df_w_MM_lists,"Genes01percentVarLD_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1AFD_hallerigenes_df_w_MM_lists,"Genes01percentAFD_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)
write.table(Perc1AFDabs_hallerigenes_df_w_MM_lists,"Genes01percentAFDabs_hallerigenome_TDmean_KromKosi.txt", sep="\t", row.names=F,quote=F)

###############################
#Tajima's D and Fay and Wu's H#
###############################
str(TajimasD_Krom)
str(TajimasD_Kosi)
#241,858 windows Krom, 223,801 windows Kosi

#################################
#kill windows bigger than 100 kb#
#################################
TajimasD_Krom<-TajimasD_Krom[TajimasD_Krom$End_pos-TajimasD_Krom$Start_pos<100000,]
TajimasD_Kosi<-TajimasD_Kosi[TajimasD_Kosi$End_pos-TajimasD_Kosi$Start_pos<100000,]
write.table(TajimasD_Krom,"TD_Krom.csv",col.names=TRUE,sep="\t",row.names=FALSE)
write.table(TajimasD_Kosi,"TD_Kosi.csv",col.names=TRUE,sep="\t",row.names=FALSE)
#241,855 windows Krom, 223,797 windows Kosi)

#############################################
#calculate average and median of window size#
#############################################
mean(TajimasD_Krom$End_pos-TajimasD_Krom$Start_pos)
#680.8791
median(TajimasD_Krom$End_pos-TajimasD_Krom$Start_pos)
#341

#quantiles
quantile(TajimasD_Krom$End_pos-TajimasD_Krom$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#1%  10%  25%  75%  90%  99% 
#51.00  109.00  182.00  698.00 1318.00 5859.46

mean(TajimasD_Kosi$End_pos-TajimasD_Kosi$Start_pos)
#734.4091
median(TajimasD_Kosi$End_pos-TajimasD_Kosi$Start_pos)
#361

#quantiles
quantile(TajimasD_Kosi$End_pos-TajimasD_Kosi$Start_pos,probs=c(0.01,0.1,0.25,0.75,0.9,0.99))
#  1%  10%  25%  75%  90%  99% 
#51  113  190  750 1448 6393

options(java.parameters = "-Xmx12000m")

require(xlsx)
write.xlsx2(Perc1Fst_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="Fst_01%",col.names=TRUE,row.names=FALSE)
write.xlsx2(Perc1DD_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="DD_01%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Perc1Nielsen_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="Nielsen_01%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Perc1Dxy_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="Dxy_01%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Perc1Flk_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="Flk_01%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Perc1VarLD_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="VarLD_01%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Perc1AFD_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="AFD_01%",col.names=TRUE,row.names=FALSE,append=TRUE)
write.xlsx2(Perc1AFDabs_hallerigenes_df_w_MM_lists,"Genes_01percent_KromKosi.xlsx",sheetName="AFDabs_01%",col.names=TRUE,row.names=FALSE,append=TRUE)


#select values for % outliers#
FSTPerc<-rbind(names(quantile(wholedata$Fst,probs = seq(0.99,1,0.001))),quantile(wholedata$Fst,probs = seq(0.99,1,0.001)))
AFDabsPerc<-rbind(names(quantile(wholedata$AFDabs,probs = seq(0.99,1,0.001))),quantile(wholedata$AFDabs,probs = seq(0.99,1,0.001)))
NIELSENPerc<-rbind(names(quantile(wholedata$Nielsen,probs = seq(0.99,1,0.001))),quantile(wholedata$Nielsen,probs = seq(0.99,1,0.001)))
DXYPerc<-rbind(names(quantile(wholedata$Dxy,probs = seq(0.99,1,0.001))),quantile(wholedata$Dxy,probs = seq(0.99,1,0.001)))
VARLDPerc<-rbind(names(quantile(wholedata$VarLD,probs = seq(0.99,1,0.001))),quantile(wholedata$VarLD,probs = seq(0.99,1,0.001)))
FLKPerc<-rbind(names(quantile(wholedata$Flk,probs = seq(0.99,1,0.001))),quantile(wholedata$Flk,probs = seq(0.99,1,0.001)))
DDPerc<-rbind(names(quantile(wholedata$DD,probs = seq(0,0.01,0.001))),quantile(wholedata$DD,probs = seq(0,0.01,0.001)))

Quantiles<-rbind(FSTPerc,AFDabsPerc,NIELSENPerc,DXYPerc,VARLDPerc,FLKPerc,DDPerc)
Quantiles2<-cbind(c("FST","FST","AFDabs","AFDabs","NIELSEN","NIELSEN","DXY","DXY","VARLD","VARLD","FLK","FLK","DD","DD"),Quantiles)

write.xlsx2(Quantiles2,"Quantiles_KromKosi.xlsx",sheetName="Quantiles",col.names=FALSE,row.names=FALSE)



jpeg("Fst_Flk_KromKosi.jpeg", width=18, height=18, units="cm", res=1000)
plot(wholedata$Fst~wholedata$Flk)
abline(lm(wholedata$Fst~wholedata$Flk),col="red")
dev.off()



