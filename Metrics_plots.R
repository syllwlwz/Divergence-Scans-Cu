

###############
#Metrics plots#
###############

#KK
#change dir
setwd("F:/Cu_project/Genome_scans/Genom_scans_Cu/Krom-Kosi")

options(java.parameters = "-Xmx12g")
require(xlsx)

RNASeqR<-read.xlsx2("F://Cu_project//RNA-Seq//Sample_counts//Deseq2_medofratios_stricter_Mapmanadded.xlsx",6,header=T)
RNASeqS<-read.xlsx2("F://Cu_project//RNA-Seq//Sample_counts//Deseq2_medofratios_stricter_Mapmanadded.xlsx",7,header=T)

Candidates2<-read.xlsx2("Genes_KromKosi_wVarLD.xlsx",1,header=TRUE)

Candidates5<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates3<-Candidates5[!is.na(Candidates5$Gene),]
Candidates4<-Candidates3[!Candidates3$Gene=="",]
Candidates<-Candidates4[Candidates4$Ah_ID%in%RNASeqR$Ah_ID|Candidates4$Ah_ID%in%RNASeqS$Ah_ID,]
#31
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0004954"|Candidates4$Ah_ID=="AHAL_G0003865"|Candidates4$Ah_ID=="AHAL_G0007270"|Candidates4$Ah_ID=="AHAL_G0006303"|Candidates4$Ah_ID=="AHAL_G0006286"|Candidates4$Ah_ID=="AHAL_G0010056"|Candidates4$Ah_ID=="AHAL_G0025694"|Candidates4$Ah_ID=="AHAL_G0030993"|Candidates4$Ah_ID=="AHAL_G0007270"|Candidates4$Ah_ID=="AHAL_G0025599"|Candidates4$Ah_ID=="AHAL_G0012884"|Candidates4$Ah_ID=="AHAL_G0001510"|Candidates4$Ah_ID=="AHAL_G0014338"|Candidates4$Ah_ID=="AHAL_G0014700"|Candidates4$Ah_ID=="AHAL_G0015690"|Candidates4$Ah_ID=="AHAL_G0002108"|Candidates4$Ah_ID=="AHAL_G0006748"|Candidates4$Ah_ID=="AHAL_G0007128"|Candidates4$Ah_ID=="AHAL_G0017516"
|Candidates4$Ah_ID=="AHAL_G0024590"|Candidates4$Ah_ID=="AHAL_G0008885"|Candidates4$Ah_ID=="AHAL_G0024044"|Candidates4$Ah_ID=="AHAL_G0015076"|Candidates4$Ah_ID=="AHAL_G0012967"|Candidates4$Ah_ID=="AHAL_G0012903"|Candidates4$Ah_ID=="AHAL_G0012433"|Candidates4$Ah_ID=="AHAL_G0012418"|Candidates4$Ah_ID=="AHAL_G0011147"|Candidates4$Ah_ID=="AHAL_G0008318"|Candidates4$Ah_ID=="AHAL_G0017626"|Candidates4$Ah_ID=="AHAL_G0017746"|Candidates4$Ah_ID=="AHAL_G0024590"|Candidates4$Ah_ID=="AHAL_G0028807"|Candidates4$Ah_ID=="AHAL_G0021363"|Candidates4$Ah_ID=="AHAL_G0030810"|Candidates4$Ah_ID=="AHAL_G0030485"|Candidates4$Ah_ID=="AHAL_G0030123"|Candidates4$Ah_ID=="AHAL_G0029902"|Candidates4$Ah_ID=="AHAL_G0006377",]
test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

#convergent candidates
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0012418"|Candidates4$Ah_ID=="AHAL_G0015692"|Candidates4$Ah_ID=="AHAL_G0021044"|Candidates4$Ah_ID=="AHAL_G0025603"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0009828"|Candidates4$Ah_ID=="AHAL_G0027258"|Candidates4$Ah_ID=="AHAL_G0010679"|Candidates4$Ah_ID=="AHAL_G0013166"|Candidates4$Ah_ID=="AHAL_G0024026",]
test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

Candidates2<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//Genes_NossPais_wVarLD.xlsx",1,header=TRUE)

Candidates5<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates3<-Candidates5[!is.na(Candidates5$Gene),]
Candidates4<-Candidates3[!Candidates3$Gene=="",]

Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0021099"|Candidates4$Ah_ID=="AHAL_G0021010"|Candidates4$Ah_ID=="AHAL_G0023747"|Candidates4$Ah_ID=="AHAL_G0002633"|Candidates4$Ah_ID=="AHAL_G0031553"|Candidates4$Ah_ID=="AHAL_G0031551"|Candidates4$Ah_ID=="AHAL_GAHAL_G0001510"|Candidates4$Ah_ID=="AHAL_G0020094"|Candidates4$Ah_ID=="AHAL_G0023512"|Candidates4$Ah_ID=="AHAL_G0031217"|Candidates4$Ah_ID=="AHAL_G0015312"|Candidates4$Ah_ID=="AHAL_G0009746"|Candidates4$Ah_ID=="AHAL_G0016514"|Candidates4$Ah_ID=="AHAL_G0030227"|Candidates4$Ah_ID=="AHAL_G0021925"|Candidates4$Ah_ID=="AHAL_G0024384"|Candidates4$Ah_ID=="AHAL_G0010752"|Candidates4$Ah_ID=="AHAL_G0011186"|Candidates4$Ah_ID=="AHAL_G0011222"|Candidates4$Ah_ID=="AHAL_G0015887"|Candidates4$Ah_ID=="AHAL_G0016505"|Candidates4$Ah_ID=="AHAL_G0009786"|Candidates4$Ah_ID=="AHAL_G0009226"|Candidates4$Ah_ID=="AHAL_G0000831"|Candidates4$Ah_ID=="AHAL_G0001230"|Candidates4$Ah_ID=="AHAL_G0006042",]
test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

#convergent
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0012418"|Candidates4$Ah_ID=="AHAL_G0015692"|Candidates4$Ah_ID=="AHAL_G0021044"|Candidates4$Ah_ID=="AHAL_G0025603"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0009828"|Candidates4$Ah_ID=="AHAL_G0027258"|Candidates4$Ah_ID=="AHAL_G0010679"|Candidates4$Ah_ID=="AHAL_G0013166"|Candidates4$Ah_ID=="AHAL_G0024026",]
test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)



Candidates2<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Lang-Best//Genes_LangBest_wVarLD.xlsx",1,header=TRUE)

Candidates5<-Candidates2[!duplicated(Candidates2$Gene),]
Candidates3<-Candidates5[!is.na(Candidates5$Gene),]
Candidates4<-Candidates3[!Candidates3$Gene=="",]

Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0001278"|Candidates4$Ah_ID=="AHAL_G0009430"|Candidates4$Ah_ID=="AHAL_G0009431"|Candidates4$Ah_ID=="AHAL_G0006543"|Candidates4$Ah_ID=="AHAL_G0019263"|Candidates4$Ah_ID=="AHAL_G0019770"
|Candidates4$Ah_ID=="AHAL_G0020401"|Candidates4$Ah_ID=="AHAL_G0025606"|Candidates4$Ah_ID=="AHAL_G0001072"|Candidates4$Ah_ID=="AHAL_G0015970"|Candidates4$Ah_ID=="AHAL_G0019595"|Candidates4$Ah_ID=="AHAL_G0028240"|Candidates4$Ah_ID=="AHAL_G0022603"|Candidates4$Ah_ID=="AHAL_G0031052"|Candidates4$Ah_ID=="AHAL_G0018180"|Candidates4$Ah_ID=="AHAL_G0004980"|Candidates4$Ah_ID=="AHAL_G0003038"|Candidates4$Ah_ID=="AHAL_G0018167"|Candidates4$Ah_ID=="AHAL_G0008064"|Candidates4$Ah_ID=="AHAL_G0009794"|Candidates4$Ah_ID=="AHAL_G0016498"|Candidates4$Ah_ID=="AHAL_G0024495"|Candidates4$Ah_ID=="AHAL_G0028079"|Candidates4$Ah_ID=="AHAL_G0028319"|Candidates4$Ah_ID==""|Candidates4$Ah_ID==""|Candidates4$Ah_ID=="AHAL_G0028320"|Candidates4$Ah_ID=="AHAL_G0023419"|Candidates4$Ah_ID=="AHAL_G0025407"|Candidates4$Ah_ID=="AHAL_G0032108"|Candidates4$Ah_ID=="AHAL_G0032173"|Candidates4$Ah_ID=="AHAL_G0006377",]
test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


#convergent
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0012418"|Candidates4$Ah_ID=="AHAL_G0015692"|Candidates4$Ah_ID=="AHAL_G0021044"|Candidates4$Ah_ID=="AHAL_G0025603"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0009828"|Candidates4$Ah_ID=="AHAL_G0027258"|Candidates4$Ah_ID=="AHAL_G0010679"|Candidates4$Ah_ID=="AHAL_G0013166"|Candidates4$Ah_ID=="AHAL_G0024026",]
test2<-as.character(Candidates$Gene)
Candidates$Gene<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)




AFdataKK<-read.table("AllpostUGtests_KromKosi.csv",sep="\t",header=TRUE,encoding="UTF-8")
names(AFdataKK)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_p_value","VarLD","AFD","AFDabs","Pi_Krom","Pi_Kosi","AFKrom","AFKosi")
test2KK<-as.character(AFdataKK$Scaffold)
AFdataKK$Scaffold<-test2KK

genelist<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(genelist)<-c("Scaffold","Source","Type","Start_pos","End_pos","Score","Strand","Phase","Attributes","Ah_ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)
genelist<-genelist[!duplicated(genelist[,-9]),]

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed
test3<-as.character(Candidates$Chr)
Candidates$Chr<-test3

winmiddleKK<-c((AFdataKK$Start_pos+AFdataKK$End_pos)/2)
AFdataKK<-cbind(AFdataKK,winmiddleKK)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestKK", j, sep ="")
 	AF1KK<-AFdataKK[AFdataKK$Scaffold==Candidates$Chr[j],]
	AFofinterestKK<-AF1KK[AF1KK$winmiddleKK>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1KK$winmiddleKK<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterestKK)
	}
rm(AFofinterestKK)
list<-ls(pattern="^AFofinterestKK")

TD_Krom<-read.csv("TD_Krom.csv",header=TRUE,sep="\t")
TD_Kosi<-read.csv("TD_Kosi.csv",header=TRUE,sep="\t")
require(vegan)
TD_Krom$TajimasDKrom<-decostand(TD_Krom$TajimasDKrom,"standardize")
TD_Kosi$TajimasDKosi<-decostand(TD_Kosi$TajimasDKosi,"standardize")

test2Krom<-as.character(TD_Krom$Scaffold)
TD_Krom$Scaffold<-test2Krom
winmiddleTDKrom<-c((TD_Krom$Start_pos+TD_Krom$End_pos)/2)
TD_Krom<-cbind(TD_Krom,winmiddleTDKrom)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKromofinterest", j, sep ="")
 	AF1Krom<-TD_Krom[TD_Krom$Scaffold==Candidates$Chr[j],]
	TDKromofinterest<-AF1Krom[AF1Krom$winmiddleTDKrom>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Krom$winmiddleTDKrom<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKromofinterest)
	}
rm(TDKromofinterest)
list<-ls(pattern="^TDKromofinterest")

test2Kosi<-as.character(TD_Kosi$Scaffold)
TD_Kosi$Scaffold<-test2Kosi
winmiddleTDKosi<-c((TD_Kosi$Start_pos+TD_Kosi$End_pos)/2)
TD_Kosi<-cbind(TD_Kosi,winmiddleTDKosi)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDKosiofinterest", j, sep ="")
 	AF1Kosi<-TD_Kosi[TD_Kosi$Scaffold==Candidates$Chr[j],]
	TDKosiofinterest<-AF1Kosi[AF1Kosi$winmiddleTDKosi>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Kosi$winmiddleTDKosi<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDKosiofinterest)
	}
rm(TDKosiofinterest)
list<-ls(pattern="^TDKosiofinterest")

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)

#Quantiles:
DD_abs_percKK<-quantile(AFdataKK$DD,0.001)
Fst_percKK<-quantile(AFdataKK$Fst,0.999)
Nielsen_percKK<-quantile(AFdataKK$Nielsen,0.999)
Dxy_percKK<-quantile(AFdataKK$Dxy,0.999)
Flk_percKK<-quantile(AFdataKK$Flk,0.999)
VarLD_percKK<-quantile(AFdataKK$VarLD,0.999)

for (j in 1:nrow(Candidates))
	{nam <- paste("Genesofinterest", j, sep ="")
 	Genes1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Genesofinterest<-Genes1[Genes1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Genes1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Genes1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Genesofinterest)
	}

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

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

Sweed_Kosi<-dataSweed_Kosi

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

Sweed_Krom<-dataSweed_Krom

test2KK<-as.character(tolower(Sweed_Krom$Scaffold))
Sweed_Krom$Scaffold<-test2KK
for (j in 1:nrow(Candidates))
	{nam <- paste("SweedofinterestKK", j, sep ="")
 	AF1KK<-Sweed_Krom[Sweed_Krom$Scaffold==Candidates$Chr[j],]
	interval<-which(AF1KK$Position>=(Candidates$gene_start[j]-windowsize[j])&AF1KK$Position<=(Candidates$gene_end[j]+windowsize[j]))
	interval2<-c(min(interval)-1,interval,max(interval)+1)
	SweedofinterestKK<-AF1KK[interval2,]
	assign(nam, SweedofinterestKK)
	}
rm(SweedofinterestKK)
list<-ls(pattern="^SweedofinterestKK")

test2Kosi<-as.character(tolower(Sweed_Kosi$Scaffold))
Sweed_Kosi$Scaffold<-test2Kosi
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Kosiofinterest", j, sep ="")
 	AF1Kosi<-Sweed_Kosi[Sweed_Kosi$Scaffold==Candidates$Chr[j],]
	Sweed_Kosiofinterest<-AF1Kosi[AF1Kosi$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Kosi$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Kosiofinterest)
	}
rm(Sweed_Kosiofinterest)
list<-ls(pattern="^Sweed_Kosiofinterest")


genelist2<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(genelist2)<-c("Scaffold","Source","Type","Start_pos","End_pos","Score","Strand","Phase","Attributes","Ah_ID")
exonlist<-na.exclude(genelist2[genelist2$Type=="exon",])
exonlist<-droplevels(exonlist)
exonlist$Scaffold<-as.character(exonlist$Scaffold)
colnames(exonlist)[10]<-"ID"
exonlist<-exonlist[!duplicated(exonlist[,-9]),]

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterest", j, sep ="")
 	Exons1<-na.exclude(exonlist[exonlist$Scaffold==Candidates$Chr[j],])
	Exonsofinterest<-Exons1[Exons1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Exons1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Exons1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Exonsofinterest)
	}
rm(Exonsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

Candidates$Ah_ID<-as.character(Candidates$Ah_ID)

for (j in 1:nrow(Candidates))
	{nam <- paste("Exonsofinterestcand", j, sep ="")
 	Exons1<-exonlist[exonlist$Scaffold==Candidates$Chr[j],]
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Ah_ID[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)

#NossPais

AFdataNP<-read.table("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
names(AFdataNP)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_p_value","VarLD","AFD","AFDabs","Pi_Noss","Pi_Pais","AF_Noss","AF_Pais")
test2NP<-as.character(AFdataNP$Scaffold)
AFdataNP$Scaffold<-test2NP

winmiddleNP<-c((AFdataNP$Start_pos+AFdataNP$End_pos)/2)
AFdataNP<-cbind(AFdataNP,winmiddleNP)
names(AFdataNP)[17]<-"winmiddleNP"
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestNP", j, sep ="")
 	AF1NP<-AFdataNP[AFdataNP$Scaffold==Candidates$Chr[j],]
	AFofinterestNP<-AF1NP[AF1NP$winmiddleNP>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1NP$winmiddleNP<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterestNP)
	}
rm(AFofinterestNP)
list<-ls(pattern="^AFofinterestNP")

TD_Noss<-read.csv("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//TD_Noss.csv",header=TRUE,sep="\t")
TD_Pais<-read.csv("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//TD_Pais.csv",header=TRUE,sep="\t")
require(vegan)
TD_Noss$TajimasDNoss<-decostand(TD_Noss$TajimasDNoss,"standardize")
TD_Pais$TajimasDPais<-decostand(TD_Pais$TajimasDPais,"standardize")

test2Noss<-as.character(TD_Noss$Scaffold)
TD_Noss$Scaffold<-test2Noss
winmiddleTDNoss<-c((TD_Noss$Start_pos+TD_Noss$End_pos)/2)
TD_Noss<-cbind(TD_Noss,winmiddleTDNoss)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDNossofinterest", j, sep ="")
 	AF1Noss<-TD_Noss[TD_Noss$Scaffold==Candidates$Chr[j],]
	TDNossofinterest<-AF1Noss[AF1Noss$winmiddleTDNoss>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Noss$winmiddleTDNoss<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDNossofinterest)
	}
rm(TDNossofinterest)
list<-ls(pattern="^TDNossofinterest")

test2Pais<-as.character(TD_Pais$Scaffold)
TD_Pais$Scaffold<-test2Pais
winmiddleTDPais<-c((TD_Pais$Start_pos+TD_Pais$End_pos)/2)
TD_Pais<-cbind(TD_Pais,winmiddleTDPais)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDPaisofinterest", j, sep ="")
 	AF1Pais<-TD_Pais[TD_Pais$Scaffold==Candidates$Chr[j],]
	TDPaisofinterest<-AF1Pais[AF1Pais$winmiddleTDPais>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Pais$winmiddleTDPais<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDPaisofinterest)
	}
rm(TDPaisofinterest)
list<-ls(pattern="^TDPaisofinterest")

#Quantiles:
DD_abs_percNP<-quantile(AFdataNP$DD,0.001)
Fst_percNP<-quantile(AFdataNP$Fst,0.999)
Nielsen_percNP<-quantile(AFdataNP$Nielsen,0.999)
Dxy_percNP<-quantile(AFdataNP$Dxy,0.999)
Flk_percNP<-quantile(AFdataNP$Flk,0.999)
VarLD_percNP<-quantile(AFdataNP$VarLD,0.999)

#change dir
setwd("F:/Cu_project/Genome_scans/Genom_scans_Cu/Noss-Pais")
filelist = list.files(pattern = "SweeD_Report.Pais_SweeD_allSNPs*") 
require(gtools)
filelist <- mixedsort(filelist)
scaff1<-read.table("SweeD_Report.Pais_SweeD_allSNPs.chr1",sep="\t",header=T)
datafr = do.call(rbind,lapply(filelist[-1],function(fn)read.table(fn,header=F,sep="\t",colClasses =c("numeric","numeric","numeric") ))) #merge all count files
chr<-strsplit(filelist,"\\.")
chr1<-lapply(chr, tail, n = 1L)
chr2<-data.frame(matrix(unlist(chr1),nrow=length(chr1),byrow=T))[,1]
#datanrow = do.call(rbind,lapply(filelist,function(fn)nrow(read.table(fn,header=FALSE, sep="\t"))))
#all 37500
Scaffold<-rep(as.character(chr2),each=37500)
colnames(datafr)<-colnames(scaff1)
datafr2<-rbind(scaff1,datafr)
dataSweed_Pais<-data.frame(Scaffold,datafr2)
names(dataSweed_Pais)<-c("Scaffold","Position","Likelihood","Alpha")

Sweed_Pais<-dataSweed_Pais

filelist = list.files(pattern = "SweeD_Report.Noss_SweeD_allSNPs*") 
require(gtools)
filelist <- mixedsort(filelist)
scaff1<-read.table("SweeD_Report.Noss_SweeD_allSNPs.chr1",sep="\t",header=T)
datafr = do.call(rbind,lapply(filelist[-1],function(fn)read.table(fn,header=F,sep="\t",colClasses =c("numeric","numeric","numeric") ))) #merge all count files
chr<-strsplit(filelist,"\\.")
chr1<-lapply(chr, tail, n = 1L)
chr2<-data.frame(matrix(unlist(chr1),nrow=length(chr1),byrow=T))[,1]
#datanrow = do.call(rbind,lapply(filelist,function(fn)nrow(read.table(fn,header=FALSE, sep="\t"))))
#all 37500
Scaffold<-rep(as.character(chr2),each=37500)
colnames(datafr)<-colnames(scaff1)
datafr2<-rbind(scaff1,datafr)
dataSweed_Noss<-data.frame(Scaffold,datafr2)
names(dataSweed_Noss)<-c("Scaffold","Position","Likelihood","Alpha")

SweedNP<-dataSweed_Noss


test2NP<-as.character(tolower(SweedNP$Scaffold))
SweedNP$Scaffold<-test2NP
for (j in 1:nrow(Candidates))
	{nam <- paste("SweedofinterestNP", j, sep ="")
 	AF1NP<-SweedNP[SweedNP$Scaffold==Candidates$Chr[j],]
	interval<-which(AF1NP$Position>=(Candidates$gene_start[j]-windowsize[j])&AF1NP$Position<=(Candidates$gene_end[j]+windowsize[j]))
	interval2<-c(min(interval)-1,interval,max(interval)+1)
	SweedofinterestNP<-AF1NP[interval2,]
	assign(nam, SweedofinterestNP)
	}
rm(SweedofinterestNP)
list<-ls(pattern="^SweedofinterestNP")

test2Pais<-as.character(tolower(Sweed_Pais$Scaffold))
Sweed_Pais$Scaffold<-test2Pais
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Paisofinterest", j, sep ="")
 	AF1Pais<-Sweed_Pais[Sweed_Pais$Scaffold==Candidates$Chr[j],]
	Sweed_Paisofinterest<-AF1Pais[AF1Pais$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Pais$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Paisofinterest)
	}
rm(Sweed_Paisofinterest)
list<-ls(pattern="^Sweed_Paisofinterest")


#LangBest
#change dir
setwd("F:/Cu_project/Genome_scans/Genom_scans_Cu/Lang-Best")

AFdataLB<-read.table("AllpostUGtests.csv",sep="\t",header=TRUE,encoding="UTF-8")
names(AFdataLB)<-c("Scaffold","Start_pos","End_pos","Fst","DD","Nielsen","Dxy","Flk","Flk_p_value","VarLD","AFD","AFDabs","Pi_Lang","Pi_Best","AF_Lang","AF_Best")
test2LB<-as.character(AFdataLB$Scaffold)
AFdataLB$Scaffold<-test2LB

winmiddleLB<-c((AFdataLB$Start_pos+AFdataLB$End_pos)/2)
AFdataLB<-cbind(AFdataLB,winmiddleLB)
for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterestLB", j, sep ="")
 	AF1LB<-AFdataLB[AFdataLB$Scaffold==Candidates$Chr[j],]
	AFofinterestLB<-AF1LB[AF1LB$winmiddleLB>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1LB$winmiddleLB<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, AFofinterestLB)
	}
rm(AFofinterestLB)
list<-ls(pattern="^AFofinterestLB")

TD_Lang<-read.csv("TD_Lang.csv",header=TRUE,sep="\t")
TD_Best<-read.csv("TD_Best.csv",header=TRUE,sep="\t")
require(vegan)
TD_Lang$TajimasDLang<-decostand(TD_Lang$TajimasDLang,"standardize")
TD_Best$TajimasDBest<-decostand(TD_Best$TajimasDBest,"standardize")

test2Lang<-as.character(TD_Lang$Scaffold)
TD_Lang$Scaffold<-test2Lang
winmiddleTDLang<-c((TD_Lang$Start_pos+TD_Lang$End_pos)/2)
TD_Lang<-cbind(TD_Lang,winmiddleTDLang)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDLangofinterest", j, sep ="")
 	AF1Lang<-TD_Lang[TD_Lang$Scaffold==Candidates$Chr[j],]
	TDLangofinterest<-AF1Lang[AF1Lang$winmiddleTDLang>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Lang$winmiddleTDLang<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDLangofinterest)
	}
rm(TDLangofinterest)
list<-ls(pattern="^TDLangofinterest")

test2Best<-as.character(TD_Best$Scaffold)
TD_Best$Scaffold<-test2Best
winmiddleTDBest<-c((TD_Best$Start_pos+TD_Best$End_pos)/2)
TD_Best<-cbind(TD_Best,winmiddleTDBest)
for (j in 1:nrow(Candidates))
	{nam <- paste("TDBestofinterest", j, sep ="")
 	AF1Best<-TD_Best[TD_Best$Scaffold==Candidates$Chr[j],]
	TDBestofinterest<-AF1Best[AF1Best$winmiddleTDBest>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Best$winmiddleTDBest<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, TDBestofinterest)
	}
rm(TDBestofinterest)
list<-ls(pattern="^TDBestofinterest")

#Quantiles:
DD_abs_percLB<-quantile(AFdataLB$DD,0.001)
Fst_percLB<-quantile(AFdataLB$Fst,0.999)
Nielsen_percLB<-quantile(AFdataLB$Nielsen,0.999)
Dxy_percLB<-quantile(AFdataLB$Dxy,0.999)
Flk_percLB<-quantile(AFdataLB$Flk,0.999)
VarLD_percLB<-quantile(AFdataLB$VarLD,0.999)

filelist = list.files(pattern = "SweeD_Report.Best_SweeD_allSNPs*") 
require(gtools)
filelist <- mixedsort(filelist)
scaff1<-read.table("SweeD_Report.Best_SweeD_allSNPs.chr1",sep="\t",header=T)
datafr = do.call(rbind,lapply(filelist[-1],function(fn)read.table(fn,header=F,sep="\t",colClasses =c("numeric","numeric","numeric") ))) #merge all count files
chr<-strsplit(filelist,"\\.")
chr1<-lapply(chr, tail, n = 1L)
chr2<-data.frame(matrix(unlist(chr1),nrow=length(chr1),byrow=T))[,1]
#datanrow = do.call(rbind,lapply(filelist,function(fn)nrow(read.table(fn,header=FALSE, sep="\t"))))
#all 37500
Scaffold<-rep(as.character(chr2),each=37500)
colnames(datafr)<-colnames(scaff1)
datafr2<-rbind(scaff1,datafr)
dataSweed_Best<-data.frame(Scaffold,datafr2)
names(dataSweed_Best)<-c("Scaffold","Position","Likelihood","Alpha")

Sweed_Best<-dataSweed_Best

filelist = list.files(pattern = "SweeD_Report.Lang_SweeD_allSNPs*") 
require(gtools)
filelist <- mixedsort(filelist)
scaff1<-read.table("SweeD_Report.Lang_SweeD_allSNPs.chr1",sep="\t",header=T)
datafr = do.call(rbind,lapply(filelist[-1],function(fn)read.table(fn,header=F,sep="\t",colClasses =c("numeric","numeric","numeric") ))) #merge all count files
chr<-strsplit(filelist,"\\.")
chr1<-lapply(chr, tail, n = 1L)
chr2<-data.frame(matrix(unlist(chr1),nrow=length(chr1),byrow=T))[,1]
#datanrow = do.call(rbind,lapply(filelist,function(fn)nrow(read.table(fn,header=FALSE, sep="\t"))))
#all 37500
Scaffold<-rep(as.character(chr2),each=37500)
colnames(datafr)<-colnames(scaff1)
datafr2<-rbind(scaff1,datafr)
dataSweed_Lang<-data.frame(Scaffold,datafr2)
names(dataSweed_Lang)<-c("Scaffold","Position","Likelihood","Alpha")

SweedLB<-dataSweed_Lang


test2LB<-as.character(tolower(SweedLB$Scaffold))
SweedLB$Scaffold<-test2LB
for (j in 1:nrow(Candidates))
	{nam <- paste("SweedofinterestLB", j, sep ="")
 	AF1LB<-SweedLB[SweedLB$Scaffold==Candidates$Chr[j],]
	interval<-which(AF1LB$Position>=(Candidates$gene_start[j]-windowsize[j])&AF1LB$Position<=(Candidates$gene_end[j]+windowsize[j]))
	interval2<-c(min(interval)-1,interval,max(interval)+1)
	SweedofinterestLB<-AF1LB[interval2,]
	assign(nam, SweedofinterestLB)
	}
rm(SweedofinterestLB)
list<-ls(pattern="^SweedofinterestLB")

test2Best<-as.character(tolower(Sweed_Best$Scaffold))
Sweed_Best$Scaffold<-test2Best
for (j in 1:nrow(Candidates))
	{nam <- paste("Sweed_Bestofinterest", j, sep ="")
 	AF1Best<-Sweed_Best[Sweed_Best$Scaffold==Candidates$Chr[j],]
	Sweed_Bestofinterest<-AF1Best[AF1Best$Position>=(Candidates$gene_start[j]-windowsize[j]-100000)&AF1Best$Position<=(Candidates$gene_end[j]+windowsize[j]+100000),]
	assign(nam, Sweed_Bestofinterest)
	}
rm(Sweed_Bestofinterest)
list<-ls(pattern="^Sweed_Bestofinterest")
setwd("F:/Cu_project/Genome_scans/")
d<-Candidates
for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("Overlap_any_Metrics_halleri_combined_plot_black",j,"_",Candidates$Ah_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=8, height=10,paper="special")
	par(mfcol=c(8,3))
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,5,3,0))
	par(mgp=c(5,1,0))


	plot(get(paste("AFofinterestKK",j,sep=""))$DD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$DD),max(AFdataKK$DD)+0.15),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataKK$DD-0.3),Candidates$gene_end[j],max(AFdataKK$DD)+0.4,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$DD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="purple")
	abline(h=DD_abs_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$DD)+0.12,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdataKK$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$DD)+0.05,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Fst~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Fst),max(AFdataKK$Fst)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataKK$Fst)-0.7,Candidates$gene_end[j],max(AFdataKK$Fst)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Fst~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="hotpink")
	abline(h=Fst_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Fst)+0.55,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdataKK$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Fst)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Nielsen~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Nielsen),max(AFdataKK$Nielsen)+150),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataKK$Nielsen)-200,Candidates$gene_end[j],max(AFdataKK$Nielsen)+200,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Nielsen~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="red")
	abline(h=Nielsen_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Nielsen)+125,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdataKK$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Nielsen)+25,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Dxy~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Dxy),max(AFdataKK$Dxy)+0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkorange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkorange",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataKK$Dxy-0.9),Candidates$gene_end[j],max(AFdataKK$Dxy)+0.9,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Dxy~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="darkorange")
	abline(h=Dxy_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Dxy)+0.7,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdataKK$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Dxy)+0.3,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$VarLD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$VarLD),max(AFdataKK$VarLD)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow3",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataKK$VarLD)-5,Candidates$gene_end[j],max(AFdataKK$VarLD)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$VarLD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="yellow3")
	abline(h=VarLD_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$VarLD)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$VarLD)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$VarLD)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$VarLD)+2,x1=Candidates$gene_end[j],y1=max(AFdataKK$VarLD)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$VarLD)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$VarLD)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$VarLD)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$VarLD)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Flk~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Flk),max(AFdataKK$Flk)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="blue",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataKK$Flk)-5,Candidates$gene_end[j],max(AFdataKK$Flk)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Flk~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="blue")
	abline(h=Flk_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Flk)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdataKK$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Flk)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("TDKromofinterest",j,sep=""))$TajimasDKrom~get(paste("TDKromofinterest",j,sep=""))$winmiddleTDKrom,ylab="",ylim=c(min(TD_Krom$TajimasDKrom),max(TD_Krom$TajimasDKrom)+3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen",xaxt="n")
	rect(Candidates$gene_start[j],min(TD_Krom$TajimasDKrom)-60,Candidates$gene_end[j],max(TD_Krom$TajimasDKrom)+60,col="lightgrey",border = NA)	
	lines(get(paste("TDKromofinterest",j,sep=""))$TajimasDKrom~get(paste("TDKromofinterest",j,sep=""))$winmiddleTDKrom,lwd=2,col="darkgreen")
	lines(get(paste("TDKosiofinterest",j,sep=""))$TajimasDKosi~get(paste("TDKosiofinterest",j,sep=""))$winmiddleTDKosi,lwd=2,col="lightgreen")
#	text(x=genemiddle[[1]][j],y=max(TD_Krom$TajimasDKrom)+2,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Krom$TajimasDKrom)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Krom$TajimasDKrom)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Krom$TajimasDKrom)+1,x1=Candidates$gene_end[j],y1=max(TD_Krom$TajimasDKrom)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Krom$TajimasDKrom)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Krom$TajimasDKrom)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Krom$TajimasDKrom)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Krom$TajimasDKrom)+1,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("SweedofinterestKK",j,sep=""))$Likelihood~get(paste("SweedofinterestKK",j,sep=""))$Position,ylab="",ylim=c(min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4",xaxt="n")
	rect(Candidates$gene_start[j],min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood,get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)-60,Candidates$gene_end[j],max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood,get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)+60,col="lightgrey",border = NA)	
	lines(get(paste("SweedofinterestKK",j,sep=""))$Likelihood~get(paste("SweedofinterestKK",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kosiofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
#	text(x=genemiddle[[1]][j],y=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)


	plot(get(paste("AFofinterestNP",j,sep=""))$DD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$DD),max(AFdataNP$DD)+0.15),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataNP$DD-0.3),Candidates$gene_end[j],max(AFdataNP$DD)+0.4,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$DD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="purple")
	abline(h=DD_abs_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$DD)+0.12,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdataNP$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$DD)+0.05,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Fst~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Fst),max(AFdataNP$Fst)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataNP$Fst)-0.7,Candidates$gene_end[j],max(AFdataNP$Fst)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Fst~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="hotpink")
	abline(h=Fst_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Fst)+0.55,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdataNP$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Fst)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Nielsen~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Nielsen),max(AFdataNP$Nielsen)+150),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataNP$Nielsen)-200,Candidates$gene_end[j],max(AFdataNP$Nielsen)+200,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Nielsen~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="red")
	abline(h=Nielsen_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Nielsen)+125,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdataNP$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Nielsen)+25,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Dxy~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Dxy),max(AFdataNP$Dxy)+0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkorange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkorange",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataNP$Dxy-0.9),Candidates$gene_end[j],max(AFdataNP$Dxy)+0.9,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Dxy~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="darkorange")
	abline(h=Dxy_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Dxy)+0.7,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdataNP$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Dxy)+0.3,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$VarLD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$VarLD),max(AFdataNP$VarLD)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow3",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataNP$VarLD)-5,Candidates$gene_end[j],max(AFdataNP$VarLD)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$VarLD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="yellow3")
	abline(h=VarLD_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$VarLD)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$VarLD)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$VarLD)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$VarLD)+2,x1=Candidates$gene_end[j],y1=max(AFdataNP$VarLD)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$VarLD)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$VarLD)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$VarLD)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$VarLD)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Flk~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Flk),max(AFdataNP$Flk)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="blue",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataNP$Flk)-5,Candidates$gene_end[j],max(AFdataNP$Flk)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Flk~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="blue")
	abline(h=Flk_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Flk)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdataNP$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Flk)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDNossofinterest",j,sep=""))$TajimasDNoss~get(paste("TDNossofinterest",j,sep=""))$winmiddleTDNoss,ylab="",ylim=c(min(TD_Noss$TajimasDNoss),max(TD_Noss$TajimasDNoss)+3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen",xaxt="n")
	rect(Candidates$gene_start[j],min(TD_Noss$TajimasDNoss)-60,Candidates$gene_end[j],max(TD_Noss$TajimasDNoss)+60,col="lightgrey",border = NA)	
	lines(get(paste("TDNossofinterest",j,sep=""))$TajimasDNoss~get(paste("TDNossofinterest",j,sep=""))$winmiddleTDNoss,lwd=2,col="darkgreen")
	lines(get(paste("TDPaisofinterest",j,sep=""))$TajimasDPais~get(paste("TDPaisofinterest",j,sep=""))$winmiddleTDPais,lwd=2,col="lightgreen")
#	text(x=genemiddle[[1]][j],y=max(TD_Noss$TajimasDNoss)+2,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Noss$TajimasDNoss)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Noss$TajimasDNoss)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Noss$TajimasDNoss)+1,x1=Candidates$gene_end[j],y1=max(TD_Noss$TajimasDNoss)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Noss$TajimasDNoss)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Noss$TajimasDNoss)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Noss$TajimasDNoss)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Noss$TajimasDNoss)+1,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("SweedofinterestNP",j,sep=""))$Likelihood~get(paste("SweedofinterestNP",j,sep=""))$Position,ylab="",ylim=c(min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4",xaxt="n")
	rect(Candidates$gene_start[j],min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood,get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)-60,Candidates$gene_end[j],max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood,get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)+60,col="lightgrey",border = NA)	
	lines(get(paste("SweedofinterestNP",j,sep=""))$Likelihood~get(paste("SweedofinterestNP",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Paisofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
#	text(x=genemiddle[[1]][j],y=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)


	plot(get(paste("AFofinterestLB",j,sep=""))$DD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$DD),max(AFdataLB$DD)+0.15),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataLB$DD-0.3),Candidates$gene_end[j],max(AFdataLB$DD)+0.4,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$DD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="purple")
	abline(h=DD_abs_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$DD)+0.12,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdataLB$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$DD)+0.05,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Fst~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Fst),max(AFdataLB$Fst)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataLB$Fst)-0.7,Candidates$gene_end[j],max(AFdataLB$Fst)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Fst~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="hotpink")
	abline(h=Fst_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Fst)+0.55,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdataLB$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Fst)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Nielsen~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Nielsen),max(AFdataLB$Nielsen)+150),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataLB$Nielsen)-200,Candidates$gene_end[j],max(AFdataLB$Nielsen)+200,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Nielsen~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="red")
	abline(h=Nielsen_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Nielsen)+125,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdataLB$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Nielsen)+25,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Dxy~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Dxy),max(AFdataLB$Dxy)+0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkorange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkorange",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataLB$Dxy-0.9),Candidates$gene_end[j],max(AFdataLB$Dxy)+0.9,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Dxy~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="darkorange")
	abline(h=Dxy_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Dxy)+0.7,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdataLB$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Dxy)+0.3,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$VarLD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$VarLD),max(AFdataLB$VarLD)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow3",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataLB$VarLD)-5,Candidates$gene_end[j],max(AFdataLB$VarLD)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$VarLD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="yellow3")
	abline(h=VarLD_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$VarLD)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$VarLD)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$VarLD)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$VarLD)+2,x1=Candidates$gene_end[j],y1=max(AFdataLB$VarLD)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$VarLD)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$VarLD)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$VarLD)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$VarLD)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Flk~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Flk),max(AFdataLB$Flk)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="blue",xaxt="n")
	rect(Candidates$gene_start[j],min(AFdataLB$Flk)-5,Candidates$gene_end[j],max(AFdataLB$Flk)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Flk~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="blue")
	abline(h=Flk_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Flk)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdataLB$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Flk)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}




	plot(get(paste("TDLangofinterest",j,sep=""))$TajimasDLang~get(paste("TDLangofinterest",j,sep=""))$winmiddleTDLang,ylab="",ylim=c(min(TD_Lang$TajimasDLang),max(TD_Lang$TajimasDLang)+3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen",xaxt="n")
	rect(Candidates$gene_start[j],min(TD_Lang$TajimasDLang)-60,Candidates$gene_end[j],max(TD_Lang$TajimasDLang)+60,col="lightgrey",border = NA)	
	lines(get(paste("TDLangofinterest",j,sep=""))$TajimasDLang~get(paste("TDLangofinterest",j,sep=""))$winmiddleTDLang,lwd=2,col="darkgreen")
	lines(get(paste("TDBestofinterest",j,sep=""))$TajimasDBest~get(paste("TDBestofinterest",j,sep=""))$winmiddleTDBest,lwd=2,col="lightgreen")
#	text(x=genemiddle[[1]][j],y=max(TD_Lang$TajimasDLang)+2,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Lang$TajimasDLang)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Lang$TajimasDLang)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Lang$TajimasDLang)+1,x1=Candidates$gene_end[j],y1=max(TD_Lang$TajimasDLang)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Lang$TajimasDLang)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Lang$TajimasDLang)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Lang$TajimasDLang)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Lang$TajimasDLang)+1,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("SweedofinterestLB",j,sep=""))$Likelihood~get(paste("SweedofinterestLB",j,sep=""))$Position,ylab="",ylim=c(min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4",xaxt="n")
	rect(Candidates$gene_start[j],min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood,get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)-60,Candidates$gene_end[j],max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood,get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)+60,col="lightgrey",border = NA)	
	lines(get(paste("SweedofinterestLB",j,sep=""))$Likelihood~get(paste("SweedofinterestLB",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Bestofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
#	text(x=genemiddle[[1]][j],y=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)



	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(DD)),side=2,line=1,outer=TRUE,cex=1.3,adj=0.95,col="purple")
	mtext(text=expression(bold(F[ST])),side=2,line=1,outer=TRUE,cex=1.3,adj=0.82,col="hotpink")
	mtext(text=expression(bold("2dSFS")),side=2,line=1,outer=TRUE,cex=1.3,adj=0.7,col="red")
	mtext(text=expression(bold(d[XY])),side=2,line=1,outer=TRUE,cex=1.3,adj=0.55,col="darkorange")
	mtext(text=expression(bold(Flk)),side=2,line=1,outer=TRUE,cex=1.3,adj=0.3,col="blue")
	mtext(text=expression(bold("VarLD")),side=2,line=1,outer=TRUE,cex=1.3,adj=0.42,col="yellow3")
	mtext(text=expression(bold("TD M")),side=2,line=1,outer=TRUE,cex=1.3,adj=0.16,col="darkgreen")
	mtext(text=expression(bold("SweeD M")),side=2,line=1,outer=TRUE,cex=1.3,adj=0,col="lightsalmon4")
	mtext(text=expression(bold("TD NM")),side=2,line=3,outer=TRUE,cex=1.3,adj=0.16,col="lightgreen")
	mtext(text=expression(bold("SweeD NM")),side=2,line=3,outer=TRUE,cex=1.3,adj=0,col="lightsalmon")
	mtext(text=expression(bold("Krom-Kosi")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.1)
	mtext(text=expression(bold("Noss-Pais")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	mtext(text=expression(bold("Lang-Best")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.9)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}
###############################################################
#KromKosi only#
###############################################################

for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("KK_any_Metrics_halleri_combined_plot_black",j,"_",Candidates$Ah_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=7, height=10,paper="special")
	par(mfcol=c(8,1))
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,6,3,0))
	par(mgp=c(5,1,0))


	plot(get(paste("AFofinterestKK",j,sep=""))$DD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$DD),max(AFdataKK$DD)+0.15),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestKK",j,sep=""))$DD),max(get(paste("AFofinterestKK",j,sep=""))$DD+0.15),2),1))
	rect(Candidates$gene_start[j],min(AFdataKK$DD-0.3),Candidates$gene_end[j],max(AFdataKK$DD)+0.4,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$DD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="purple")
	abline(h=DD_abs_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$DD)+0.12,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdataKK$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$DD)+0.05,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Fst~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Fst),max(AFdataKK$Fst)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestKK",j,sep=""))$Fst),max(get(paste("AFofinterestKK",j,sep=""))$Fst+0.6),2),0))
	rect(Candidates$gene_start[j],min(AFdataKK$Fst)-0.7,Candidates$gene_end[j],max(AFdataKK$Fst)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Fst~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="hotpink")
	abline(h=Fst_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Fst)+0.55,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdataKK$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Fst)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Nielsen~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Nielsen),max(AFdataKK$Nielsen)+150),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestKK",j,sep=""))$Nielsen),max(get(paste("AFofinterestKK",j,sep=""))$Nielsen+150),2),0))
	rect(Candidates$gene_start[j],min(AFdataKK$Nielsen)-200,Candidates$gene_end[j],max(AFdataKK$Nielsen)+200,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Nielsen~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="red")
	abline(h=Nielsen_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Nielsen)+125,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdataKK$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Nielsen)+25,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Dxy~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Dxy),max(AFdataKK$Dxy)+0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkorange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkorange",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestKK",j,sep=""))$Dxy),max(get(paste("AFofinterestKK",j,sep=""))$Dxy+0.8),2),0))
	rect(Candidates$gene_start[j],min(AFdataKK$Dxy-0.9),Candidates$gene_end[j],max(AFdataKK$Dxy)+0.9,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Dxy~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="darkorange")
	abline(h=Dxy_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Dxy)+0.7,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdataKK$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Dxy)+0.3,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$VarLD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$VarLD),max(AFdataKK$VarLD)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow3",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestKK",j,sep=""))$VarLD),max(get(paste("AFofinterestKK",j,sep=""))$VarLD+0.5),2),0))
	rect(Candidates$gene_start[j],min(AFdataKK$VarLD)-5,Candidates$gene_end[j],max(AFdataKK$VarLD)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$VarLD~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="yellow3")
	abline(h=VarLD_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$VarLD)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$VarLD)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$VarLD)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$VarLD)+2,x1=Candidates$gene_end[j],y1=max(AFdataKK$VarLD)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$VarLD)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$VarLD)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$VarLD)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$VarLD)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestKK",j,sep=""))$Flk~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,ylab="",ylim=c(min(AFdataKK$Flk),max(AFdataKK$Flk)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="blue",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestKK",j,sep=""))$Flk),max(get(paste("AFofinterestKK",j,sep=""))$Flk+5),2),0))
	rect(Candidates$gene_start[j],min(AFdataKK$Flk)-5,Candidates$gene_end[j],max(AFdataKK$Flk)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestKK",j,sep=""))$Flk~get(paste("AFofinterestKK",j,sep=""))$winmiddleKK,lwd=2,col="blue")
	abline(h=Flk_percKK,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataKK$Flk)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataKK$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdataKK$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataKK$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataKK$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataKK$Flk)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}




	plot(get(paste("TDKromofinterest",j,sep=""))$TajimasDKrom~get(paste("TDKromofinterest",j,sep=""))$winmiddleTDKrom,ylab="",ylim=c(min(TD_Krom$TajimasDKrom),max(TD_Krom$TajimasDKrom)+3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen",xaxt="n",las=2,yaxp=round(c(min(TD_Krom$TajimasDKrom),max(TD_Krom$TajimasDKrom)+3,2),0))
	rect(Candidates$gene_start[j],min(TD_Krom$TajimasDKrom)-60,Candidates$gene_end[j],max(TD_Krom$TajimasDKrom)+60,col="lightgrey",border = NA)	
	lines(get(paste("TDKromofinterest",j,sep=""))$TajimasDKrom~get(paste("TDKromofinterest",j,sep=""))$winmiddleTDKrom,lwd=2,col="darkgreen")
	lines(get(paste("TDKosiofinterest",j,sep=""))$TajimasDKosi~get(paste("TDKosiofinterest",j,sep=""))$winmiddleTDKosi,lwd=2,col="lightgreen")
#	text(x=genemiddle[[1]][j],y=max(TD_Krom$TajimasDKrom)+2,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Krom$TajimasDKrom)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Krom$TajimasDKrom)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Krom$TajimasDKrom)+1,x1=Candidates$gene_end[j],y1=max(TD_Krom$TajimasDKrom)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Krom$TajimasDKrom)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Krom$TajimasDKrom)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Krom$TajimasDKrom)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Krom$TajimasDKrom)+1,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("SweedofinterestKK",j,sep=""))$Likelihood~get(paste("SweedofinterestKK",j,sep=""))$Position,ylab="",ylim=c(min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4",xaxt="n",las=2,yaxp=round(c(min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/4,2),0))

	rect(Candidates$gene_start[j],min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood,get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)-60,Candidates$gene_end[j],max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood,get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)+60,col="lightgrey",border = NA)	
	lines(get(paste("SweedofinterestKK",j,sep=""))$Likelihood~get(paste("SweedofinterestKK",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Kosiofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
#	text(x=genemiddle[[1]][j],y=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),max(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestKK",j,sep=""))$Likelihood),min(get(paste("Sweed_Kosiofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)



	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(DD)),side=2,line=2,outer=TRUE,cex=1.3,adj=0.95,col="purple")
	mtext(text=expression(bold(F[ST])),side=2,line=2,outer=TRUE,cex=1.3,adj=0.82,col="hotpink")
	mtext(text=expression(bold("2dSFS")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.7,col="red")
	mtext(text=expression(bold(d[XY])),side=2,line=2,outer=TRUE,cex=1.3,adj=0.55,col="darkorange")
	mtext(text=expression(bold(Flk)),side=2,line=2,outer=TRUE,cex=1.3,adj=0.3,col="blue")
	mtext(text=expression(bold("VarLD")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.42,col="yellow3")
	mtext(text=expression(bold("TD M")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.16,col="darkgreen")
	mtext(text=expression(bold("SweeD M")),side=2,line=2,outer=TRUE,cex=1.3,adj=0,col="lightsalmon4")
	mtext(text=expression(bold("TD NM")),side=2,line=4,outer=TRUE,cex=1.3,adj=0.16,col="lightgreen")
	mtext(text=expression(bold("SweeD NM")),side=2,line=4,outer=TRUE,cex=1.3,adj=0,col="lightsalmon")
	mtext(text=expression(bold("Krom-Kosi")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}




#########################################################################
#NossPais only#
#########################################################################
for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("NP_any_Metrics_halleri_combined_plot_black",j,"_",Candidates$Ah_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=7, height=10,paper="special")
	par(mfcol=c(8,1))
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,6,3,0))
	par(mgp=c(5,1,0))



	plot(get(paste("AFofinterestNP",j,sep=""))$DD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$DD),max(AFdataNP$DD)+0.15),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestNP",j,sep=""))$DD),max(get(paste("AFofinterestNP",j,sep=""))$DD+0.15),2),1))
	rect(Candidates$gene_start[j],min(AFdataNP$DD-0.3),Candidates$gene_end[j],max(AFdataNP$DD)+0.4,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$DD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="purple")
	abline(h=DD_abs_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$DD)+0.12,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdataNP$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$DD)+0.05,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Fst~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Fst),max(AFdataNP$Fst)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestNP",j,sep=""))$Fst),max(get(paste("AFofinterestNP",j,sep=""))$Fst+0.6),2),0))
	rect(Candidates$gene_start[j],min(AFdataNP$Fst)-0.7,Candidates$gene_end[j],max(AFdataNP$Fst)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Fst~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="hotpink")
	abline(h=Fst_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Fst)+0.55,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdataNP$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Fst)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Nielsen~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Nielsen),max(AFdataNP$Nielsen)+150),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestNP",j,sep=""))$Nielsen),max(get(paste("AFofinterestNP",j,sep=""))$Nielsen+150),2),0))
	rect(Candidates$gene_start[j],min(AFdataNP$Nielsen)-200,Candidates$gene_end[j],max(AFdataNP$Nielsen)+200,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Nielsen~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="red")
	abline(h=Nielsen_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Nielsen)+125,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdataNP$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Nielsen)+25,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Dxy~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Dxy),max(AFdataNP$Dxy)+0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkorange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkorange",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestNP",j,sep=""))$Dxy),max(get(paste("AFofinterestNP",j,sep=""))$Dxy+0.8),2),0))
	rect(Candidates$gene_start[j],min(AFdataNP$Dxy-0.9),Candidates$gene_end[j],max(AFdataNP$Dxy)+0.9,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Dxy~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="darkorange")
	abline(h=Dxy_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Dxy)+0.7,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdataNP$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Dxy)+0.3,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$VarLD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$VarLD),max(AFdataNP$VarLD)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow3",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestNP",j,sep=""))$VarLD),max(get(paste("AFofinterestNP",j,sep=""))$VarLD+0.5),2),0))
	rect(Candidates$gene_start[j],min(AFdataNP$VarLD)-5,Candidates$gene_end[j],max(AFdataNP$VarLD)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$VarLD~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="yellow3")
	abline(h=VarLD_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$VarLD)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$VarLD)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$VarLD)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$VarLD)+2,x1=Candidates$gene_end[j],y1=max(AFdataNP$VarLD)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$VarLD)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$VarLD)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$VarLD)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$VarLD)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestNP",j,sep=""))$Flk~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,ylab="",ylim=c(min(AFdataNP$Flk),max(AFdataNP$Flk)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="blue",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestNP",j,sep=""))$Flk),max(get(paste("AFofinterestNP",j,sep=""))$Flk+5),2),0))
	rect(Candidates$gene_start[j],min(AFdataNP$Flk)-5,Candidates$gene_end[j],max(AFdataNP$Flk)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestNP",j,sep=""))$Flk~get(paste("AFofinterestNP",j,sep=""))$winmiddleNP,lwd=2,col="blue")
	abline(h=Flk_percNP,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataNP$Flk)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataNP$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdataNP$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataNP$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataNP$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataNP$Flk)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}



	plot(get(paste("TDNossofinterest",j,sep=""))$TajimasDNoss~get(paste("TDNossofinterest",j,sep=""))$winmiddleTDNoss,ylab="",ylim=c(min(TD_Noss$TajimasDNoss),max(TD_Noss$TajimasDNoss)+3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen",xaxt="n",las=2,yaxp=round(c(min(TD_Noss$TajimasDNoss),max(TD_Noss$TajimasDNoss)+3,2),0))
	rect(Candidates$gene_start[j],min(TD_Noss$TajimasDNoss)-60,Candidates$gene_end[j],max(TD_Noss$TajimasDNoss)+60,col="lightgrey",border = NA)	
	lines(get(paste("TDNossofinterest",j,sep=""))$TajimasDNoss~get(paste("TDNossofinterest",j,sep=""))$winmiddleTDNoss,lwd=2,col="darkgreen")
	lines(get(paste("TDPaisofinterest",j,sep=""))$TajimasDPais~get(paste("TDPaisofinterest",j,sep=""))$winmiddleTDPais,lwd=2,col="lightgreen")
#	text(x=genemiddle[[1]][j],y=max(TD_Noss$TajimasDNoss)+2,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Noss$TajimasDNoss)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Noss$TajimasDNoss)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Noss$TajimasDNoss)+1,x1=Candidates$gene_end[j],y1=max(TD_Noss$TajimasDNoss)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Noss$TajimasDNoss)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Noss$TajimasDNoss)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Noss$TajimasDNoss)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Noss$TajimasDNoss)+1,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("SweedofinterestNP",j,sep=""))$Likelihood~get(paste("SweedofinterestNP",j,sep=""))$Position,ylab="",ylim=c(min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4",xaxt="n",las=2,yaxp=round(c(min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/4,2),0))

	rect(Candidates$gene_start[j],min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood,get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)-60,Candidates$gene_end[j],max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood,get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)+60,col="lightgrey",border = NA)	
	lines(get(paste("SweedofinterestNP",j,sep=""))$Likelihood~get(paste("SweedofinterestNP",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Paisofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
#	text(x=genemiddle[[1]][j],y=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),max(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestNP",j,sep=""))$Likelihood),min(get(paste("Sweed_Paisofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)


	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(DD)),side=2,line=2,outer=TRUE,cex=1.3,adj=0.95,col="purple")
	mtext(text=expression(bold(F[ST])),side=2,line=2,outer=TRUE,cex=1.3,adj=0.82,col="hotpink")
	mtext(text=expression(bold("2dSFS")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.7,col="red")
	mtext(text=expression(bold(d[XY])),side=2,line=2,outer=TRUE,cex=1.3,adj=0.55,col="darkorange")
	mtext(text=expression(bold(Flk)),side=2,line=2,outer=TRUE,cex=1.3,adj=0.3,col="blue")
	mtext(text=expression(bold("VarLD")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.42,col="yellow3")
	mtext(text=expression(bold("TD M")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.16,col="darkgreen")
	mtext(text=expression(bold("SweeD M")),side=2,line=2,outer=TRUE,cex=1.3,adj=0,col="lightsalmon4")
	mtext(text=expression(bold("TD NM")),side=2,line=4,outer=TRUE,cex=1.3,adj=0.16,col="lightgreen")
	mtext(text=expression(bold("SweeD NM")),side=2,line=4,outer=TRUE,cex=1.3,adj=0,col="lightsalmon")
	mtext(text=expression(bold("Noss-Pais")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}



########################################################################
#LangBest only#
########################################################################
for (j in 1:nrow(d))
	{
	arrowdir<-ifelse(get(paste("Genesofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("LB_any_Metrics_halleri_combined_plot_black",j,"_",Candidates$Ah_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=7, height=10,paper="special")
	par(mfcol=c(8,1))
	par(mar=c(0,2,1,1)+0.1)
	par(oma=c(5,6,3,0))
	par(mgp=c(5,1,0))



	plot(get(paste("AFofinterestLB",j,sep=""))$DD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$DD),max(AFdataLB$DD)+0.15),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="purple",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="purple",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestLB",j,sep=""))$DD),max(get(paste("AFofinterestLB",j,sep=""))$DD+0.15),2),1))
	rect(Candidates$gene_start[j],min(AFdataLB$DD-0.3),Candidates$gene_end[j],max(AFdataLB$DD)+0.4,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$DD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="purple")
	abline(h=DD_abs_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$DD)+0.12,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$DD)+0.05,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$DD)+0.05,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$DD)+0.05,x1=Candidates$gene_end[j],y1=max(AFdataLB$DD)+0.05,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$DD)+0.05,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$DD)+0.05,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$DD)+0.05,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$DD)+0.05,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Fst~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Fst),max(AFdataLB$Fst)+0.6),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="hotpink",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="hotpink",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestLB",j,sep=""))$Fst),max(get(paste("AFofinterestLB",j,sep=""))$Fst+0.6),2),0))
	rect(Candidates$gene_start[j],min(AFdataLB$Fst)-0.7,Candidates$gene_end[j],max(AFdataLB$Fst)+0.7,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Fst~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="hotpink")
	abline(h=Fst_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Fst)+0.55,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Fst)+0.2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Fst)+0.2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Fst)+0.2,x1=Candidates$gene_end[j],y1=max(AFdataLB$Fst)+0.2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Fst)+0.2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Fst)+0.2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Fst)+0.2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Fst)+0.2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Nielsen~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Nielsen),max(AFdataLB$Nielsen)+150),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="red",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="red",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestLB",j,sep=""))$Nielsen),max(get(paste("AFofinterestLB",j,sep=""))$Nielsen+150),2),0))
	rect(Candidates$gene_start[j],min(AFdataLB$Nielsen)-200,Candidates$gene_end[j],max(AFdataLB$Nielsen)+200,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Nielsen~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="red")
	abline(h=Nielsen_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Nielsen)+125,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Nielsen)+25,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Nielsen)+25,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Nielsen)+25,x1=Candidates$gene_end[j],y1=max(AFdataLB$Nielsen)+25,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Nielsen)+25,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Nielsen)+25,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Nielsen)+25,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Nielsen)+25,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Dxy~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Dxy),max(AFdataLB$Dxy)+0.8),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkorange",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkorange",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestLB",j,sep=""))$Dxy),max(get(paste("AFofinterestLB",j,sep=""))$Dxy+0.8),2),0))
	rect(Candidates$gene_start[j],min(AFdataLB$Dxy-0.9),Candidates$gene_end[j],max(AFdataLB$Dxy)+0.9,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Dxy~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="darkorange")
	abline(h=Dxy_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Dxy)+0.7,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Dxy)+0.3,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Dxy)+0.3,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Dxy)+0.3,x1=Candidates$gene_end[j],y1=max(AFdataLB$Dxy)+0.3,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Dxy)+0.3,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Dxy)+0.3,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Dxy)+0.3,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Dxy)+0.3,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$VarLD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$VarLD),max(AFdataLB$VarLD)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="yellow3",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestLB",j,sep=""))$VarLD),max(get(paste("AFofinterestLB",j,sep=""))$VarLD+0.5),2),0))
	rect(Candidates$gene_start[j],min(AFdataLB$VarLD)-5,Candidates$gene_end[j],max(AFdataLB$VarLD)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$VarLD~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="yellow3")
	abline(h=VarLD_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$VarLD)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$VarLD)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$VarLD)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$VarLD)+2,x1=Candidates$gene_end[j],y1=max(AFdataLB$VarLD)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$VarLD)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$VarLD)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$VarLD)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$VarLD)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}

	plot(get(paste("AFofinterestLB",j,sep=""))$Flk~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,ylab="",ylim=c(min(AFdataLB$Flk),max(AFdataLB$Flk)+5),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="blue",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="blue",xaxt="n",las=2,yaxp=round(c(min(get(paste("AFofinterestLB",j,sep=""))$Flk),max(get(paste("AFofinterestLB",j,sep=""))$Flk+5),2),0))
	rect(Candidates$gene_start[j],min(AFdataLB$Flk)-5,Candidates$gene_end[j],max(AFdataLB$Flk)+8,col="lightgrey",border = NA)	
	lines(get(paste("AFofinterestLB",j,sep=""))$Flk~get(paste("AFofinterestLB",j,sep=""))$winmiddleLB,lwd=2,col="blue")
	abline(h=Flk_percLB,lty=2,col="black")
#	text(x=genemiddle[[1]][j],y=max(AFdataLB$Flk)+4,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Flk)+2,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Flk)+2,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(AFdataLB$Flk)+2,x1=Candidates$gene_end[j],y1=max(AFdataLB$Flk)+2,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Flk)+2,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(AFdataLB$Flk)+2,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(AFdataLB$Flk)+2,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(AFdataLB$Flk)+2,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("TDLangofinterest",j,sep=""))$TajimasDLang~get(paste("TDLangofinterest",j,sep=""))$winmiddleTDLang,ylab="",ylim=c(min(TD_Lang$TajimasDLang),max(TD_Lang$TajimasDLang)+3),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="darkgreen",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="darkgreen",xaxt="n",las=2,yaxp=round(c(min(TD_Lang$TajimasDLang),max(TD_Lang$TajimasDLang)+3,2),0))
	rect(Candidates$gene_start[j],min(TD_Lang$TajimasDLang)-60,Candidates$gene_end[j],max(TD_Lang$TajimasDLang)+60,col="lightgrey",border = NA)	
	lines(get(paste("TDLangofinterest",j,sep=""))$TajimasDLang~get(paste("TDLangofinterest",j,sep=""))$winmiddleTDLang,lwd=2,col="darkgreen")
	lines(get(paste("TDBestofinterest",j,sep=""))$TajimasDBest~get(paste("TDBestofinterest",j,sep=""))$winmiddleTDBest,lwd=2,col="lightgreen")
#	text(x=genemiddle[[1]][j],y=max(TD_Lang$TajimasDLang)+2,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Lang$TajimasDLang)+1,x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(TD_Lang$TajimasDLang)+1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(TD_Lang$TajimasDLang)+1,x1=Candidates$gene_end[j],y1=max(TD_Lang$TajimasDLang)+1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(TD_Lang$TajimasDLang)+1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(TD_Lang$TajimasDLang)+1,code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(TD_Lang$TajimasDLang)+1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(TD_Lang$TajimasDLang)+1,code=arrowdir2[j],length=0,col="red",lwd=3)
		}


	plot(get(paste("SweedofinterestLB",j,sep=""))$Likelihood~get(paste("SweedofinterestLB",j,sep=""))$Position,ylab="",ylim=c(min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/4),cex=1.5,cex.lab=1.5,cex.axis=1.5,type="l",lwd=2,col="lightsalmon4",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),col.lab="lightsalmon4",xaxt="n",las=2,yaxp=round(c(min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)),
	+max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),
+min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/4,2),0))

	rect(Candidates$gene_start[j],min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood,get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)-60,Candidates$gene_end[j],max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood,get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)+60,col="lightgrey",border = NA)	
	lines(get(paste("SweedofinterestLB",j,sep=""))$Likelihood~get(paste("SweedofinterestLB",j,sep=""))$Position,lwd=2,col="lightsalmon4")
	lines(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood~get(paste("Sweed_Bestofinterest",j,sep=""))$Position,lwd=2,col="lightsalmon")
#	text(x=genemiddle[[1]][j],y=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/5),labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Genesofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Genesofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Genesofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),x1=Candidates$gene_end[j],
	y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="grey",lwd=3)
		}
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),
		x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))+((max(max(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),max(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood))-min(min(get(paste("SweedofinterestLB",j,sep=""))$Likelihood),min(get(paste("Sweed_Bestofinterest",j,sep=""))$Likelihood)))/12),code=1,length=0,col="red",lwd=3)
		}

	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)



	mtext(text=expression(bold(Scaffold~position~(kb))),side=1,line=4,outer=TRUE,cex=1.5)

	mtext(text=expression(bold(DD)),side=2,line=2,outer=TRUE,cex=1.3,adj=0.95,col="purple")
	mtext(text=expression(bold(F[ST])),side=2,line=2,outer=TRUE,cex=1.3,adj=0.82,col="hotpink")
	mtext(text=expression(bold("2dSFS")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.7,col="red")
	mtext(text=expression(bold(d[XY])),side=2,line=2,outer=TRUE,cex=1.3,adj=0.55,col="darkorange")
	mtext(text=expression(bold(Flk)),side=2,line=2,outer=TRUE,cex=1.3,adj=0.3,col="blue")
	mtext(text=expression(bold("VarLD")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.42,col="yellow3")
	mtext(text=expression(bold("TD M")),side=2,line=2,outer=TRUE,cex=1.3,adj=0.16,col="darkgreen")
	mtext(text=expression(bold("SweeD M")),side=2,line=2,outer=TRUE,cex=1.3,adj=0,col="lightsalmon4")
	mtext(text=expression(bold("TD NM")),side=2,line=4,outer=TRUE,cex=1.3,adj=0.16,col="lightgreen")
	mtext(text=expression(bold("SweeD NM")),side=2,line=4,outer=TRUE,cex=1.3,adj=0,col="lightsalmon")
	mtext(text=expression(bold("Lang-Best")),side=3,line=0,outer=TRUE,cex=1.5,adj=0.5)
	
	axis(side=1,cex=1.5,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex=1.5,cex.lab=1.5,cex.axis=1.5)

	dev.off()
	}


