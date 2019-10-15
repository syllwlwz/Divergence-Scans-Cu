Cnmops
require(cn.mops)

BAMFilesKosi <- list.files(pattern="^Kosi(.)*bam$")
BAMFilesKrom <- list.files(pattern="^Krom(.)*bam$")
BAMFilesLang <- list.files(pattern="^Lan(.)*bam$")
BAMFilesBest <- list.files(pattern="^Best(.)*bam$")
BAMFilesNoss <- list.files(pattern="^Noss(.)*bam$")
BAMFilesPais <- list.files(pattern="^Pais(.)*bam$")

#Window length set to: 500
#C=m*RL/WL
#C: coverage
#m: average number of reads per bin ~ 50-100
#RL: Read length
#WL: window length

#WL=m*RL/C
#WL=100*125/19 ~658
#or WL=100*150/15~1000
#m should be between 50 and 100 https://support.bioconductor.org/p/75969/ 25.01.16

require(Rsamtools)
test<-scanBamHeader(BAMFilesWulm)
scaffolds<-names(test[[1]]$targets)[1:8]

trace("getReadCountsFromBAM",edit=TRUE)
#change sl <- as.integer(unique(sl)), delete unique()

Lang <- getReadCountsFromBAM(BAMFiles=BAMFilesLang,sampleNames=c("Lan3_1","Lang18","Lang21","Lang28","Lang3","Lang4","Lang5","Lange1","Lange12","Lange2","Lange4","Lange6","Lange9"),refSeqName=sort(scaffolds),WL=500) 
Best <- getReadCountsFromBAM(BAMFiles=BAMFilesBest,sampleNames=c("Best10","Best17","Best18","Best19","Best20","Best33","Best36","Best4","Best40","Best41","Best6","Best8"),refSeqName=sort(scaffolds),WL=500) 
LB <- referencecn.mops(cases=Lang,controls=Best, minWidth = 2, minReadCount=2)
resCNMOPSLB <- calcIntegerCopyNumbers(LB)

Noss <- getReadCountsFromBAM(BAMFiles=BAMFilesNoss,sampleNames=c("Noss10","Noss21","Noss22","Noss23","Noss24","Noss25","Noss26","Noss27","Noss28","Noss3","Noss5","Noss6"),refSeqName=sort(scaffolds),WL=500) 
Pais <- getReadCountsFromBAM(BAMFiles=BAMFilesPais,sampleNames=c("Pais10","Pais11","Paos12","Pais22","Pais25","Pais26","Pais27","Pais29","Pais31","Pais4","Pais5","Pais9"),refSeqName=sort(scaffolds),WL=500) 
NP <- referencecn.mops(cases=Noss,controls=Pais, minWidth = 2, minReadCount=2)
resCNMOPSNP <- calcIntegerCopyNumbers(NP)

Krom <- getReadCountsFromBAM(BAMFiles=BAMFilesKrom,sampleNames=c("Krom_002_h04","Krom_002_h05","Krom_002_h06","Krom_002_h07","Krom_002_h08","Krom_002_h10","Krom_003_h08","Krom_003_h10","Krom21","Krom23","Krom24","Krom28"),refSeqName=sort(scaffolds),WL=500) 
Kosi <- getReadCountsFromBAM(BAMFiles=BAMFilesKosi,sampleNames=c("Kosi_002_h01","Kosi_002_h02","Kosi_002_h03","Kosi_002_h05","Kosi_002_h06","Kosi_002_h07","Kosi_002_h08","Kosi_003_h01","Kosi_003_h09","Kosi26","Kosi31","Kosi34"),refSeqName=sort(scaffolds),WL=500) 
KK <- referencecn.mops(cases=Krom,controls=Kosi, minWidth = 2, minReadCount=2)
resCNMOPSKK <- calcIntegerCopyNumbers(KK)

BL <- referencecn.mops(cases=Best,controls=Lang, minWidth = 2, minReadCount=2)
resCNMOPSBL <- calcIntegerCopyNumbers(BL)
PN <- referencecn.mops(cases=Pais,controls=Noss, minWidth = 2, minReadCount=2)
resCNMOPSPN <- calcIntegerCopyNumbers(PN)
KoKr <- referencecn.mops(cases=Kosi,controls=Krom, minWidth = 2, minReadCount=2)
resCNMOPSKoKr <- calcIntegerCopyNumbers(KoKr)


#########################################################
#overlap between individual data and regions to get mean#
#########################################################
#after refrencecn.mops and CalcIntegerCopyNumbers
#########################################################
#overlap between individual data and regions to get mean#
#########################################################
dupl_cov<-function(arg1){
	arg2<-arg1
	arg2[[2]]<-round(arg2[[2]],digits=2)
	arg2[[3]]<-round(arg2[[3]],digits=2)
	if (any(duplicated(arg2[[1]]))){
	duplist<-c(as.character(arg2[[1]][duplicated(arg2[[1]])]))
	for (i in duplist){
		arg2[[2]][arg2[[1]]==i]<-mean(arg2[[2]][arg2[[1]]==i])
		arg2[[3]][arg2[[1]]==i]<-median(arg2[[3]][arg2[[1]]==i])
		arg2[[4]][arg2[[1]]==i]<-names(which.max(table(arg2[[4]][arg2[[1]]==i])))
		}
	arg2[[2]]<-arg2[[2]][!duplicated(arg2[[1]])]
	arg2[[3]]<-arg2[[3]][!duplicated(arg2[[1]])]
	arg2[[4]]<-arg2[[4]][!duplicated(arg2[[1]])]
	arg2[[1]]<-as.character(arg2[[1]][!duplicated(arg2[[1]])])
	}
	return(arg2)
}


#after refrencecn.mops and CalcIntegerCopyNumbers
regionsKK<-cnvr(resCNMOPSKK)
singleKK<-cnvs(resCNMOPSKK)
testKK<-findOverlaps(regionsKK,singleKK)
test2KK<- DataFrame(splitAsList(singleKK$sampleName[subjectHits(testKK)], queryHits(testKK)),splitAsList(singleKK$median[subjectHits(testKK)], queryHits(testKK)),splitAsList(singleKK$mean[subjectHits(testKK)], queryHits(testKK)),splitAsList(singleKK$CN[subjectHits(testKK)], queryHits(testKK)))
colnames(test2KK)<-c("sampleNames","mean","median","CN")
test7KK<-apply(as.data.frame(test2KK),1,dupl_cov)
mcols(regionsKK)<-as.data.frame(as.matrix(test7KK))

dataKK<-as.data.frame(regionsKK)

data_change<-function(arg1){
	arg2<-data.frame(matrix(ncol=6,nrow=0))
	arg2[[1,1]]<-as.character(sapply(arg1[[6]][1],paste0,collapse=",")) 
	arg2[[1,2]]<-sapply(arg1[[6]][2],paste0,collapse=",") 
	arg2[[1,3]]<-sapply(arg1[[6]][3],paste0,collapse=",") 
	arg2[[1,4]]<-sapply(arg1[[6]][4],paste0,collapse=",") 
	arg2[[1,5]]<-names(which.max(table(arg1[[6]][4])))
	arg2[[1,6]]<-max(table(arg1[[6]][4]))/lengths(arg1[[6]][4])
	return(arg2)
}

adddataKK<-do.call("rbind",apply(dataKK,1,data_change))
dataKK<-dataKK[,-6]
frequKK<-cbind(dataKK,adddataKK)

names(frequKK)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KK<-frequKK[count.fields(textConnection(frequKK$sampleNames),sep=",")>=5,]

require(xlsx)
write.xlsx2(frequ3KK,"CnvregionswithmeanhalleriKromKosi_minRead2_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE)
frequbigKK<-read.xlsx2("CnvregionswithmeanhalleriKromKosi_minRead2_WL500_minL2_min5samples.xlsx",1,header=TRUE)
frequbigKK[,2]<-as.numeric(as.character(frequbigKK[,2]))
frequbigKK[,3]<-as.numeric(as.character(frequbigKK[,3]))

#########################################
#find overlaps between windows and genes#
#########################################

library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)

## Load GFF files and add promotor region to genes

### Step 1: import genes from GFF files

gff<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(gff)<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes","Ah_ID")
#gff2<-gff
gff<-gff[gff$Type=="gene",]
gff_GRange<-GRanges(seqnames=tolower(gff$Scaffold),ranges=IRanges(start=gff$Start,end=gff$End))
HalGenes <- gff_GRange
values(HalGenes)<-gff$Ah_ID

###################################################
#find overlaps between windows and A. halleri genes#
###################################################
frequbigKK_GRange<-GRanges(seqnames=tolower(frequbigKK$Scaffold),ranges=IRanges(start=frequbigKK$start,end=frequbigKK$end))
values(frequbigKK_GRange)<-frequbigKK[,6:11]

testKK<-findOverlaps(frequbigKK_GRange,HalGenes)
test2KK<-pintersect(frequbigKK_GRange[queryHits(testKK)],HalGenes[subjectHits(testKK)])
test3KK<-as.data.frame(test2KK)
frequbigKK_halgenes=mergeByOverlaps(frequbigKK_GRange,HalGenes,type=c("any"))
frequbigKK_halgenes_df=data.frame(as.data.frame(frequbigKK_halgenes$frequbigKK_GRange),as.data.frame(frequbigKK_halgenes$HalGenes),test3KK$width)
#frequbigKK_halgenes_df<-frequbigKK_halgenes_df[,-c(6,7)]
colnames(frequbigKK_halgenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKK_halgenes_df<-frequbigKK_halgenes_df[(frequbigKK_halgenes_df$Widthofoverlap>=(0.9*frequbigKK_halgenes_df$gene_size)),]
write.xlsx2(frequbigKK_halgenes_df,"KK_temp.xlsx",col.names=TRUE,row.names=FALSE)


regionsLB<-cnvr(resCNMOPSLB)
singleLB<-cnvs(resCNMOPSLB)
testLB<-findOverlaps(regionsLB,singleLB)
test2LB<- DataFrame(splitAsList(singleLB$sampleName[subjectHits(testLB)], queryHits(testLB)),splitAsList(singleLB$median[subjectHits(testLB)], queryHits(testLB)),splitAsList(singleLB$mean[subjectHits(testLB)], queryHits(testLB)),splitAsList(singleLB$CN[subjectHits(testLB)], queryHits(testLB)))
colnames(test2LB)<-c("sampleNames","mean","median","CN")
test7LB<-apply(as.data.frame(test2LB),1,dupl_cov)
mcols(regionsLB)<-as.data.frame(as.matrix(test7LB))

dataLB<-as.data.frame(regionsLB)

adddataLB<-do.call("rbind",apply(dataLB,1,data_change))
dataLB<-dataLB[,-6]
frequLB<-cbind(dataLB,adddataLB)

names(frequLB)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3LB<-frequLB[count.fields(textConnection(frequLB$sampleNames),sep=",")>=5,]

require(xlsx)
write.xlsx2(frequ3LB,"CnvregionswithmeanhalleriLangBest_minRead2_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE)
frequbigLB<-read.xlsx2("CnvregionswithmeanhalleriLangBest_minRead2_WL500_minL2_min5samples.xlsx",1,header=TRUE)
frequbigLB[,2]<-as.numeric(as.character(frequbigLB[,2]))
frequbigLB[,3]<-as.numeric(as.character(frequbigLB[,3]))

###################################################
#find overlaps between windows and A. halleri genes#
###################################################
frequbigLB_GRange<-GRanges(seqnames=tolower(frequbigLB$Scaffold),ranges=IRanges(start=frequbigLB$start,end=frequbigLB$end))
values(frequbigLB_GRange)<-frequbigLB[,6:11]

testLB<-findOverlaps(frequbigLB_GRange,HalGenes)
test2LB<-pintersect(frequbigLB_GRange[queryHits(testLB)],HalGenes[subjectHits(testLB)])
test3LB<-as.data.frame(test2LB)
frequbigLB_halgenes=mergeByOverlaps(frequbigLB_GRange,HalGenes,type=c("any"))
frequbigLB_halgenes_df=data.frame(as.data.frame(frequbigLB_halgenes$frequbigLB_GRange),as.data.frame(frequbigLB_halgenes$HalGenes),test3LB$width)
#frequbigLB_halgenes_df<-frequbigLB_halgenes_df[,-c(6,7)]
colnames(frequbigLB_halgenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigLB_halgenes_df<-frequbigLB_halgenes_df[(frequbigLB_halgenes_df$Widthofoverlap>=(0.9*frequbigLB_halgenes_df$gene_size)),]
write.xlsx2(frequbigLB_halgenes_df,"LB_temp.xlsx",col.names=TRUE,row.names=FALSE)



regionsNP<-cnvr(resCNMOPSNP)
singleNP<-cnvs(resCNMOPSNP)
testNP<-findOverlaps(regionsNP,singleNP)
test2NP<- DataFrame(splitAsList(singleNP$sampleName[subjectHits(testNP)], queryHits(testNP)),splitAsList(singleNP$median[subjectHits(testNP)], queryHits(testNP)),splitAsList(singleNP$mean[subjectHits(testNP)], queryHits(testNP)),splitAsList(singleNP$CN[subjectHits(testNP)], queryHits(testNP)))
colnames(test2NP)<-c("sampleNames","mean","median","CN")
test7NP<-apply(as.data.frame(test2NP),1,dupl_cov)
mcols(regionsNP)<-as.data.frame(as.matrix(test7NP))

dataNP<-as.data.frame(regionsNP)

adddataNP<-do.call("rbind",apply(dataNP,1,data_change))
dataNP<-dataNP[,-6]
frequNP<-cbind(dataNP,adddataNP)

names(frequNP)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3NP<-frequNP[count.fields(textConnection(frequNP$sampleNames),sep=",")>=5,]

require(xlsx)
write.xlsx2(frequ3NP,"CnvregionswithmeanhalleriNossPais_minRead2_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE)
frequbigNP<-read.xlsx2("CnvregionswithmeanhalleriNossPais_minRead2_WL500_minL2_min5samples.xlsx",1,header=TRUE)
frequbigNP[,2]<-as.numeric(as.character(frequbigNP[,2]))
frequbigNP[,3]<-as.numeric(as.character(frequbigNP[,3]))

###################################################
#find overlaps between windows and A. halleri genes#
###################################################
frequbigNP_GRange<-GRanges(seqnames=tolower(frequbigNP$Scaffold),ranges=IRanges(start=frequbigNP$start,end=frequbigNP$end))
values(frequbigNP_GRange)<-frequbigNP[,6:11]

testNP<-findOverlaps(frequbigNP_GRange,HalGenes)
test2NP<-pintersect(frequbigNP_GRange[queryHits(testNP)],HalGenes[subjectHits(testNP)])
test3NP<-as.data.frame(test2NP)
frequbigNP_halgenes=mergeByOverlaps(frequbigNP_GRange,HalGenes,type=c("any"))
frequbigNP_halgenes_df=data.frame(as.data.frame(frequbigNP_halgenes$frequbigNP_GRange),as.data.frame(frequbigNP_halgenes$HalGenes),test3NP$width)
#frequbigNP_halgenes_df<-frequbigNP_halgenes_df[,-c(6,7)]
colnames(frequbigNP_halgenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigNP_halgenes_df<-frequbigNP_halgenes_df[(frequbigNP_halgenes_df$Widthofoverlap>=(0.9*frequbigNP_halgenes_df$gene_size)),]
write.xlsx2(frequbigNP_halgenes_df,"NP_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigNP_halgenes_df<-read.xlsx2("NP_temp.xlsx",1,header=T)
frequbigKK_halgenes_df<-read.xlsx2("KK_temp.xlsx",1,header=T)
frequbigLB_halgenes_df<-read.xlsx2("LB_temp.xlsx",1,header=T)

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

frequbigKK_halgenes_df_w_orthogroup= merge(OG_Ahalleri,frequbigKK_halgenes_df,by.y="Gene", by.x="Ah_ID",all.y=TRUE)
frequbigKK_halgenes_df_w_MM=merge(frequbigKK_halgenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
frequbigKK_halgenes_df_w_MM_lists=merge(frequbigKK_halgenes_df_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
frequbigKK_halgenes_df_w_MM_lists<-frequbigKK_halgenes_df_w_MM_lists[,c(2,10:26,1,3:9,27:70)]
colnames(frequbigKK_halgenes_df_w_MM_lists)[27]<-"Mapman_category"

frequbigLB_halgenes_df_w_orthogroup= merge(OG_Ahalleri,frequbigLB_halgenes_df,by.y="Gene", by.x="Ah_ID",all.y=TRUE)
frequbigLB_halgenes_df_w_MM=merge(frequbigLB_halgenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
frequbigLB_halgenes_df_w_MM_lists=merge(frequbigLB_halgenes_df_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
frequbigLB_halgenes_df_w_MM_lists<-frequbigLB_halgenes_df_w_MM_lists[,c(2,10:26,1,3:9,27:70)]
colnames(frequbigLB_halgenes_df_w_MM_lists)[27]<-"Mapman_category"

frequbigNP_halgenes_df_w_orthogroup= merge(OG_Ahalleri,frequbigNP_halgenes_df,by.y="Gene", by.x="Ah_ID",all.y=TRUE)
frequbigNP_halgenes_df_w_MM=merge(frequbigNP_halgenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
frequbigNP_halgenes_df_w_MM_lists=merge(frequbigNP_halgenes_df_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
frequbigNP_halgenes_df_w_MM_lists<-frequbigNP_halgenes_df_w_MM_lists[,c(2,10:26,1,3:9,27:70)]
colnames(frequbigNP_halgenes_df_w_MM_lists)[27]<-"Mapman_category"

options(java.parameters = "-Xmx12000m")

write.xlsx2(frequbigKK_halgenes_df_w_MM_lists,"GeneshalleriKromKosicnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",sheetName="GeneshalleriKromKosicnmops",col.names=TRUE,row.names=FALSE,quote=F)
write.xlsx2(frequbigLB_halgenes_df_w_MM_lists,"GeneshalleriLangBestcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",sheetName="GeneshalleriLangBestcnmops",col.names=TRUE,row.names=FALSE,quote=F)
write.xlsx2(frequbigNP_halgenes_df_w_MM_lists,"GeneshalleriNossPaiscnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",sheetName="GeneshalleriNossPaiscnmops",col.names=TRUE,row.names=FALSE,quote=F)

write.table(frequbigKK_halgenes_df_w_MM_lists,"GeneshalleriKromKosicnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.table",row.names=FALSE,quote=F)
write.table(frequbigLB_halgenes_df_w_MM_lists,"GeneshalleriLangBestcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.table",row.names=FALSE,quote=F)
write.table(frequbigNP_halgenes_df_w_MM_lists,"GeneshalleriNossPaiscnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.table",row.names=FALSE,quote=F)


#################################################
#Pais to Noss#
#################################################

regionsPN<-cnvr(resCNMOPSPN)
singlePN<-cnvs(resCNMOPSPN)
testPN<-findOverlaps(regionsPN,singlePN)
test2PN<- DataFrame(splitAsList(singlePN$sampleName[subjectHits(testPN)], queryHits(testPN)),splitAsList(singlePN$median[subjectHits(testPN)], queryHits(testPN)),splitAsList(singlePN$mean[subjectHits(testPN)], queryHits(testPN)),splitAsList(singlePN$CN[subjectHits(testPN)], queryHits(testPN)))
colnames(test2PN)<-c("sampleNames","mean","median","CN")
test7PN<-apply(as.data.frame(test2PN),1,dupl_cov)
mcols(regionsPN)<-as.data.frame(as.matrix(test7PN))

dataPN<-as.data.frame(regionsPN)

adddataPN<-do.call("rbind",apply(dataPN,1,data_change))
dataPN<-dataPN[,-6]
frequPN<-cbind(dataPN,adddataPN)

names(frequPN)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3PN<-frequPN[count.fields(textConnection(frequPN$sampleNames),sep=",")>=5,]

require(xlsx)
write.xlsx2(frequ3PN,"CnvregionswithmeanhalleriPaisNoss_minRead2_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE)
frequbigPN<-read.xlsx2("CnvregionswithmeanhalleriPaisNoss_minRead2_WL500_minL2_min5samples.xlsx",1,header=TRUE)
frequbigPN[,2]<-as.numeric(as.character(frequbigPN[,2]))
frequbigPN[,3]<-as.numeric(as.character(frequbigPN[,3]))

###################################################
#find overlaps between windows and A. halleri genes#
###################################################
frequbigPN_GRange<-GRanges(seqnames=tolower(frequbigPN$Scaffold),ranges=IRanges(start=frequbigPN$start,end=frequbigPN$end))
values(frequbigPN_GRange)<-frequbigPN[,6:11]

testPN<-findOverlaps(frequbigPN_GRange,HalGenes)
test2PN<-pintersect(frequbigPN_GRange[queryHits(testPN)],HalGenes[subjectHits(testPN)])
test3PN<-as.data.frame(test2PN)
frequbigPN_halgenes=mergeByOverlaps(frequbigPN_GRange,HalGenes,type=c("any"))
frequbigPN_halgenes_df=data.frame(as.data.frame(frequbigPN_halgenes$frequbigPN_GRange),as.data.frame(frequbigPN_halgenes$HalGenes),test3PN$width)
#frequbigPN_halgenes_df<-frequbigPN_halgenes_df[,-c(6,7)]
colnames(frequbigPN_halgenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigPN_halgenes_df<-frequbigPN_halgenes_df[(frequbigPN_halgenes_df$Widthofoverlap>=(0.9*frequbigPN_halgenes_df$gene_size)),]
write.xlsx2(frequbigPN_halgenes_df,"PN_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigPN_halgenes_df<-read.xlsx2("PN_temp.xlsx",1,header=TRUE)

frequbigPN_halgenes_df_w_orthogroup= merge(OG_Ahalleri,frequbigPN_halgenes_df,by.y="Gene", by.x="Ah_ID",all.y=TRUE)
frequbigPN_halgenes_df_w_MM=merge(frequbigPN_halgenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
frequbigPN_halgenes_df_w_MM_lists=merge(frequbigPN_halgenes_df_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
frequbigPN_halgenes_df_w_MM_lists<-frequbigPN_halgenes_df_w_MM_lists[,c(2,10:26,1,3:9,27:70)]
colnames(frequbigPN_halgenes_df_w_MM_lists)[27]<-"Mapman_category"

write.xlsx2(frequbigPN_halgenes_df_w_MM_lists,"GeneshalleriPaisNosscnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",sheetName="GeneshalleriPaisNosscnmops",col.names=TRUE,row.names=FALSE,quote=F)
write.table(frequbigPN_halgenes_df_w_MM_lists,"GeneshalleriPaisNosscnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.table",row.names=FALSE,quote=F)


#################################################
#Best to Lang#
#################################################

regionsBL<-cnvr(resCNMOPSBL)
singleBL<-cnvs(resCNMOPSBL)
testBL<-findOverlaps(regionsBL,singleBL)
test2BL<- DataFrame(splitAsList(singleBL$sampleName[subjectHits(testBL)], queryHits(testBL)),splitAsList(singleBL$median[subjectHits(testBL)], queryHits(testBL)),splitAsList(singleBL$mean[subjectHits(testBL)], queryHits(testBL)),splitAsList(singleBL$CN[subjectHits(testBL)], queryHits(testBL)))
colnames(test2BL)<-c("sampleNames","mean","median","CN")
test7BL<-apply(as.data.frame(test2BL),1,dupl_cov)
mcols(regionsBL)<-as.data.frame(as.matrix(test7BL))

dataBL<-as.data.frame(regionsBL)

adddataBL<-do.call("rbind",apply(dataBL,1,data_change))
dataBL<-dataBL[,-6]
frequBL<-cbind(dataBL,adddataBL)

names(frequBL)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3BL<-frequBL[count.fields(textConnection(frequBL$sampleNames),sep=",")>=5,]

require(xlsx)
write.xlsx2(frequ3BL,"CnvregionswithmeanhalleriBestLang_minRead2_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE)
frequbigBL<-read.xlsx2("CnvregionswithmeanhalleriBestLang_minRead2_WL500_minL2_min5samples.xlsx",1,header=TRUE)
frequbigBL[,2]<-as.numeric(as.character(frequbigBL[,2]))
frequbigBL[,3]<-as.numeric(as.character(frequbigBL[,3]))

###################################################
#find overlaps between windows and A. halleri genes#
###################################################
frequbigBL_GRange<-GRanges(seqnames=tolower(frequbigBL$Scaffold),ranges=IRanges(start=frequbigBL$start,end=frequbigBL$end))
values(frequbigBL_GRange)<-frequbigBL[,6:11]

testBL<-findOverlaps(frequbigBL_GRange,HalGenes)
test2BL<-pintersect(frequbigBL_GRange[queryHits(testBL)],HalGenes[subjectHits(testBL)])
test3BL<-as.data.frame(test2BL)
frequbigBL_halgenes=mergeByOverlaps(frequbigBL_GRange,HalGenes,type=c("any"))
frequbigBL_halgenes_df=data.frame(as.data.frame(frequbigBL_halgenes$frequbigBL_GRange),as.data.frame(frequbigBL_halgenes$HalGenes),test3BL$width)
#frequbigBL_halgenes_df<-frequbigBL_halgenes_df[,-c(6,7)]
colnames(frequbigBL_halgenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigBL_halgenes_df<-frequbigBL_halgenes_df[(frequbigBL_halgenes_df$Widthofoverlap>=(0.9*frequbigBL_halgenes_df$gene_size)),]
write.xlsx2(frequbigBL_halgenes_df,"BL_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigBL_halgenes_df<-read.xlsx2("BL_temp.xlsx",1,header=TRUE)

frequbigBL_halgenes_df_w_orthogroup= merge(OG_Ahalleri,frequbigBL_halgenes_df,by.y="Gene", by.x="Ah_ID",all.y=TRUE)
frequbigBL_halgenes_df_w_MM=merge(frequbigBL_halgenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
frequbigBL_halgenes_df_w_MM_lists=merge(frequbigBL_halgenes_df_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
frequbigBL_halgenes_df_w_MM_lists<-frequbigBL_halgenes_df_w_MM_lists[,c(2,10:26,1,3:9,27:70)]
colnames(frequbigBL_halgenes_df_w_MM_lists)[27]<-"Mapman_category"

write.xlsx2(frequbigBL_halgenes_df_w_MM_lists,"GeneshalleriBestLangcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",sheetName="GeneshalleriBestLangcnmops",col.names=TRUE,row.names=FALSE,quote=F)
write.table(frequbigBL_halgenes_df_w_MM_lists,"GeneshalleriBestLangcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.table",row.names=FALSE,quote=F)


#################################################
#Kosi to Krom#
#################################################

regionsKoKr<-cnvr(resCNMOPSKoKr)
singleKoKr<-cnvs(resCNMOPSKoKr)
testKoKr<-findOverlaps(regionsKoKr,singleKoKr)
test2KoKr<- DataFrame(splitAsList(singleKoKr$sampleName[subjectHits(testKoKr)], queryHits(testKoKr)),splitAsList(singleKoKr$median[subjectHits(testKoKr)], queryHits(testKoKr)),splitAsList(singleKoKr$mean[subjectHits(testKoKr)], queryHits(testKoKr)),splitAsList(singleKoKr$CN[subjectHits(testKoKr)], queryHits(testKoKr)))
colnames(test2KoKr)<-c("sampleNames","mean","median","CN")
test7KoKr<-apply(as.data.frame(test2KoKr),1,dupl_cov)
mcols(regionsKoKr)<-as.data.frame(as.matrix(test7KoKr))

dataKoKr<-as.data.frame(regionsKoKr)

adddataKoKr<-do.call("rbind",apply(dataKoKr,1,data_change))
dataKoKr<-dataKoKr[,-6]
frequKoKr<-cbind(dataKoKr,adddataKoKr)

names(frequKoKr)<-c("Scaffold","start","end","width","strand","sampleNames","mean","median","CNsingle","CN_class","Percentagesupport")

frequ3KoKr<-frequKoKr[count.fields(textConnection(frequKoKr$sampleNames),sep=",")>=5,]

require(xlsx)
write.xlsx2(frequ3KoKr,"CnvregionswithmeanhalleriKosiKrom_minRead2_WL500_minL2_min5samples.xlsx",col.names=TRUE,row.names=FALSE)
frequbigKoKr<-read.xlsx2("CnvregionswithmeanhalleriKosiKrom_minRead2_WL500_minL2_min5samples.xlsx",1,header=TRUE)
frequbigKoKr[,2]<-as.numeric(as.character(frequbigKoKr[,2]))
frequbigKoKr[,3]<-as.numeric(as.character(frequbigKoKr[,3]))

###################################################
#find overlaps between windows and A. halleri genes#
###################################################
frequbigKoKr_GRange<-GRanges(seqnames=tolower(frequbigKoKr$Scaffold),ranges=IRanges(start=frequbigKoKr$start,end=frequbigKoKr$end))
values(frequbigKoKr_GRange)<-frequbigKoKr[,6:11]

testKoKr<-findOverlaps(frequbigKoKr_GRange,HalGenes)
test2KoKr<-pintersect(frequbigKoKr_GRange[queryHits(testKoKr)],HalGenes[subjectHits(testKoKr)])
test3KoKr<-as.data.frame(test2KoKr)
frequbigKoKr_halgenes=mergeByOverlaps(frequbigKoKr_GRange,HalGenes,type=c("any"))
frequbigKoKr_halgenes_df=data.frame(as.data.frame(frequbigKoKr_halgenes$frequbigKoKr_GRange),as.data.frame(frequbigKoKr_halgenes$HalGenes),test3KoKr$width)
#frequbigKoKr_halgenes_df<-frequbigKoKr_halgenes_df[,-c(6,7)]
colnames(frequbigKoKr_halgenes_df)=c("Chr", "Copy_start", "Copy_end", "Copy_width", "Copy_strand", "sampleNames","meanCN","medianCN","CNsingle","CN_class","Percentagesupport","Chr_gene", "gene_start", "gene_end", "gene_size", "gene_strand", "Gene","Widthofoverlap" )

frequbigKoKr_halgenes_df<-frequbigKoKr_halgenes_df[(frequbigKoKr_halgenes_df$Widthofoverlap>=(0.9*frequbigKoKr_halgenes_df$gene_size)),]
write.xlsx2(frequbigKoKr_halgenes_df,"KoKr_temp.xlsx",col.names=TRUE,row.names=FALSE)

frequbigKoKr_halgenes_df<-read.xlsx2("KoKr_temp.xlsx",1,header=TRUE)

frequbigKoKr_halgenes_df_w_orthogroup= merge(OG_Ahalleri,frequbigKoKr_halgenes_df,by.y="Gene", by.x="Ah_ID",all.y=TRUE)
frequbigKoKr_halgenes_df_w_MM=merge(frequbigKoKr_halgenes_df_w_orthogroup,Mapman_df,by.x="Ath_ID",by.y="IDENTIFIER",all.x=TRUE)
frequbigKoKr_halgenes_df_w_MM_lists=merge(frequbigKoKr_halgenes_df_w_MM,All_lists,by.x="Ath_ID",by.y="AGI_Code",all.x=T)
frequbigKoKr_halgenes_df_w_MM_lists<-frequbigKoKr_halgenes_df_w_MM_lists[,c(2,10:26,1,3:9,27:70)]
colnames(frequbigKoKr_halgenes_df_w_MM_lists)[27]<-"Mapman_category"

write.xlsx2(frequbigKoKr_halgenes_df_w_MM_lists,"GeneshalleriKosiKromcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",sheetName="GeneshalleriKosiKromcnmops",col.names=TRUE,row.names=FALSE,quote=F)
write.table(frequbigKoKr_halgenes_df_w_MM_lists,"GeneshalleriKosiKromcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.table",row.names=FALSE,quote=F)



options(java.parameters = "-Xmx12000m")
require("xlsx")

Kosi<-read.xlsx2("GeneshalleriKosiKromcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",1)
Krom<-read.xlsx2("GeneshalleriKromKosicnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",1)

Krom_unique<-Krom[!duplicated(Krom$Ath_ID,Krom$Ah_ID,Krom$Copy_start,Krom$Copy_end),]
Kosi_unique<-Kosi[!duplicated(Kosi$Ath_ID,Kosi$Ah_ID,Kosi$Copy_start,Kosi$Copy_end),]

Krom_unique$Ah_ID<-as.character(Krom_unique$Ah_ID)
Kosi_unique$Ah_ID<-as.character(Kosi_unique$Ah_ID)
Krom_unique$CN_class<-as.character(Krom_unique$CN_class)
Kosi_unique$CN_class<-as.character(Kosi_unique$CN_class)

Kromonly<-Krom_unique[!((Krom_unique$Ah_ID%in%Kosi_unique$Ah_ID)&(Krom_unique$CN_class==Kosi_unique$CN_class[match(Krom_unique$Ah_ID,Kosi_unique$Ah_ID)])),]
Kosionly<-Kosi_unique[!((Kosi_unique$Ah_ID%in%Krom_unique$Ah_ID)&(Kosi_unique$CN_class==Krom_unique$CN_class[match(Kosi_unique$Ah_ID,Krom_unique$Ah_ID)])),]

write.table(Kromonly,"Krom_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(Kosionly,"Kosi_cnvs_only_strict.table",row.names=F,sep="\t")

Best<-read.xlsx2("GeneshalleriBestLangcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",1)
Lang<-read.xlsx2("GeneshalleriLangBestcnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",1)

Lang_unique<-Lang[!duplicated(Lang$Ath_ID,Lang$Ah_ID,Lang$Copy_start,Lang$Copy_end),]
Best_unique<-Best[!duplicated(Best$Ath_ID,Best$Ah_ID,Best$Copy_start,Best$Copy_end),]

Lang_unique$Ah_ID<-as.character(Lang_unique$Ah_ID)
Best_unique$Ah_ID<-as.character(Best_unique$Ah_ID)
Lang_unique$CN_class<-as.character(Lang_unique$CN_class)
Best_unique$CN_class<-as.character(Best_unique$CN_class)

Langonly<-Lang_unique[!((Lang_unique$Ah_ID%in%Best_unique$Ah_ID)&(Lang_unique$CN_class==Best_unique$CN_class[match(Lang_unique$Ah_ID,Best_unique$Ah_ID)])),]
Bestonly<-Best_unique[!((Best_unique$Ah_ID%in%Lang_unique$Ah_ID)&(Best_unique$CN_class==Lang_unique$CN_class[match(Best_unique$Ah_ID,Lang_unique$Ah_ID)])),]

write.table(Langonly,"Lang_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(Bestonly,"Best_cnvs_only_strict.table",row.names=F,sep="\t")

Pais<-read.xlsx2("GeneshalleriPaisNosscnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",1)
Noss<-read.xlsx2("GeneshalleriNossPaiscnvregions_minRead2_WL500_minL2_minoverlap90percentofgene_atleast5samples.xlsx",1)

Noss_unique<-Noss[!duplicated(Noss$Ath_ID,Noss$Ah_ID,Noss$Copy_start,Noss$Copy_end),]
Pais_unique<-Pais[!duplicated(Pais$Ath_ID,Pais$Ah_ID,Pais$Copy_start,Pais$Copy_end),]

Noss_unique$Ah_ID<-as.character(Noss_unique$Ah_ID)
Pais_unique$Ah_ID<-as.character(Pais_unique$Ah_ID)
Noss_unique$CN_class<-as.character(Noss_unique$CN_class)
Pais_unique$CN_class<-as.character(Pais_unique$CN_class)

Nossonly<-Noss_unique[!((Noss_unique$Ah_ID%in%Pais_unique$Ah_ID)&(Noss_unique$CN_class==Pais_unique$CN_class[match(Noss_unique$Ah_ID,Pais_unique$Ah_ID)])),]
Paisonly<-Pais_unique[!((Pais_unique$Ah_ID%in%Noss_unique$Ah_ID)&(Pais_unique$CN_class==Noss_unique$CN_class[match(Pais_unique$Ah_ID,Noss_unique$Ah_ID)])),]

write.table(Nossonly,"Noss_cnvs_only_strict.table",row.names=F,sep="\t")
write.table(Paisonly,"Pais_cnvs_only_strict.table",row.names=F,sep="\t")




