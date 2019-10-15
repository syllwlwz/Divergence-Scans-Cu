options(java.parameters = "-Xmx12000m")
require(xlsx)

KK_classical<-read.xlsx("Genes_KromKosi_wVarLD.xlsx",1,header=T)
KK_high<-read.xlsx("Genes_KromKosi_wVarLD.xlsx",2,header=T)
KK_Indels<-read.xlsx("Genes_KromKosi_wVarLD.xlsx",3,header=T)
KK_CNVs<-read.xlsx("Genes_KromKosi_wVarLD.xlsx",4,header=T)

#fix wrong gene start and end and strand
genelist<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(genelist)<-c("Scaffold","Source","Type","Start_pos","End_pos","Score","Strand","Phase","Attributes","Ah_ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)
genelist<-genelist[!duplicated(genelist[,-9]),]

KK_classical$gene_start<-genelist$Start_pos[match(KK_classical$Ah_ID,genelist$Ah_ID)]
KK_classical$gene_end<-genelist$End_pos[match(KK_classical$Ah_ID,genelist$Ah_ID)]
KK_classical$gene_strand<-genelist$Strand[match(KK_classical$Ah_ID,genelist$Ah_ID)]


dup<-duplicated(KK_classical$Ah_ID)
KK_classical$Rank_in_outlier_metric<-as.numeric(as.character(KK_classical$Rank_in_outlier_metric))
duplist<-unique(KK_classical$Ah_ID[dup])
KK_classical_uniquewindows<-KK_classical[!KK_classical$Ah_ID%in%duplist,]
for (i in 1:length(duplist))
	{Ahgene<-KK_classical[KK_classical$Ah_ID==duplist[i],]
	KK_classical_uniquewindows<-rbind(KK_classical_uniquewindows,Ahgene[Ahgene$Rank_in_outlier_metric==min(Ahgene$Rank_in_outlier_metric),])
	}

str(KK_classical_uniquewindows)
#Mapman duplicates + same rank

#remove Mapman duplicates and same rank duplicates
KK_classical_uniquewindows$Mapman_category<-as.character(KK_classical_uniquewindows$Mapman_category)

KK_classical_uniquewindows2<-data.frame()
for (i in unique(KK_classical_uniquewindows$Ah_ID[duplicated(KK_classical_uniquewindows$Ah_ID)]))
	{KK_classical_uniquewindows2<-rbind(KK_classical_uniquewindows2,KK_classical_uniquewindows[KK_classical_uniquewindows$Ah_ID==i,][1,])
	KK_classical_uniquewindows2$Mapman_category[nrow(KK_classical_uniquewindows2)]<-paste(unique(KK_classical_uniquewindows$Mapman_category[KK_classical_uniquewindows$Ah_ID==i]),sep="/",collapse="/")
	}

KK_GS<-rbind(KK_classical_uniquewindows[!KK_classical_uniquewindows$Ah_ID%in%KK_classical_uniquewindows$Ah_ID[duplicated(KK_classical_uniquewindows$Ah_ID)],],KK_classical_uniquewindows2)




#remove Mapman duplicates for indels
KK_Indels$Mapman_category<-as.character(KK_Indels$Mapman_category)

KK_Indels2<-data.frame()
for (i in unique(KK_Indels$Ah_ID[duplicated(KK_Indels$Ah_ID)]))
	{KK_Indels2<-rbind(KK_Indels2,KK_Indels[KK_Indels$Ah_ID==i,])
	KK_Indels2$Mapman_category[KK_Indels2$Ah_ID==i]<-paste(unique(KK_Indels$Mapman_category[KK_Indels$Ah_ID==i]),sep="/",collapse="/")
	}

KK_INDEL<-rbind(KK_Indels[!KK_Indels$Ah_ID%in%KK_Indels$Ah_ID[duplicated(KK_Indels$Ah_ID)],],KK_Indels2)
KK_INDEL<-KK_INDEL[!duplicated(KK_INDEL),]
KK_INDEL<-na.omit(KK_INDEL)

#remove Mapman duplicates for high effect SNPs
KK_high$Mapman_category<-as.character(KK_high$Mapman_category)

KK_high2<-data.frame()
for (i in unique(KK_high$Ah_ID[duplicated(KK_high$Ah_ID)]))
	{KK_high2<-rbind(KK_high2,KK_high[KK_high$Ah_ID==i,])
	KK_high2$Mapman_category[KK_high2$Ah_ID==i]<-paste(unique(KK_high$Mapman_category[KK_high$Ah_ID==i]),sep="/",collapse="/")
	}

KK_HIGH<-rbind(KK_high[!KK_high$Ah_ID%in%KK_high$Ah_ID[duplicated(KK_high$Ah_ID)],],KK_high2)
KK_HIGH<-KK_HIGH[!duplicated(KK_HIGH),]
KK_HIGH<-na.omit(KK_HIGH)

#remove Mapman duplicates for high effect SNPs
KK_CNVs$Mapman_category<-as.character(KK_CNVs$Mapman_category)

KK_CNVs2<-data.frame()
for (i in unique(KK_CNVs$Ah_ID[duplicated(KK_CNVs$Ah_ID)]))
	{KK_CNVs2<-rbind(KK_CNVs2,KK_CNVs[KK_CNVs$Ah_ID==i,])
	KK_CNVs2$Mapman_category[KK_CNVs2$Ah_ID==i]<-paste(unique(KK_CNVs$Mapman_category[KK_CNVs$Ah_ID==i]),sep="/",collapse="/")
	}

KK_CNV<-rbind(KK_CNVs[!KK_CNVs$Ah_ID%in%KK_CNVs$Ah_ID[duplicated(KK_CNVs$Ah_ID)],],KK_CNVs2)
KK_CNV<-KK_CNV[!duplicated(KK_CNV),]
KK_CNV<-na.omit(KK_CNV)

write.xlsx(na.omit(KK_GS),"Genes_KromKosi_dedup.xlsx",sheetName="Genome_scans",col.names=TRUE,row.names=FALSE)
write.xlsx(KK_INDEL,"Genes_KromKosi_dedup.xlsx",sheetName="High_effect_indels",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx(KK_HIGH,"Genes_KromKosi_dedup.xlsx",sheetName="High_effect_SNPs",col.names=TRUE,row.names=FALSE,append=T)
write.xlsx(KK_CNV,"Genes_KromKosi_dedup.xlsx",sheetName="Copy_number_variations",col.names=TRUE,row.names=FALSE,append=T)



length(unique(KK_CNV$Ah_ID[grep("Kosi",KK_CNV$sampleNames)]))
#261
length(unique(KK_CNV$Ah_ID[grep("Krom",KK_CNV$sampleNames)]))
#119


Kosi_CNVs<-KK_CNV[grep("Kosi",KK_CNV$sampleNames),]
length(unique(Kosi_CNVs$Ah_ID[Kosi_CNVs$ALT=="<DEL>"]))
#247
length(unique(Kosi_CNVs$Ah_ID[Kosi_CNVs$ALT=="<DUP>"]))
#14

Krom_CNVs<-KK_CNV[grep("Krom",KK_CNV$sampleNames),]
length(unique(Krom_CNVs$Ah_ID[Krom_CNVs$ALT=="<DEL>"]))
#110
length(unique(Krom_CNVs$Ah_ID[Krom_CNVs$ALT=="<DUP>"]))
#9


























