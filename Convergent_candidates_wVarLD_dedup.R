options(java.parameters = "-Xmx12000m")
require(xlsx)


KK_GS<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Krom-Kosi//Genes_KromKosi_wVarLD.xlsx",1,header=T)
KK_HigheffectSNPs<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Krom-Kosi//Genes_KromKosi_wVarLD.xlsx",2,header=T)
KK_HigheffectIndels<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Krom-Kosi//Genes_KromKosi_wVarLD.xlsx",3,header=T)

LB_GS<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Lang-Best//Genes_LangBest_wVarLD.xlsx",1,header=T)
LB_HigheffectSNPs<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Lang-Best//Genes_LangBest_wVarLD.xlsx",2,header=T)
LB_HigheffectIndels<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Lang-Best//Genes_LangBest_wVarLD.xlsx",3,header=T)

NP_GS<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//Genes_NossPais_wVarLD.xlsx",1,header=T)
NP_HigheffectSNPs<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//Genes_NossPais_wVarLD.xlsx",2,header=T)
NP_HigheffectIndels<-read.xlsx2("F://Cu_project//Genome_scans//Genom_scans_Cu//Noss-Pais//Genes_NossPais_wVarLD.xlsx",3,header=T)

KK_LB_GS<-merge(KK_GS,LB_GS,by="Ahal_Gene")
KK_LB_GS<-KK_LB_GS[!duplicated(KK_LB_GS),]

KK_NP_GS<-merge(KK_GS,NP_GS,by="Ahal_Gene")
KK_NP_GS<-KK_NP_GS[!duplicated(KK_NP_GS),]

LB_NP_GS<-merge(LB_GS,NP_GS,by="Ahal_Gene")
LB_NP_GS<-LB_NP_GS[!duplicated(LB_NP_GS),]

GS_only<-rbind(KK_LB_GS[,c(1:7,26:30)],KK_NP_GS[,c(1:7,26:30)],LB_NP_GS[,c(1:7,26:30)])
GS_only<-GS_only[!duplicated(GS_only),]

write.xlsx2(GS_only,"Convergent_candidates_wVarLD.xlsx",sheetName="Genome_scans",row.names=F)

colnames(KK_HigheffectIndels)[c(1:7,22:24,26)]<-colnames(KK_GS)[c(1:7,26:28,30)]
colnames(KK_HigheffectSNPs)[c(1:7,13,16:17,19)]<-colnames(KK_GS)[c(1:7,26:28,30)]
KK_GS_HEvar<-rbind(KK_GS[,c(1:7,26:28,30)],KK_HigheffectSNPs[,c(1:7,13,16:17,19)],KK_HigheffectIndels[,c(1:7,22:24,26)])

colnames(LB_HigheffectIndels)[c(1:7,22:24,26)]<-colnames(LB_GS)[c(1:7,26:28,30)]
colnames(LB_HigheffectSNPs)[c(1:7,13,16:17,19)]<-colnames(LB_GS)[c(1:7,26:28,30)]
LB_GS_HEvar<-rbind(LB_GS[,c(1:7,26:28,30)],LB_HigheffectSNPs[,c(1:7,13,16:17,19)],LB_HigheffectIndels[,c(1:7,22:24,26)])

colnames(NP_HigheffectIndels)[c(1:7,22:24,26)]<-colnames(NP_GS)[c(1:7,26:28,30)]
colnames(NP_HigheffectSNPs)[c(1:7,13,16:17,19)]<-colnames(NP_GS)[c(1:7,26:28,30)]
NP_GS_HEvar<-rbind(NP_GS[,c(1:7,26:28,30)],NP_HigheffectSNPs[,c(1:7,13,16:17,19)],NP_HigheffectIndels[,c(1:7,22:24,26)])

KK_LB_GS_HEvar<-merge(KK_GS_HEvar,LB_GS_HEvar,by="Ahal_Gene")
KK_LB_GS_HEvar<-KK_LB_GS_HEvar[!duplicated(KK_LB_GS_HEvar),]

KK_NP_GS_HEvar<-merge(KK_GS_HEvar,NP_GS_HEvar,by="Ahal_Gene")
KK_NP_GS_HEvar<-KK_NP_GS_HEvar[!duplicated(KK_NP_GS_HEvar),]

LB_NP_GS_HEvar<-merge(LB_GS_HEvar,NP_GS_HEvar,by="Ahal_Gene")
LB_NP_GS_HEvar<-LB_NP_GS_HEvar[!duplicated(LB_NP_GS_HEvar),]

GS_HEvar<-rbind(KK_LB_GS_HEvar,KK_NP_GS_HEvar,LB_NP_GS_HEvar)
GS_HEvar<-GS_HEvar[!duplicated(GS_HEvar),]

write.xlsx2(GS_HEvar[,c(1:11)],"Convergent_candidates_wVarLD.xlsx",sheetName="Genome_scans_higheffectvar",row.names=F,append=T)

#CNVs
KK_CNVs<-read.xlsx2("F://Cu_project//CNVS_Indels//CNVs_overlap_strict.xlsx",1,header=T)
LB_CNVs<-read.xlsx2("F://Cu_project//CNVS_Indels//CNVs_overlap_strict.xlsx",3,header=T)
NP_CNVs<-read.xlsx2("F://Cu_project//CNVS_Indels//CNVs_overlap_strict.xlsx",2,header=T)

KK_LB_CNVs<-merge(KK_CNVs,LB_CNVs,by="Ahal_Gene")
KK_LB_CNVs<-KK_LB_CNVs[!duplicated(KK_LB_CNVs),]

KK_NP_CNVs<-merge(KK_CNVs,NP_CNVs,by="Ahal_Gene")
KK_NP_CNVs<-KK_NP_CNVs[!duplicated(KK_NP_CNVs),]

LB_NP_CNVs<-merge(LB_CNVs,NP_CNVs,by="Ahal_Gene")
LB_NP_CNVs<-LB_NP_CNVs[!duplicated(LB_NP_CNVs),]

CNVs<-rbind(KK_LB_CNVs,KK_NP_CNVs,LB_NP_CNVs)
CNVs<-CNVs[!duplicated(CNVs),]

write.xlsx2(CNVs,"Convergent_candidates_wVarLD.xlsx",sheetName="Copy_number_variants",row.names=F,append=T)

Gene<-read.xlsx2("Convergent_candidates_wVarLD.xlsx",1,header=T)
Geneplus<-read.xlsx2("Convergent_candidates_wVarLD.xlsx",2,header=T)
Geneplusonly<-Geneplus[!Geneplus$Ahal_Gene%in%Gene$Ahal_Gene,]
write.table(Geneplusonly,"Convergent_genes_GSplusHEvaronly_wVarLD.table",sep="\t",row.names=F)

######################################################
#MapMan#
######################################################
GS_only_Mapman<-GS_only[,c(2,9)]
GS_only_Mapman<-GS_only_Mapman[!GS_only_Mapman[,1]=="",]
GS_only_Mapman<-GS_only_Mapman[!duplicated(GS_only_Mapman$Gene.x),]

write.table(GS_only_Mapman,"GS_only_Mapman_wVarLD.txt",sep="\t",row.names=F,col.names=F,quote=F)

GS_HEvar_Mapman<-GS_HEvar[,c(2,9)]
GS_HEvar_Mapman<-GS_HEvar_Mapman[!GS_HEvar_Mapman[,1]=="",]
GS_HEvar_Mapman<-GS_HEvar_Mapman[!duplicated(GS_HEvar_Mapman$Gene.x),]

write.table(GS_HEvar_Mapman,"GS_HEvar_Mapman_wVarLD.txt",sep="\t",row.names=F,col.names=F,quote=F)

KK_GS<-KK_GS[!duplicated(KK_GS$Ahal_Gene),]
KK_GS_Mapman<-KK_GS[,c(1,21)]
KK_GS_Mapman<-KK_GS_Mapman[!KK_GS_Mapman[,1]=="",]
KK_GS_Mapman<-KK_GS_Mapman[!duplicated(KK_GS_Mapman$Gene),]
write.table(KK_GS_Mapman,"KK_GS_Mapman_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

LB_GS<-LB_GS[!duplicated(LB_GS$Ahal_Gene),]
LB_GS_Mapman<-LB_GS[,c(1,21)]
LB_GS_Mapman<-LB_GS_Mapman[!LB_GS_Mapman[,1]=="",]
LB_GS_Mapman<-LB_GS_Mapman[!duplicated(LB_GS_Mapman$Gene),]
write.table(LB_GS_Mapman,"LB_GS_Mapman_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

NP_GS<-NP_GS[!duplicated(NP_GS$Ahal_Gene),]
NP_GS_Mapman<-NP_GS[,c(1,21)]
NP_GS_Mapman<-NP_GS_Mapman[!NP_GS_Mapman[,1]=="",]
NP_GS_Mapman<-NP_GS_Mapman[!duplicated(NP_GS_Mapman$Gene),]
write.table(NP_GS_Mapman,"NP_GS_Mapman_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

KK_GS_HEvar<-KK_GS_HEvar[!duplicated(KK_GS_HEvar$Ahal_Gene),]
KK_GS_HEvar_Mapman<-KK_GS_HEvar[,c(1,9)]
KK_GS_HEvar_Mapman<-KK_GS_HEvar_Mapman[!KK_GS_HEvar_Mapman[,1]=="",]
KK_GS_HEvar_Mapman<-KK_GS_HEvar_Mapman[!duplicated(KK_GS_HEvar_Mapman$Gene),]
write.table(KK_GS_HEvar_Mapman,"KK_GS_HEvar_Mapman_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

LB_GS_HEvar<-LB_GS_HEvar[!duplicated(LB_GS_HEvar$Ahal_Gene),]
LB_GS_HEvar_Mapman<-LB_GS_HEvar[,c(1,9)]
LB_GS_HEvar_Mapman<-LB_GS_HEvar_Mapman[!LB_GS_HEvar_Mapman[,1]=="",]
LB_GS_HEvar_Mapman<-LB_GS_HEvar_Mapman[!duplicated(LB_GS_HEvar_Mapman$Gene),]
write.table(LB_GS_HEvar_Mapman,"LB_GS_HEvar_Mapman_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

NP_GS_HEvar<-NP_GS_HEvar[!duplicated(NP_GS_HEvar$Ahal_Gene),]
NP_GS_HEvar_Mapman<-NP_GS_HEvar[,c(1,9)]
NP_GS_HEvar_Mapman<-NP_GS_HEvar_Mapman[!NP_GS_HEvar_Mapman[,1]=="",]
NP_GS_HEvar_Mapman<-NP_GS_HEvar_Mapman[!duplicated(NP_GS_HEvar_Mapman$Gene),]
write.table(NP_GS_HEvar_Mapman,"NP_GS_HEvar_Mapman_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)


Thal_ortholist=read.delim("F://Cu_project//Genome_scans//Genom_scans_Cu//Krom-Kosi//OG_Athaliana.table", header=TRUE)
Thal_ortholist$Ath<-trimws(as.character(Thal_ortholist$Ath))
Thal_ortholist$Ath<-substr(as.character(Thal_ortholist$Ath),1,9)

#25447
All<-data.frame(Thal_ortholist$Ath,1:length(Thal_ortholist$Ath))
All<-All[!duplicated(All[,1]),]
#25447
write.table(All,"All_GS_wVarLD_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)


require(xlsx)
All_GS_n<-nrow(read.table("All_GS_wVarLD_dedup.txt",sep="\t",header=F))
All_GS_Mapman<-read.xlsx2("MapMan_All_wVarLD_dedup.xls",1,header=T)
All_GS<-All_GS_Mapman[All_GS_Mapman$bin %in% 1:38,]
All_GS<-droplevels(All_GS)
All_GS<-All_GS[order(All_GS$bin),]

KK_GS_n<-nrow(read.table("KK_GS_Mapman_wVarLD_dedup.txt",sep="\t",header=F))
KK_GS_Mapman<-read.xlsx2("MapMan_KK_GS_wVarLD_dedup.xls",1,header=T)
KK_GS<-KK_GS_Mapman[KK_GS_Mapman$bin %in% 1:38,]
KK_GS<-droplevels(KK_GS)
KK_GS<-KK_GS[order(KK_GS$bin),]

x=KK_GS[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KK_GS_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KK_GS_n))/(dat[2,1]/All_GS_n))
}

KK_GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KK_GS_fisher$q_value<-qvalue(KK_GS_fisher$p.value,lambda=0)$qvalues
KK_GS_fisher_sig<-KK_GS_fisher[KK_GS_fisher$q_value<=0.05,]
write.table(KK_GS_fisher_sig, "KK_GS_fisher_wVarLD_dedup.txt",sep="\t",row.names=F,quote=F)

LB_GS_n<-nrow(read.table("LB_GS_Mapman_wVarLD_dedup.txt",sep="\t",header=F))
LB_GS_Mapman<-read.xlsx2("MapMan_LB_GS_wVarLD_dedup.xls",1,header=T)
LB_GS<-LB_GS_Mapman[LB_GS_Mapman$bin %in% 1:38,]
LB_GS<-droplevels(LB_GS)
LB_GS<-LB_GS[order(LB_GS$bin),]

x=LB_GS[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- LB_GS_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(LB_GS_n))/(dat[2,1]/All_GS_n))
}

LB_GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
LB_GS_fisher$q_value<-qvalue(LB_GS_fisher$p.value)$qvalues
LB_GS_fisher_sig<-LB_GS_fisher[LB_GS_fisher$q_value<=0.05,]
write.table(LB_GS_fisher_sig, "LB_GS_fisher_wVarLD_dedup.txt",sep="\t",row.names=F,quote=F)

NP_GS_n<-nrow(read.table("NP_GS_Mapman_wVarLD_dedup.txt",sep="\t",header=F))
NP_GS_Mapman<-read.xlsx2("MapMan_NP_GS_wVarLD_dedup.xls",1,header=T)
NP_GS<-NP_GS_Mapman[NP_GS_Mapman$bin %in% 1:38,]
NP_GS<-droplevels(NP_GS)
NP_GS<-NP_GS[order(NP_GS$bin),]

x=NP_GS[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- NP_GS_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(NP_GS_n))/(dat[2,1]/All_GS_n))
}

NP_GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
NP_GS_fisher$q_value<-qvalue(NP_GS_fisher$p.value)$qvalues
NP_GS_fisher_sig<-NP_GS_fisher[NP_GS_fisher$q_value<=0.05,]
write.table(NP_GS_fisher_sig, "NP_GS_fisher_wVarLD_dedup.txt",sep="\t",row.names=F,quote=F)

###############################################################################################
#Genome scans and high effect variants#
###############################################################################################

KK_GS_HEvar_n<-nrow(read.table("KK_GS_HEvar_Mapman_wVarLD.txt",sep="\t",header=F))
KK_GS_HEvar_Mapman<-read.xlsx2("MapMan_KK_GS_HEvar_wVarLD.xls",1,header=T)
KK_GS_HEvar<-KK_GS_HEvar_Mapman[KK_GS_HEvar_Mapman$bin %in% 1:38,]
KK_GS_HEvar<-droplevels(KK_GS_HEvar)
KK_GS_HEvar<-KK_GS_HEvar[order(KK_GS_HEvar$bin),]

x=KK_GS_HEvar[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KK_GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KK_GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

KK_GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KK_GS_HEvar_fisher$q_value<-qvalue(KK_GS_HEvar_fisher$p.value)$qvalues
KK_GS_HEvar_fisher_sig<-KK_GS_HEvar_fisher[KK_GS_HEvar_fisher$q_value<=0.05,]
write.table(KK_GS_HEvar_fisher_sig, "KK_GS_HEvar_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)

LB_GS_HEvar_n<-nrow(read.table("LB_GS_HEvar_Mapman_wVarLD.txt",sep="\t",header=F))
LB_GS_HEvar_Mapman<-read.xlsx2("MapMan_LB_GS_HEvar_wVarLD.xls",1,header=T)
LB_GS_HEvar<-LB_GS_HEvar_Mapman[LB_GS_HEvar_Mapman$bin %in% 1:38,]
LB_GS_HEvar<-droplevels(LB_GS_HEvar)
LB_GS_HEvar<-LB_GS_HEvar[order(LB_GS_HEvar$bin),]

x=LB_GS_HEvar[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- LB_GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(LB_GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

LB_GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
LB_GS_HEvar_fisher$q_value<-qvalue(LB_GS_HEvar_fisher$p.value)$qvalues
LB_GS_HEvar_fisher_sig<-LB_GS_HEvar_fisher[LB_GS_HEvar_fisher$q_value<=0.05,]
write.table(LB_GS_HEvar_fisher_sig, "LB_GS_HEvar_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)

NP_GS_HEvar_n<-nrow(read.table("NP_GS_HEvar_Mapman_wVarLD.txt",sep="\t",header=F))
NP_GS_HEvar_Mapman<-read.xlsx2("MapMan_NP_GS_HEvar_wVarLD.xls",1,header=T)
NP_GS_HEvar<-NP_GS_HEvar_Mapman[NP_GS_HEvar_Mapman$bin %in% 1:38,]
NP_GS_HEvar<-droplevels(NP_GS_HEvar)
NP_GS_HEvar<-NP_GS_HEvar[order(NP_GS_HEvar$bin),]

x=NP_GS_HEvar[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- NP_GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(NP_GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

NP_GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
NP_GS_HEvar_fisher$q_value<-qvalue(NP_GS_HEvar_fisher$p.value)$qvalues
NP_GS_HEvar_fisher_sig<-NP_GS_HEvar_fisher[NP_GS_HEvar_fisher$q_value<=0.05,]
write.table(NP_GS_HEvar_fisher_sig, "NP_GS_HEvar_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)





#test all subcategories
KK_GS_Mapman<-KK_GS_Mapman[order(KK_GS_Mapman$bin),]
LB_GS_Mapman<-LB_GS_Mapman[order(LB_GS_Mapman$bin),]
NP_GS_Mapman<-NP_GS_Mapman[order(NP_GS_Mapman$bin),]
All_GS_Mapman<-All_GS_Mapman[order(All_GS_Mapman$bin),]

x=KK_GS_Mapman[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS_Mapman
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KK_GS_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KK_GS_n))/(dat[2,1]/All_GS_n))
}

KK_GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KK_GS_fisher$q_value<-qvalue(KK_GS_fisher$p.value)$qvalues
KK_GS_fisher_sig<-KK_GS_fisher[KK_GS_fisher$q_value<=0.05,]
#write.table(KK_GS_fisher_sig, "KK_GS_fisher_sub_wVarLD.txt",sep="\t",row.names=F,quote=F)

x=LB_GS_Mapman[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS_Mapman
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- LB_GS_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(LB_GS_n))/(dat[2,1]/All_GS_n))
}

LB_GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
LB_GS_fisher$q_value<-qvalue(LB_GS_fisher$p.value)$qvalues
LB_GS_fisher_sig<-LB_GS_fisher[LB_GS_fisher$q_value<=0.05,]
#write.table(LB_GS_fisher_sig, "LB_GS_fisher_sub_wVarLD.txt",sep="\t",row.names=F,quote=F)

x=NP_GS_Mapman[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS_Mapman
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- NP_GS_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(NP_GS_n))/(dat[2,1]/All_GS_n))
}

NP_GS_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
NP_GS_fisher$q_value<-qvalue(NP_GS_fisher$p.value)$qvalues
NP_GS_fisher_sig<-NP_GS_fisher[NP_GS_fisher$q_value<=0.05,]
write.table(NP_GS_fisher_sig, "NP_GS_fisher_sub_wVarLD.txt",sep="\t",row.names=F,quote=F)

x=KK_GS_HEvar_Mapman[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS_Mapman
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KK_GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KK_GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

KK_GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KK_GS_HEvar_fisher$q_value<-qvalue(KK_GS_HEvar_fisher$p.value)$qvalues
KK_GS_HEvar_fisher_sig<-KK_GS_HEvar_fisher[KK_GS_HEvar_fisher$q_value<=0.05,]
write.table(KK_GS_HEvar_fisher_sig, "KK_GS_HEvar_fisher_sub_wVarLD.txt",sep="\t",row.names=F,quote=F)

x=LB_GS_HEvar_Mapman[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS_Mapman
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- LB_GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(LB_GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

LB_GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
LB_GS_HEvar_fisher$q_value<-qvalue(LB_GS_HEvar_fisher$p.value)$qvalues
LB_GS_HEvar_fisher_sig<-LB_GS_HEvar_fisher[LB_GS_HEvar_fisher$q_value<=0.05,]
#write.table(LB_GS_HEvar_fisher_sig, "LB_GS_HEvar_fisher_sub_wVarLD.txt",sep="\t",row.names=F,quote=F)

x=NP_GS_HEvar_Mapman[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS_Mapman
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- NP_GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(NP_GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

NP_GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
NP_GS_HEvar_fisher$q_value<-qvalue(NP_GS_HEvar_fisher$p.value)$qvalues
NP_GS_HEvar_fisher_sig<-NP_GS_HEvar_fisher[NP_GS_HEvar_fisher$q_value<=0.05,]
#write.table(NP_GS_HEvar_fisher_sig, "NP_GS_HEvar_fisher_sub_wVarLD.txt",sep="\t",row.names=F,quote=F)

GS_only_n<-nrow(read.table("GS_only_Mapman_wVarLD.txt",sep="\t",header=F))
GS_only_Mapman<-read.xlsx2("MapMan_GS_only_wVarLD.xls",1,header=T)
GS_only<-GS_only_Mapman[GS_only_Mapman$bin %in% 1:38,]
GS_only<-droplevels(GS_only)

x=GS_only[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- GS_only_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(GS_only_n))/(dat[2,1]/All_GS_n))
}

GS_only_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
GS_only_fisher$q_value<-qvalue(GS_only_fisher$p.value,lambda=0)$qvalues
GS_only_fisher_sig<-GS_only_fisher[GS_only_fisher$q_value<=0.05,]
#write.table(GS_only_fisher_sig, "GS_only_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)

GS_HEvar_n<-nrow(read.table("GS_HEvar_Mapman_wVarLD.txt",sep="\t",header=F))
GS_HEvar_Mapman<-read.xlsx2("MapMan_GS_HEvar_wVarLD.xls",1,header=T)
GS_HEvar<-GS_HEvar_Mapman[GS_HEvar_Mapman$bin %in% 1:38,]
GS_HEvar<-droplevels(GS_HEvar)

x=GS_HEvar[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
GS_HEvar_fisher$q_value<-qvalue(GS_HEvar_fisher$p.value)$qvalues
GS_HEvar_fisher_sig<-GS_HEvar_fisher[GS_HEvar_fisher$q_value<=0.05,]
#write.table(GS_HEvar_fisher_sig, "GS_HEvar_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)


#all bins
GS_only_n<-nrow(read.table("GS_only_Mapman_wVarLD.txt",sep="\t",header=F))
GS_only_Mapman<-read.xlsx2("MapMan_GS_only_wVarLD.xls",1,header=T)

x=GS_only[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- GS_only_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(GS_only_n))/(dat[2,1]/All_GS_n))
}

GS_only_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
GS_only_fisher$q_value<-qvalue(GS_only_fisher$p.value,lambda=0)$qvalues
GS_only_fisher_sig<-GS_only_fisher[GS_only_fisher$q_value<=0.05,]
#write.table(GS_only_fisher_sig, "GS_only_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)

GS_HEvar_n<-nrow(read.table("GS_HEvar_Mapman_wVarLD.txt",sep="\t",header=F))
GS_HEvar_Mapman<-read.xlsx2("MapMan_GS_HEvar_wVarLD.xls",1,header=T)

x=GS_HEvar[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_GS
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- GS_HEvar_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_GS_n-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(GS_HEvar_n))/(dat[2,1]/All_GS_n))
}

GS_HEvar_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
GS_HEvar_fisher$q_value<-qvalue(GS_HEvar_fisher$p.value)$qvalues
GS_HEvar_fisher_sig<-GS_HEvar_fisher[GS_HEvar_fisher$q_value<=0.05,]
#write.table(GS_HEvar_fisher_sig, "GS_HEvar_fisher_wVarLD.txt",sep="\t",row.names=F,quote=F)

