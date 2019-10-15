require(xlsx)
#Mapman input files to get input gene dataset because of overlap between categories
All_n_S<-nrow(read.table("All_cpm_Mapman_stricter_S_dedup.txt",sep="\t",header=F))
All_n_R<-nrow(read.table("All_cpm_Mapman_stricter_R_dedup.txt",sep="\t",header=F))
All_Mapman_R<-read.xlsx2("MapMan_All_stricter_R.xls",1,header=T)
All_Mapman_S<-read.xlsx2("MapMan_All_stricter_S.xls",1,header=T)
All_R<-All_Mapman_R[All_Mapman_R$bin %in% 1:38,]
All_R<-droplevels(All_R)
All_S<-All_Mapman_S[All_Mapman_S$bin %in% 1:38,]
All_S<-droplevels(All_S)
All_R<-All_R[order(All_R$bin),]
All_S<-All_S[order(All_S$bin),]


KromKosi_S_normal_Cu_up_n<-nrow(read.table("Krom_Kosi_S_normal_Cu_Mapman_stricter_up_dedup.txt",sep="\t",header=F))
KromKosi_S_normal_Cu_down_n<-nrow(read.table("Krom_Kosi_S_normal_Cu_Mapman_stricter_down_dedup.txt",sep="\t",header=F))
KromKosi_R_normal_Cu_up_n<-nrow(read.table("Krom_Kosi_R_normal_Cu_Mapman_stricter_up_dedup.txt",sep="\t",header=F))
KromKosi_R_normal_Cu_down_n<-nrow(read.table("Krom_Kosi_R_normal_Cu_Mapman_stricter_down_dedup.txt",sep="\t",header=F))

KromKosi_S_normal_Cu_Mapman_up<-read.xlsx2("MapMan_Krom_Kosi_S_normal_Cu_stricter_up.xls",1,header=T)
KromKosi_S_normal_Cu_Mapman_down<-read.xlsx2("MapMan_Krom_Kosi_S_normal_Cu_stricter_down.xls",1,header=T)
KromKosi_R_normal_Cu_Mapman_up<-read.xlsx2("MapMan_Krom_Kosi_R_normal_Cu_stricter_up.xls",1,header=T)
KromKosi_R_normal_Cu_Mapman_down<-read.xlsx2("MapMan_Krom_Kosi_R_normal_Cu_stricter_down.xls",1,header=T)
KromKosi_S_normal_Cu_up<-KromKosi_S_normal_Cu_Mapman_up[KromKosi_S_normal_Cu_Mapman_up$bin %in% 1:38,]
KromKosi_S_normal_Cu_up<-droplevels(KromKosi_S_normal_Cu_up)
KromKosi_R_normal_Cu_up<-KromKosi_R_normal_Cu_Mapman_up[KromKosi_R_normal_Cu_Mapman_up$bin %in% 1:38,]
KromKosi_R_normal_Cu_up<-droplevels(KromKosi_R_normal_Cu_up)
KromKosi_S_normal_Cu_down<-KromKosi_S_normal_Cu_Mapman_down[KromKosi_S_normal_Cu_Mapman_down$bin %in% 1:38,]
KromKosi_S_normal_Cu_down<-droplevels(KromKosi_S_normal_Cu_down)
KromKosi_R_normal_Cu_down<-KromKosi_R_normal_Cu_Mapman_down[KromKosi_R_normal_Cu_Mapman_down$bin %in% 1:38,]
KromKosi_R_normal_Cu_down<-droplevels(KromKosi_R_normal_Cu_down)

KromKosi_S_high_Cu_up_n<-nrow(read.table("Krom_Kosi_S_high_Cu_Mapman_stricter_up_dedup.txt",sep="\t",header=F))
KromKosi_S_high_Cu_down_n<-nrow(read.table("Krom_Kosi_S_high_Cu_Mapman_stricter_down_dedup.txt",sep="\t",header=F))
KromKosi_R_high_Cu_up_n<-nrow(read.table("Krom_Kosi_R_high_Cu_Mapman_stricter_up_dedup.txt",sep="\t",header=F))
KromKosi_R_high_Cu_down_n<-nrow(read.table("Krom_Kosi_R_high_Cu_Mapman_stricter_down_dedup.txt",sep="\t",header=F))
KromKosi_S_high_Cu_Mapman_up<-read.xlsx2("MapMan_Krom_Kosi_S_high_Cu_stricter_up.xls",1,header=T)
KromKosi_S_high_Cu_Mapman_down<-read.xlsx2("MapMan_Krom_Kosi_S_high_Cu_stricter_down.xls",1,header=T)
KromKosi_R_high_Cu_Mapman_up<-read.xlsx2("MapMan_Krom_Kosi_R_high_Cu_stricter_up.xls",1,header=T)
KromKosi_R_high_Cu_Mapman_down<-read.xlsx2("MapMan_Krom_Kosi_R_high_Cu_stricter_down.xls",1,header=T)
KromKosi_S_high_Cu_up<-KromKosi_S_high_Cu_Mapman_up[KromKosi_S_high_Cu_Mapman_up$bin %in% 1:38,]
KromKosi_S_high_Cu_up<-droplevels(KromKosi_S_high_Cu_up)
KromKosi_R_high_Cu_up<-KromKosi_R_high_Cu_Mapman_up[KromKosi_R_high_Cu_Mapman_up$bin %in% 1:38,]
KromKosi_R_high_Cu_up<-droplevels(KromKosi_R_high_Cu_up)
KromKosi_S_high_Cu_down<-KromKosi_S_high_Cu_Mapman_down[KromKosi_S_high_Cu_Mapman_down$bin %in% 1:38,]
KromKosi_S_high_Cu_down<-droplevels(KromKosi_S_high_Cu_down)
KromKosi_R_high_Cu_down<-KromKosi_R_high_Cu_Mapman_down[KromKosi_R_high_Cu_Mapman_down$bin %in% 1:38,]
KromKosi_R_high_Cu_down<-droplevels(KromKosi_R_high_Cu_down)


KromKosi_S_high_Cu_up<-KromKosi_S_high_Cu_up[order(KromKosi_S_high_Cu_up$bin),]
KromKosi_R_high_Cu_up<-KromKosi_R_high_Cu_up[order(KromKosi_R_high_Cu_up$bin),]
KromKosi_S_high_Cu_down<-KromKosi_S_high_Cu_down[order(KromKosi_S_high_Cu_down$bin),]
KromKosi_R_high_Cu_down<-KromKosi_R_high_Cu_down[order(KromKosi_R_high_Cu_down$bin),]
KromKosi_S_normal_Cu_up<-KromKosi_S_normal_Cu_up[order(KromKosi_S_normal_Cu_up$bin),]
KromKosi_R_normal_Cu_up<-KromKosi_R_normal_Cu_up[order(KromKosi_R_normal_Cu_up$bin),]
KromKosi_S_normal_Cu_down<-KromKosi_S_normal_Cu_down[order(KromKosi_S_normal_Cu_down$bin),]
KromKosi_R_normal_Cu_down<-KromKosi_R_normal_Cu_down[order(KromKosi_R_normal_Cu_down$bin),]


#KromKosi_S_normal_Cu_Mapman_up[grep("^26",KromKosi_S_normal_Cu_Mapman_up$bin),]

x=KromKosi_S_normal_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_normal_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KromKosi_S_normal_Cu_up_n))/(dat[2,1]/All_n_S))
}

KromKosi_S_normal_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KromKosi_S_normal_Cu_up_fisher$q_value<-qvalue(KromKosi_S_normal_Cu_up_fisher$p.value)$qvalues
KromKosi_S_normal_Cu_up_fisher_sig<-KromKosi_S_normal_Cu_up_fisher[KromKosi_S_normal_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_S_normal_Cu_up_fisher_sig, "KromKosi_S_normal_Cu_up_fisher_stricter_allS_dedup.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_normal_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_normal_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KromKosi_S_normal_Cu_down_n))/(dat[2,1]/All_n_S))
}

KromKosi_S_normal_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KromKosi_S_normal_Cu_down_fisher$q_value<-qvalue(KromKosi_S_normal_Cu_down_fisher$p.value)$qvalues
KromKosi_S_normal_Cu_down_fisher_sig<-KromKosi_S_normal_Cu_down_fisher[KromKosi_S_normal_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_S_normal_Cu_down_fisher_sig, "KromKosi_S_normal_Cu_down_fisher_stricter_allS_dedup.txt",sep="\t",row.names=F,quote=F)


x=KromKosi_R_normal_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
p.value <- c()
Fold_enrichment=c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_normal_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_normal_Cu_up_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_normal_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_normal_Cu_up_fisher$q_value<-qvalue(KromKosi_R_normal_Cu_up_fisher$p.value)$qvalues
KromKosi_R_normal_Cu_up_fisher_sig<-KromKosi_R_normal_Cu_up_fisher[KromKosi_R_normal_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_R_normal_Cu_up_fisher_sig, "KromKosi_R_normal_Cu_up_fisher_stricter_allR_dedup.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_R_normal_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
p.value <- c()
Fold_enrichment=c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_normal_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_normal_Cu_down_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_normal_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_normal_Cu_down_fisher$q_value<-qvalue(KromKosi_R_normal_Cu_down_fisher$p.value)$qvalues
KromKosi_R_normal_Cu_down_fisher_sig<-KromKosi_R_normal_Cu_down_fisher[KromKosi_R_normal_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_R_normal_Cu_down_fisher_sig, "KromKosi_R_normal_Cu_down_fisher_stricter_allR_dedup.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_high_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_high_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_S_high_Cu_up_n)/(dat[2,1]/All_n_S))
}

KromKosi_S_high_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_S_high_Cu_up_fisher$q_value<-qvalue(KromKosi_S_high_Cu_up_fisher$p.value)$qvalues
KromKosi_S_high_Cu_up_fisher_sig<-KromKosi_S_high_Cu_up_fisher[KromKosi_S_high_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_S_high_Cu_up_fisher_sig, "KromKosi_S_high_Cu_up_fisher_stricter_allS_dedup.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_high_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_high_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_S_high_Cu_down_n)/(dat[2,1]/All_n_S))
}

KromKosi_S_high_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_S_high_Cu_down_fisher$q_value<-qvalue(KromKosi_S_high_Cu_down_fisher$p.value)$qvalues
KromKosi_S_high_Cu_down_fisher_sig<-KromKosi_S_high_Cu_down_fisher[KromKosi_S_high_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_S_high_Cu_down_fisher_sig, "KromKosi_S_high_Cu_down_fisher_stricter_allS_dedup.txt",sep="\t",row.names=F,quote=F)


x=KromKosi_R_high_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_high_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_high_Cu_up_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_high_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_high_Cu_up_fisher$q_value<-qvalue(KromKosi_R_high_Cu_up_fisher$p.value)$qvalues
KromKosi_R_high_Cu_up_fisher_sig<-KromKosi_R_high_Cu_up_fisher[KromKosi_R_high_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_R_high_Cu_up_fisher_sig, "KromKosi_R_high_Cu_up_fisher_stricter_allR_dedup.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_R_high_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_high_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_high_Cu_down_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_high_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_high_Cu_down_fisher$q_value<-qvalue(KromKosi_R_high_Cu_down_fisher$p.value)$qvalues
KromKosi_R_high_Cu_down_fisher_sig<-KromKosi_R_high_Cu_down_fisher[KromKosi_R_high_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_R_high_Cu_down_fisher_sig, "KromKosi_R_high_Cu_down_fisher_stricter_allR_dedup.txt",sep="\t",row.names=F,quote=F)

bpKromKosi_R_normal_Cu_down_fisher_sig<-barplot(sort(-log10(KromKosi_R_normal_Cu_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_R_normal_Cu_down_fisher_sig$q_value))),width=1.2)
bpKromKosi_R_normal_Cu_up_fisher_sig<-barplot(sort(-log10(KromKosi_R_normal_Cu_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_R_normal_Cu_up_fisher_sig$q_value))),width=1.2)
bpKromKosi_S_normal_Cu_down_fisher_sig<-barplot(sort(-log10(KromKosi_S_normal_Cu_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_S_normal_Cu_down_fisher_sig$q_value))),width=1.2)
bpKromKosi_S_normal_Cu_up_fisher_sig<-barplot(sort(-log10(KromKosi_S_normal_Cu_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_S_normal_Cu_up_fisher_sig$q_value))),width=1.2)


pdf('GOenrichment_Mapman_candidates_const_stricter_upanddownsep_qvalue.pdf', width=8, height=6,paper="special",pointsize=15)
par(oma=c(6,3,2,1))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(1,0.75),widths=c(1,0.84))
par(mar=c(2,10,1,0),xpd=T)

barplot(sort(-log10(KromKosi_R_normal_Cu_down_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,xlim=c(0,35),width=1.2,axes=F,col="blue",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Gluconeogenesis /\nGlyoxylate cycle","Sec. metabol.","Amino acid metabol.","Hormone metabol.","Stress","Misc."),at=bpKromKosi_R_normal_Cu_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=1,xlim=c(0,15),width=1.2)
#mtext(text=expression("Transcript abundance lower in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(KromKosi_S_normal_Cu_down_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,xlim=c(0,25),width=1.2,col="blue",axes=F,cex.axis=1.3)
axis(side=2,tick=F,labels=c("Hormone metabol.","Secondary metabol.","N-metabol.","Misc."),at=bpKromKosi_S_normal_Cu_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(KromKosi_R_normal_Cu_up_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,35),col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Misc.","Cell","Stress","Secondary metabol."),at=bpKromKosi_R_normal_Cu_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
#mtext(text=expression("Transcript abundance higher in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)

barplot(sort(-log10(KromKosi_S_normal_Cu_up_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,25),col="red",cex.axis=1.3)
axis(side=2,tick=F,labels="Stress",at=bpKromKosi_S_normal_Cu_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(-log[10](adjusted~p-value)),side=1,cex=1.3,outer=T,line=1)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend=c("Transcript abundance lower in Krom than Kosi","Transcript abundance higher in Krom than Kosi"),fill=c("blue","red"),xpd=T,bty="n",inset=c(0,0),cex=1.3)

dev.off()







bpKromKosi_R_high_Cu_down_fisher_sig<-barplot(sort(-log10(KromKosi_R_high_Cu_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_R_high_Cu_down_fisher_sig$q_value))),width=1.2)
bpKromKosi_R_high_Cu_up_fisher_sig<-barplot(sort(-log10(KromKosi_R_high_Cu_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_R_high_Cu_up_fisher_sig$q_value))),width=1.2)
bpKromKosi_S_high_Cu_down_fisher_sig<-barplot(sort(-log10(KromKosi_S_high_Cu_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_S_high_Cu_down_fisher_sig$q_value))),width=1.2)
bpKromKosi_S_high_Cu_up_fisher_sig<-barplot(sort(-log10(KromKosi_S_high_Cu_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_S_high_Cu_up_fisher_sig$q_value))),width=1.2)


pdf('GOenrichment_Mapman_candidates_const_stricter_upanddownsep_qvalue_high.pdf', width=8, height=6,paper="special",pointsize=15)
par(oma=c(6,3,2,1))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(1,0.75),widths=c(1,0.84))
par(mar=c(2,10,1,0),xpd=T)

barplot(sort(-log10(KromKosi_R_high_Cu_down_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,xlim=c(0,35),width=1.2,axes=F,col="blue",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Hormone metabol.","Stress","Misc."),at=bpKromKosi_R_high_Cu_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=1,xlim=c(0,15),width=1.2)
#mtext(text=expression("Transcript abundance lower in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(KromKosi_S_high_Cu_down_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,xlim=c(0,25),width=1.2,col="blue",axes=F,cex.axis=1.3)
axis(side=2,tick=F,labels=c("Stress","Misc."),at=bpKromKosi_S_high_Cu_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(KromKosi_R_high_Cu_up_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,35),col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Stress"),at=bpKromKosi_R_high_Cu_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
#mtext(text=expression("Transcript abundance higher in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)

barplot(sort(-log10(KromKosi_S_high_Cu_up_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,25),col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=c("Misc.","Stress"),at=bpKromKosi_S_high_Cu_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(-log[10](adjusted~p-value)),side=1,cex=1.3,outer=T,line=1)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend=c("Transcript abundance lower in Krom than Kosi","Transcript abundance higher in Krom than Kosi"),fill=c("blue","red"),xpd=T,bty="n",inset=c(0,0),cex=1.3)

dev.off()



#Misc
All_Mapman_R<-read.xlsx2("MapMan_All_stricter_R.xls",1,header=T)
All_Mapman_S<-read.xlsx2("MapMan_All_stricter_S.xls",1,header=T)
All_R<-All_Mapman_R[grep("^26",All_Mapman_R$bin),]
All_R<-droplevels(All_R)
All_S<-All_Mapman_S[grep("^26",All_Mapman_S$bin),]
All_S<-droplevels(All_S)

KromKosi_S_normal_Cu_Mapman_up<-read.xlsx2("MapMan_KromKosi_S_normal_Cu_up_stricter.xls",1,header=T)
KromKosi_S_normal_Cu_Mapman_down<-read.xlsx2("MapMan_KromKosi_S_normal_Cu_down_stricter.xls",1,header=T)
KromKosi_R_normal_Cu_Mapman_up<-read.xlsx2("MapMan_KromKosi_R_normal_Cu_up_stricter.xls",1,header=T)
KromKosi_R_normal_Cu_Mapman_down<-read.xlsx2("MapMan_KromKosi_R_normal_Cu_down_stricter.xls",1,header=T)
KromKosi_S_normal_Cu_up<-KromKosi_S_normal_Cu_Mapman_up[grep("^26",KromKosi_S_normal_Cu_Mapman_up$bin),]
KromKosi_S_normal_Cu_up<-droplevels(KromKosi_S_normal_Cu_up)
KromKosi_R_normal_Cu_up<-KromKosi_R_normal_Cu_Mapman_up[grep("^26",KromKosi_R_normal_Cu_Mapman_up$bin),]
KromKosi_R_normal_Cu_up<-droplevels(KromKosi_R_normal_Cu_up)
KromKosi_S_normal_Cu_down<-KromKosi_S_normal_Cu_Mapman_down[grep("^26",KromKosi_S_normal_Cu_Mapman_down$bin),]
KromKosi_S_normal_Cu_down<-droplevels(KromKosi_S_normal_Cu_down)
KromKosi_R_normal_Cu_down<-KromKosi_R_normal_Cu_Mapman_down[grep("^26",KromKosi_R_normal_Cu_Mapman_down$bin),]
KromKosi_R_normal_Cu_down<-droplevels(KromKosi_R_normal_Cu_down)

KromKosi_S_high_Cu_up_n<-nrow(read.table("Krom_Kosi_S_high_Cu_Mapman_stricter_up.txt",sep="\t",header=F))
KromKosi_S_high_Cu_down_n<-nrow(read.table("Krom_Kosi_S_high_Cu_Mapman_stricter_down.txt",sep="\t",header=F))
KromKosi_R_high_Cu_up_n<-nrow(read.table("Krom_Kosi_R_high_Cu_Mapman_stricter_up.txt",sep="\t",header=F))
KromKosi_R_high_Cu_down_n<-nrow(read.table("Krom_Kosi_R_high_Cu_Mapman_stricter_down.txt",sep="\t",header=F))
KromKosi_S_high_Cu_Mapman_up<-read.xlsx2("MapMan_KromKosi_S_high_Cu_up_stricter.xls",1,header=T)
KromKosi_S_high_Cu_Mapman_down<-read.xlsx2("MapMan_KromKosi_S_high_Cu_down_stricter.xls",1,header=T)
KromKosi_R_high_Cu_Mapman_up<-read.xlsx2("MapMan_KromKosi_R_high_Cu_up_stricter.xls",1,header=T)
KromKosi_R_high_Cu_Mapman_down<-read.xlsx2("MapMan_KromKosi_R_high_Cu_down_stricter.xls",1,header=T)
KromKosi_S_high_Cu_up<-KromKosi_S_high_Cu_Mapman_up[grep("^26",KromKosi_S_high_Cu_Mapman_up$bin),]
KromKosi_S_high_Cu_up<-droplevels(KromKosi_S_high_Cu_up)
KromKosi_R_high_Cu_up<-KromKosi_R_high_Cu_Mapman_up[grep("^26",KromKosi_R_high_Cu_Mapman_up$bin),]
KromKosi_R_high_Cu_up<-droplevels(KromKosi_R_high_Cu_up)
KromKosi_S_high_Cu_down<-KromKosi_S_high_Cu_Mapman_down[grep("^26",KromKosi_S_high_Cu_Mapman_down$bin),]
KromKosi_S_high_Cu_down<-droplevels(KromKosi_S_high_Cu_down)
KromKosi_R_high_Cu_down<-KromKosi_R_high_Cu_Mapman_down[grep("^26",KromKosi_R_high_Cu_Mapman_down$bin),]
KromKosi_R_high_Cu_down<-droplevels(KromKosi_R_high_Cu_down)


KromKosi_S_high_Cu_up<-KromKosi_S_high_Cu_up[order(KromKosi_S_high_Cu_up$bin),]
KromKosi_R_high_Cu_up<-KromKosi_R_high_Cu_up[order(KromKosi_R_high_Cu_up$bin),]
KromKosi_S_high_Cu_down<-KromKosi_S_high_Cu_down[order(KromKosi_S_high_Cu_down$bin),]
KromKosi_R_high_Cu_down<-KromKosi_R_high_Cu_down[order(KromKosi_R_high_Cu_down$bin),]
KromKosi_S_normal_Cu_up<-KromKosi_S_normal_Cu_up[order(KromKosi_S_normal_Cu_up$bin),]
KromKosi_R_normal_Cu_up<-KromKosi_R_normal_Cu_up[order(KromKosi_R_normal_Cu_up$bin),]
KromKosi_S_normal_Cu_down<-KromKosi_S_normal_Cu_down[order(KromKosi_S_normal_Cu_down$bin),]
KromKosi_R_normal_Cu_down<-KromKosi_R_normal_Cu_down[order(KromKosi_R_normal_Cu_down$bin),]


#KromKosi_S_normal_Cu_Mapman_up[grep("^26",KromKosi_S_normal_Cu_Mapman_up$bin),]

x=KromKosi_S_normal_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_normal_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KromKosi_S_normal_Cu_up_n))/(dat[2,1]/All_n_S))
}

KromKosi_S_normal_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KromKosi_S_normal_Cu_up_fisher$q_value<-qvalue(KromKosi_S_normal_Cu_up_fisher$p.value,lambda=0)$qvalues
KromKosi_S_normal_Cu_up_fisher_sig<-KromKosi_S_normal_Cu_up_fisher[KromKosi_S_normal_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_S_normal_Cu_up_fisher_sig, "KromKosi_S_normal_Cu_up_fisher_stricter_allS_Misc.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_normal_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_normal_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KromKosi_S_normal_Cu_down_n))/(dat[2,1]/All_n_S))
}

KromKosi_S_normal_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KromKosi_S_normal_Cu_down_fisher$q_value<-qvalue(KromKosi_S_normal_Cu_down_fisher$p.value)$qvalues
KromKosi_S_normal_Cu_down_fisher_sig<-KromKosi_S_normal_Cu_down_fisher[KromKosi_S_normal_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_S_normal_Cu_down_fisher_sig, "KromKosi_S_normal_Cu_down_fisher_stricter_allS.txt",sep="\t",row.names=F,quote=F)
KromKosi_S_normal_Cu_down_fisher_sig<-read.table("KromKosi_S_normal_Cu_down_fisher_stricter_allS_Misc.txt",sep="\t",header=T)


x=KromKosi_R_normal_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
p.value <- c()
Fold_enrichment=c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_normal_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_normal_Cu_up_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_normal_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_normal_Cu_up_fisher$q_value<-qvalue(KromKosi_R_normal_Cu_up_fisher$p.value)$qvalues
KromKosi_R_normal_Cu_up_fisher_sig<-KromKosi_R_normal_Cu_up_fisher[KromKosi_R_normal_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_R_normal_Cu_up_fisher_sig, "KromKosi_R_normal_Cu_up_fisher_stricter_allR.txt",sep="\t",row.names=F,quote=F)
KromKosi_R_normal_Cu_up_fisher_sig<-read.table("KromKosi_R_normal_Cu_up_fisher_stricter_allR_Misc.txt",sep="\t",header=T)

x=KromKosi_R_normal_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
p.value <- c()
Fold_enrichment=c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_normal_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_normal_Cu_down_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_normal_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_normal_Cu_down_fisher$q_value<-qvalue(KromKosi_R_normal_Cu_down_fisher$p.value)$qvalues
KromKosi_R_normal_Cu_down_fisher_sig<-KromKosi_R_normal_Cu_down_fisher[KromKosi_R_normal_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_R_normal_Cu_down_fisher_sig, "KromKosi_R_normal_Cu_down_fisher_stricter_allR.txt",sep="\t",row.names=F,quote=F)
KromKosi_R_normal_Cu_down_fisher_sig<-read.table("KromKosi_R_normal_Cu_down_fisher_stricter_allR_Misc.txt",sep="\t",header=T)


x=KromKosi_S_high_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_high_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_S_high_Cu_up_n)/(dat[2,1]/All_n_S))
}

KromKosi_S_high_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_S_high_Cu_up_fisher$q_value<-qvalue(KromKosi_S_high_Cu_up_fisher$p.value)$qvalues
KromKosi_S_high_Cu_up_fisher_sig<-KromKosi_S_high_Cu_up_fisher[KromKosi_S_high_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_S_high_Cu_up_fisher_sig, "KromKosi_S_high_Cu_up_fisher_stricter_allS_Misc.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_high_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_high_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_S_high_Cu_down_n)/(dat[2,1]/All_n_S))
}

KromKosi_S_high_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_S_high_Cu_down_fisher$q_value<-qvalue(KromKosi_S_high_Cu_down_fisher$p.value)$qvalues
KromKosi_S_high_Cu_down_fisher_sig<-KromKosi_S_high_Cu_down_fisher[KromKosi_S_high_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_S_high_Cu_down_fisher_sig, "KromKosi_S_high_Cu_down_fisher_stricter_allS_Misc.txt",sep="\t",row.names=F,quote=F)


x=KromKosi_R_high_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_high_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_high_Cu_up_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_high_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_high_Cu_up_fisher$q_value<-qvalue(KromKosi_R_high_Cu_up_fisher$p.value)$qvalues
KromKosi_R_high_Cu_up_fisher_sig<-KromKosi_R_high_Cu_up_fisher[KromKosi_R_high_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_R_high_Cu_up_fisher_sig, "KromKosi_R_high_Cu_up_fisher_stricter_allR_Misc.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_R_high_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_high_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_high_Cu_down_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_high_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_high_Cu_down_fisher$q_value<-qvalue(KromKosi_R_high_Cu_down_fisher$p.value)$qvalues
KromKosi_R_high_Cu_down_fisher_sig<-KromKosi_R_high_Cu_down_fisher[KromKosi_R_high_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_R_high_Cu_down_fisher_sig, "KromKosi_R_high_Cu_down_fisher_stricter_allR_Misc.txt",sep="\t",row.names=F,quote=F)

######
MapMan_Krom_Kosi_S_both_stricter_up_n<-nrow(read.table("Krom_Kosi_S_both_MapMan_stricter_up_dedup.txt",sep="\t",header=F))
MapMan_Krom_Kosi_S_both_stricter_down_n<-nrow(read.table("Krom_Kosi_S_both_MapMan_stricter_down_dedup.txt",sep="\t",header=F))
MapMan_Krom_Kosi_R_both_stricter_up_n<-nrow(read.table("Krom_Kosi_R_both_MapMan_stricter_up_dedup.txt",sep="\t",header=F))
MapMan_Krom_Kosi_R_both_stricter_down_n<-nrow(read.table("Krom_Kosi_R_both_MapMan_stricter_down_dedup.txt",sep="\t",header=F))
MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter_n<-nrow(read.table("Krom_Kosi_R_Diff_normal_Kromlower_MapMan_stricter_dedup.txt",sep="\t",header=F))
MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter_n<-nrow(read.table("Krom_Kosi_R_Diff_normal_Kromhigher_MapMan_stricter_dedup.txt",sep="\t",header=F))

MapMan_Krom_Kosi_S_both_stricter_up<-read.xlsx2("MapMan_Krom_Kosi_S_both_stricter_up.xls",1,header=T)
MapMan_Krom_Kosi_S_both_stricter_down<-read.xlsx2("MapMan_Krom_Kosi_S_both_stricter_down.xls",1,header=T)
MapMan_Krom_Kosi_R_both_stricter_up<-read.xlsx2("MapMan_Krom_Kosi_R_both_stricter_up.xls",1,header=T)
MapMan_Krom_Kosi_R_both_stricter_down<-read.xlsx2("MapMan_Krom_Kosi_R_both_stricter_down.xls",1,header=T)
MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter<-read.xlsx2("MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter.xls",1,header=T)
MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter<-read.xlsx2("MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter.xls",1,header=T)

MapMan_Krom_Kosi_S_both_stricter_up<-MapMan_Krom_Kosi_S_both_stricter_up[order(MapMan_Krom_Kosi_S_both_stricter_up$bin),]
MapMan_Krom_Kosi_S_both_stricter_down<-MapMan_Krom_Kosi_S_both_stricter_down[order(MapMan_Krom_Kosi_S_both_stricter_down$bin),]
MapMan_Krom_Kosi_R_both_stricter_up<-MapMan_Krom_Kosi_R_both_stricter_up[order(MapMan_Krom_Kosi_R_both_stricter_up$bin),]
MapMan_Krom_Kosi_R_both_stricter_down<-MapMan_Krom_Kosi_R_both_stricter_down[order(MapMan_Krom_Kosi_R_both_stricter_down$bin),]
MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter<-MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter[order(MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter$bin),]
MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter<-MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter[order(MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter$bin),]

MapMan_Krom_Kosi_S_both_stricter_up<-MapMan_Krom_Kosi_S_both_stricter_up[MapMan_Krom_Kosi_S_both_stricter_up$bin %in% 1:38,]
MapMan_Krom_Kosi_S_both_stricter_down<-MapMan_Krom_Kosi_S_both_stricter_down[MapMan_Krom_Kosi_S_both_stricter_down$bin %in% 1:38,]
MapMan_Krom_Kosi_R_both_stricter_up<-MapMan_Krom_Kosi_R_both_stricter_up[MapMan_Krom_Kosi_R_both_stricter_up$bin %in% 1:38,]
MapMan_Krom_Kosi_R_both_stricter_down<-MapMan_Krom_Kosi_R_both_stricter_down[MapMan_Krom_Kosi_R_both_stricter_down$bin %in% 1:38,]
MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter<-MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter[MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter$bin %in% 1:38,]
MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter<-MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter[MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter$bin %in% 1:38,]

MapMan_Krom_Kosi_S_both_stricter_up<-droplevels(MapMan_Krom_Kosi_S_both_stricter_up)
MapMan_Krom_Kosi_S_both_stricter_down<-droplevels(MapMan_Krom_Kosi_S_both_stricter_down)
MapMan_Krom_Kosi_R_both_stricter_up<-droplevels(MapMan_Krom_Kosi_R_both_stricter_up)
MapMan_Krom_Kosi_R_both_stricter_down<-droplevels(MapMan_Krom_Kosi_R_both_stricter_down)
MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter<-droplevels(MapMan_Krom_Kosi_R_Diff_normal_Kromlower_stricter)
MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter<-droplevels(MapMan_Krom_Kosi_R_Diff_normal_Kromhigher_stricter)



x=MapMan_Krom_Kosi_S_both_stricter_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- MapMan_Krom_Kosi_S_both_stricter_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(MapMan_Krom_Kosi_S_both_stricter_up_n))/(dat[2,1]/All_n_S))
}

MapMan_Krom_Kosi_S_both_stricter_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
MapMan_Krom_Kosi_S_both_stricter_up_fisher$q_value<-qvalue(MapMan_Krom_Kosi_S_both_stricter_up_fisher$p.value)$qvalues
MapMan_Krom_Kosi_S_both_stricter_up_fisher_sig<-MapMan_Krom_Kosi_S_both_stricter_up_fisher[MapMan_Krom_Kosi_S_both_stricter_up_fisher$q_value<=0.05,]
write.table(MapMan_Krom_Kosi_S_both_stricter_up_fisher_sig, "MapMan_Krom_Kosi_S_both_stricter_up_fisher_stricter_allS_dedup.txt",sep="\t",row.names=F,quote=F)

x=MapMan_Krom_Kosi_S_both_stricter_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- MapMan_Krom_Kosi_S_both_stricter_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(MapMan_Krom_Kosi_S_both_stricter_down_n))/(dat[2,1]/All_n_S))
}

MapMan_Krom_Kosi_S_both_stricter_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
MapMan_Krom_Kosi_S_both_stricter_down_fisher$q_value<-qvalue(MapMan_Krom_Kosi_S_both_stricter_down_fisher$p.value)$qvalues
MapMan_Krom_Kosi_S_both_stricter_down_fisher_sig<-MapMan_Krom_Kosi_S_both_stricter_down_fisher[MapMan_Krom_Kosi_S_both_stricter_down_fisher$q_value<=0.05,]
write.table(MapMan_Krom_Kosi_S_both_stricter_down_fisher_sig, "MapMan_Krom_Kosi_S_both_stricter_down_fisher_stricter_allS_dedup.txt",sep="\t",row.names=F,quote=F)


x=MapMan_Krom_Kosi_R_both_stricter_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- MapMan_Krom_Kosi_R_both_stricter_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(MapMan_Krom_Kosi_R_both_stricter_up_n))/(dat[2,1]/All_n_R))
}

MapMan_Krom_Kosi_R_both_stricter_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
MapMan_Krom_Kosi_R_both_stricter_up_fisher$q_value<-qvalue(MapMan_Krom_Kosi_R_both_stricter_up_fisher$p.value)$qvalues
MapMan_Krom_Kosi_R_both_stricter_up_fisher_sig<-MapMan_Krom_Kosi_R_both_stricter_up_fisher[MapMan_Krom_Kosi_R_both_stricter_up_fisher$q_value<=0.05,]
write.table(MapMan_Krom_Kosi_R_both_stricter_up_fisher_sig, "MapMan_Krom_Kosi_R_both_stricter_up_fisher_stricter_allR_dedup.txt",sep="\t",row.names=F,quote=F)

x=MapMan_Krom_Kosi_R_both_stricter_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- MapMan_Krom_Kosi_R_both_stricter_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(MapMan_Krom_Kosi_R_both_stricter_down_n))/(dat[2,1]/All_n_R))
}

MapMan_Krom_Kosi_R_both_stricter_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
MapMan_Krom_Kosi_R_both_stricter_down_fisher$q_value<-qvalue(MapMan_Krom_Kosi_R_both_stricter_down_fisher$p.value)$qvalues
MapMan_Krom_Kosi_R_both_stricter_down_fisher_sig<-MapMan_Krom_Kosi_R_both_stricter_down_fisher[MapMan_Krom_Kosi_R_both_stricter_down_fisher$q_value<=0.05,]
write.table(MapMan_Krom_Kosi_R_both_stricter_down_fisher_sig, "MapMan_Krom_Kosi_R_both_stricter_down_fisher_stricter_allR_dedup.txt",sep="\t",row.names=F,quote=F)






#Stress
All_R<-All_Mapman_R[grep("^20",All_Mapman_R$bin),]
All_R<-droplevels(All_R)
All_S<-All_Mapman_S[grep("^20",All_Mapman_S$bin),]
All_S<-droplevels(All_S)
All_R<-All_R[order(All_R$bin),]
All_S<-All_S[order(All_S$bin),]

KromKosi_S_normal_Cu_up<-KromKosi_S_normal_Cu_Mapman_up[grep("^20",KromKosi_S_normal_Cu_Mapman_up$bin),]
KromKosi_S_normal_Cu_up<-droplevels(KromKosi_S_normal_Cu_up)
KromKosi_R_normal_Cu_up<-KromKosi_R_normal_Cu_Mapman_up[grep("^20",KromKosi_R_normal_Cu_Mapman_up$bin),]
KromKosi_R_normal_Cu_up<-droplevels(KromKosi_R_normal_Cu_up)
KromKosi_S_normal_Cu_down<-KromKosi_S_normal_Cu_Mapman_down[grep("^20",KromKosi_S_normal_Cu_Mapman_down$bin),]
KromKosi_S_normal_Cu_down<-droplevels(KromKosi_S_normal_Cu_down)
KromKosi_R_normal_Cu_down<-KromKosi_R_normal_Cu_Mapman_down[grep("^20",KromKosi_R_normal_Cu_Mapman_down$bin),]
KromKosi_R_normal_Cu_down<-droplevels(KromKosi_R_normal_Cu_down)

KromKosi_S_high_Cu_up<-KromKosi_S_high_Cu_Mapman_up[grep("^20",KromKosi_S_high_Cu_Mapman_up$bin),]
KromKosi_S_high_Cu_up<-droplevels(KromKosi_S_high_Cu_up)
KromKosi_R_high_Cu_up<-KromKosi_R_high_Cu_Mapman_up[grep("^20",KromKosi_R_high_Cu_Mapman_up$bin),]
KromKosi_R_high_Cu_up<-droplevels(KromKosi_R_high_Cu_up)
KromKosi_S_high_Cu_down<-KromKosi_S_high_Cu_Mapman_down[grep("^20",KromKosi_S_high_Cu_Mapman_down$bin),]
KromKosi_S_high_Cu_down<-droplevels(KromKosi_S_high_Cu_down)
KromKosi_R_high_Cu_down<-KromKosi_R_high_Cu_Mapman_down[grep("^20",KromKosi_R_high_Cu_Mapman_down$bin),]
KromKosi_R_high_Cu_down<-droplevels(KromKosi_R_high_Cu_down)


KromKosi_S_high_Cu_up<-KromKosi_S_high_Cu_up[order(KromKosi_S_high_Cu_up$bin),]
KromKosi_R_high_Cu_up<-KromKosi_R_high_Cu_up[order(KromKosi_R_high_Cu_up$bin),]
KromKosi_S_high_Cu_down<-KromKosi_S_high_Cu_down[order(KromKosi_S_high_Cu_down$bin),]
KromKosi_R_high_Cu_down<-KromKosi_R_high_Cu_down[order(KromKosi_R_high_Cu_down$bin),]
KromKosi_S_normal_Cu_up<-KromKosi_S_normal_Cu_up[order(KromKosi_S_normal_Cu_up$bin),]
KromKosi_R_normal_Cu_up<-KromKosi_R_normal_Cu_up[order(KromKosi_R_normal_Cu_up$bin),]
KromKosi_S_normal_Cu_down<-KromKosi_S_normal_Cu_down[order(KromKosi_S_normal_Cu_down$bin),]
KromKosi_R_normal_Cu_down<-KromKosi_R_normal_Cu_down[order(KromKosi_R_normal_Cu_down$bin),]

KromKosi_S_high_Cu_up<-droplevels(KromKosi_S_high_Cu_up)
KromKosi_R_high_Cu_up<-droplevels(KromKosi_R_high_Cu_up)
KromKosi_S_high_Cu_down<-droplevels(KromKosi_S_high_Cu_down)
KromKosi_R_high_Cu_down<-droplevels(KromKosi_R_high_Cu_down)
KromKosi_S_normal_Cu_up<-droplevels(KromKosi_S_normal_Cu_up)
KromKosi_R_normal_Cu_up<-droplevels(KromKosi_R_normal_Cu_up)
KromKosi_S_normal_Cu_down<-droplevels(KromKosi_S_normal_Cu_down)
KromKosi_R_normal_Cu_down<-droplevels(KromKosi_R_normal_Cu_down)

#KromKosi_S_normal_Cu_Mapman_up[grep("^26",KromKosi_S_normal_Cu_Mapman_up$bin),]

x=KromKosi_S_normal_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_normal_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KromKosi_S_normal_Cu_up_n))/(dat[2,1]/All_n_S))
}

KromKosi_S_normal_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KromKosi_S_normal_Cu_up_fisher$q_value<-qvalue(KromKosi_S_normal_Cu_up_fisher$p.value)$qvalues
KromKosi_S_normal_Cu_up_fisher_sig<-KromKosi_S_normal_Cu_up_fisher[KromKosi_S_normal_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_S_normal_Cu_up_fisher_sig, "KromKosi_S_normal_Cu_up_fisher_stricter_allS_dedup_stress.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_normal_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_normal_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/(KromKosi_S_normal_Cu_down_n))/(dat[2,1]/All_n_S))
}

KromKosi_S_normal_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
require(qvalue)
KromKosi_S_normal_Cu_down_fisher$q_value<-qvalue(KromKosi_S_normal_Cu_down_fisher$p.value)$qvalues
KromKosi_S_normal_Cu_down_fisher_sig<-KromKosi_S_normal_Cu_down_fisher[KromKosi_S_normal_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_S_normal_Cu_down_fisher_sig, "KromKosi_S_normal_Cu_down_fisher_stricter_allS_dedup_stress.txt",sep="\t",row.names=F,quote=F)


x=KromKosi_R_normal_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
p.value <- c()
Fold_enrichment=c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_normal_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_normal_Cu_up_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_normal_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_normal_Cu_up_fisher$q_value<-qvalue(KromKosi_R_normal_Cu_up_fisher$p.value,lambda=0)$qvalues
KromKosi_R_normal_Cu_up_fisher_sig<-KromKosi_R_normal_Cu_up_fisher[KromKosi_R_normal_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_R_normal_Cu_up_fisher_sig, "KromKosi_R_normal_Cu_up_fisher_stricter_allR_dedup_stress.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_R_normal_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
p.value <- c()
Fold_enrichment=c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_normal_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_normal_Cu_down_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_normal_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_normal_Cu_down_fisher$q_value<-qvalue(KromKosi_R_normal_Cu_down_fisher$p.value,lambda=0)$qvalues
KromKosi_R_normal_Cu_down_fisher_sig<-KromKosi_R_normal_Cu_down_fisher[KromKosi_R_normal_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_R_normal_Cu_down_fisher_sig, "KromKosi_R_normal_Cu_down_fisher_stricter_allR_dedup_stress.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_high_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_high_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_S_high_Cu_up_n)/(dat[2,1]/All_n_S))
}

KromKosi_S_high_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_S_high_Cu_up_fisher$q_value<-qvalue(KromKosi_S_high_Cu_up_fisher$p.value,lambda=0)$qvalues
KromKosi_S_high_Cu_up_fisher_sig<-KromKosi_S_high_Cu_up_fisher[KromKosi_S_high_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_S_high_Cu_up_fisher_sig, "KromKosi_S_high_Cu_up_fisher_stricter_allS_dedup_stress.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_S_high_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_S
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_S_high_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_S-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_S_high_Cu_down_n)/(dat[2,1]/All_n_S))
}

KromKosi_S_high_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_S_high_Cu_down_fisher$q_value<-qvalue(KromKosi_S_high_Cu_down_fisher$p.value,lambda=0)$qvalues
KromKosi_S_high_Cu_down_fisher_sig<-KromKosi_S_high_Cu_down_fisher[KromKosi_S_high_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_S_high_Cu_down_fisher_sig, "KromKosi_S_high_Cu_down_fisher_stricter_allS_dedup_stress.txt",sep="\t",row.names=F,quote=F)


x=KromKosi_R_high_Cu_up[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_high_Cu_up_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_high_Cu_up_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_high_Cu_up_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_high_Cu_up_fisher$q_value<-qvalue(KromKosi_R_high_Cu_up_fisher$p.value,lambda=0)$qvalues
KromKosi_R_high_Cu_up_fisher_sig<-KromKosi_R_high_Cu_up_fisher[KromKosi_R_high_Cu_up_fisher$q_value<=0.05,]
write.table(KromKosi_R_high_Cu_up_fisher_sig, "KromKosi_R_high_Cu_up_fisher_stricter_allR_dedup_stress.txt",sep="\t",row.names=F,quote=F)

x=KromKosi_R_high_Cu_down[,1:3]
x[,3]<-as.numeric(as.character(x[,3]))
y=All_R
y[,3]<-as.numeric(as.character(y[,3]))

dat <- matrix(NA, 2,2)
Fold_enrichment=c()
p.value <- c()
for (i in levels(x$bin)) 
{  
  dat[1,1] <- x$elements[x$bin==i]
  dat[1,2] <- KromKosi_R_high_Cu_down_n-dat[1,1]
  dat[2,1] <- y$elements[y$bin==i]-dat[1,1]
  dat[2,2] <- All_n_R-y$elements[y$bin==i]-dat[1,2]

  fit <- fisher.test(dat,alternative="greater")#greater for overrepresentation
  p.value <- rbind(p.value,fit$p.value)
  Fold_enrichment=rbind(Fold_enrichment,(dat[1,1]/KromKosi_R_high_Cu_down_n)/(dat[2,1]/All_n_R))
}

KromKosi_R_high_Cu_down_fisher <- cbind(x, p.value,Fold_enrichment) 
KromKosi_R_high_Cu_down_fisher$q_value<-qvalue(KromKosi_R_high_Cu_down_fisher$p.value,lambda=0)$qvalues
KromKosi_R_high_Cu_down_fisher_sig<-KromKosi_R_high_Cu_down_fisher[KromKosi_R_high_Cu_down_fisher$q_value<=0.05,]
write.table(KromKosi_R_high_Cu_down_fisher_sig, "KromKosi_R_high_Cu_down_fisher_stricter_allR_dedup_stress.txt",sep="\t",row.names=F,quote=F)

bpKromKosi_R_high_Cu_down_fisher_sig<-barplot(sort(-log10(KromKosi_R_high_Cu_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_R_high_Cu_down_fisher_sig$q_value))),width=1.2)
bpKromKosi_R_high_Cu_up_fisher_sig<-barplot(sort(-log10(KromKosi_R_high_Cu_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_R_high_Cu_up_fisher_sig$q_value))),width=1.2)
bpKromKosi_S_high_Cu_down_fisher_sig<-barplot(sort(-log10(KromKosi_S_high_Cu_down_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_S_high_Cu_down_fisher_sig$q_value))),width=1.2)
bpKromKosi_S_high_Cu_up_fisher_sig<-barplot(sort(-log10(KromKosi_S_high_Cu_up_fisher_sig$q_value)),las=3,xlim=c(0,max(-log10(KromKosi_S_high_Cu_up_fisher_sig$q_value))),width=1.2)


pdf('GOenrichment_Mapman_candidates_const_stricter_upanddownsep_qvalue_high.pdf', width=8, height=6,paper="special",pointsize=15)
par(oma=c(6,3,2,1))
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),heights=c(1,0.75),widths=c(1,0.84))
par(mar=c(2,10,1,0),xpd=T)

barplot(sort(-log10(KromKosi_R_high_Cu_down_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,xlim=c(0,35),width=1.2,axes=F,col="blue",cex.axis=1.3)
axis(side=2,tick=F,labels=KromKosi_R_high_Cu_down_fisher_sig$name[order(-KromKosi_R_high_Cu_down_fisher_sig$q_value)],at=bpKromKosi_R_high_Cu_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Enriched~Mapman~categories),side=2,cex=1.3,outer=T,line=1,xlim=c(0,15),width=1.2)
#mtext(text=expression("Transcript abundance lower in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(KromKosi_S_high_Cu_down_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,xlim=c(0,25),width=1.2,col="blue",axes=F,cex.axis=1.3)
axis(side=2,tick=F,labels=KromKosi_S_high_Cu_down_fisher_sig$name[order(-KromKosi_S_high_Cu_down_fisher_sig$q_value)],at=bpKromKosi_S_high_Cu_down_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)

barplot(sort(-log10(KromKosi_R_high_Cu_up_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,35),col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=KromKosi_R_high_Cu_up_fisher_sig$name[order(-KromKosi_R_high_Cu_up_fisher_sig$q_value)],at=bpKromKosi_R_high_Cu_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(Root),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
#mtext(text=expression("Transcript abundance higher in Krom than Kosi"),side=3,cex=1.3,outer=F,line=2,adj=0.3)

barplot(sort(-log10(KromKosi_S_high_Cu_up_fisher_sig$q_value)),beside=T,las=3,horiz=T,las=1,cex.lab=1.5,width=1.2,xlim=c(0,25),col="red",cex.axis=1.3)
axis(side=2,tick=F,labels=KromKosi_S_high_Cu_up_fisher_sig$name[order(-KromKosi_S_high_Cu_up_fisher_sig$q_value)],at=bpKromKosi_S_high_Cu_up_fisher_sig,xpd=T,srt=45,las=2,font=2,cex=1.3)
mtext(text=expression(-log[10](adjusted~p-value)),side=1,cex=1.3,outer=T,line=1)
mtext(text=expression(Shoot),side=3,cex=1.3,outer=F,line=0,at=-2,0.5)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend=c("Transcript abundance lower in Krom than Kosi","Transcript abundance higher in Krom than Kosi"),fill=c("blue","red"),xpd=T,bty="n",inset=c(0,0),cex=1.3)

dev.off()



