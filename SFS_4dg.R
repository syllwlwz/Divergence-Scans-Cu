Intervals<-read.table("Ahalleri_Felix_4dg_sites.intervals.list",sep="\t")


AFdata_KK1<-read.table("HalleriKromKosiGS_new.csv",header=TRUE,sep="\t")
AFdata_KK<-AFdata_KK1[paste(AFdata_KK1$CHROM,AFdata_KK1$POS,sep=":")%in%Intervals$V1,]
AFdata_KK$AF<-AFdata_KK$AC/AFdata_KK$AN
AFdata_KK$AF.1<-AFdata_KK$AC.1/AFdata_KK$AN.1
names(AFdata_KK)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_KK=AFdata_KK[AFdata_KK$AF<max(AFdata_KK$AF) | AFdata_KK$AF.1<max(AFdata_KK$AF.1),]
d_KK_Krom<-ifelse(d_KK$AF>0.5,1-d_KK$AF,d_KK$AF)
d_KK_Kosi<-ifelse(d_KK$AF.1>0.5,1-d_KK$AF.1,d_KK$AF.1)

AFdata_NP1<-read.table("HalleriNossPaisGS_new.csv",header=TRUE,sep="\t")
AFdata_NP<-AFdata_NP1[paste(AFdata_NP1$CHROM,AFdata_NP1$POS,sep=":")%in%Intervals$V1,]
AFdata_NP$AF<-AFdata_NP$AC/AFdata_NP$AN
AFdata_NP$AF.1<-AFdata_NP$AC.1/AFdata_NP$AN.1
names(AFdata_NP)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_NP=AFdata_NP[AFdata_NP$AF<max(AFdata_NP$AF) | AFdata_NP$AF.1<max(AFdata_NP$AF.1),]
d_NP_Noss<-ifelse(d_NP$AF>0.5,1-d_NP$AF,d_NP$AF)
d_NP_Pais<-ifelse(d_NP$AF.1>0.5,1-d_NP$AF.1,d_NP$AF.1)

AFdata_LB1<-read.table("HalleriLangBestGS_new.csv",header=TRUE,sep="\t")
AFdata_LB<-AFdata_LB1[paste(AFdata_LB1$CHROM,AFdata_LB1$POS,sep=":")%in%Intervals$V1,]
AFdata_LB$AF<-AFdata_LB$AC/AFdata_LB$AN
AFdata_LB$AF.1<-AFdata_LB$AC.1/AFdata_LB$AN.1
names(AFdata_LB)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_LB=AFdata_LB[AFdata_LB$AF<max(AFdata_LB$AF) | AFdata_LB$AF.1<max(AFdata_LB$AF.1),]
d_LB_Lang<-ifelse(d_LB$AF>0.5,1-d_LB$AF,d_LB$AF)
d_LB_Best<-ifelse(d_LB$AF.1>0.5,1-d_LB$AF.1,d_LB$AF.1)

AFdata_WH<-read.table("HalleriWulmHDGS_new.csv",header=TRUE,sep="\t")
AFdata_WH$AF<-AFdata_WH$AC/AFdata_WH$AN
AFdata_WH$AF.1<-AFdata_WH$AC.1/AFdata_WH$AN.1
names(AFdata_WH)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1","AF","AF.1")
d_WH=AFdata_WH[AFdata_WH$AF<max(AFdata_WH$AF) | AFdata_WH$AF.1<max(AFdata_WH$AF.1),]
d_WH_Wulm<-ifelse(d_WH$AF>0.5,1-d_WH$AF,d_WH$AF)
d_WH_HD<-ifelse(d_WH$AF.1>0.5,1-d_WH$AF.1,d_WH$AF.1)




require(DescTools)
pos<-barplot(t(cbind(Freq(d_KK_Kosi,breaks=10)$perc,Freq(d_KK_Krom,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
pdf("SFS_KKNPLB_ordered_4dg.pdf", width=12, height=4, paper="special",pointsize=15)
par(mfrow=c(1,3))
par(mar=c(1,1,1,0)+0.1)
par(oma=c(4,5,0,1))
par(mgp=c(5,0.2,0))
options(scipen=10)
par(lwd=2)
barplot(t(cbind(Freq(d_KK_Kosi,breaks=10)$perc,Freq(d_KK_Krom,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
legend("topright",legend=c("Krom","Kosi"),pch=15,fill="white",cex=1.5,bty="n",border=c("red","black"),inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.3,cex.lab=1.3,mgp=c(5,0.75,0),las=1,at=apply(pos,2,mean)[seq(2,10,2)],labels=c(seq(0.1,0.5,0.1)))
axis(2,line=0,cex.axis=1.3,cex.lab=1.3,mgp=c(5,0.75,0),las=1)
box(lwd=1)

barplot(t(cbind(Freq(d_NP_Pais,breaks=10)$perc,Freq(d_NP_Noss,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
legend("topright",legend=c("Noss","Pais"),pch=15,fill="white",cex=1.5,bty="n",border=c("red","black"),inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.3,cex.lab=1.3,mgp=c(5,0.75,0),las=1,at=apply(pos,2,mean)[seq(2,10,2)],labels=c(seq(0.1,0.5,0.1)))
box(lwd=1)

barplot(t(cbind(Freq(d_LB_Best,breaks=10)$perc,Freq(d_LB_Lang,breaks=10)$perc)),main="",cex=1,cex.lab=1.3,cex.axis=1,mgp=c(2,1,0),axes=F,beside=T,border=c("black","red"),col="white",space=c(0.1,1),ylim=c(0,0.5),las=1)
legend("topright",legend=c("Lang","Best"),pch=15,fill="white",cex=1.5,bty="n",border=c("red","black"),inset=0.05,col=NA)
axis(1,line=0,cex.axis=1.3,cex.lab=1.3,mgp=c(5,0.75,0),las=1,at=apply(pos,2,mean)[seq(2,10,2)],labels=c(seq(0.1,0.5,0.1)))
box(lwd=1)

mtext(side=1,line=2,"Allele frequency",outer=T,cex=1.2)
mtext(side=2,line=2,"Proportion of SNPs [%]",outer=T,cex=1.2)

dev.off()


