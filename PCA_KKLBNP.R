library("SNPRelate")
vcf.fn<-"AllCucombKKLBNP_4dg_sites.vcf"
snpgdsVCF2GDS(vcf.fn,"KKLBNP.gds",method="biallelic.only",ignore.chr.prefix=c("chr","SCF"))
genofile<-snpgdsOpen("KKLBNP.gds")
KKLBNP_pca<-snpgdsPCA(genofile,num.thread=8,algorithm="exact")
summary(KKLBNP_pca)
KKLBNP_pca$eigenval
jpeg("Eigenvalues_Kaiser-Guttman-Test.jpeg", width=26, height=18, units="cm", res=1000)
ev <- KKLBNP_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
barplot (ev, main ="Eigenvalues", col="black", las=2)
abline (h=mean(ev,na.rm=T), col="red")
dev.off()

pdf("Eigenvalues_Kaiser-Guttman-Test.pdf",width=8,height=6,paper="special",pointsize=16)
ev <- KKLBNP_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
barplot (ev, main ="Eigenvalues", col="black", las=2)
abline (h=mean(ev,na.rm=T), col="red")
dev.off()



ev <- KKLBNP_pca$eigenval
ev[ev> mean(ev,na.rm=T)]
n<- length (ev)
mean(ev,na.rm=T)
colour_pops<-c(rep("seagreen3",12),rep("coral",12),rep("coral4",12),rep("darkgreen",13),rep("darkblue",12),rep("dodgerblue",12))
jpeg("PCA_KKNPLB.jpeg", width=26, height=18, units="cm", res=1000)
plot(KKLBNP_pca$eigenvect[,1]~KKLBNP_pca$eigenvect[,2], pch=17,col=colour_pops,xlab=paste("PC1 (",round(KKLBNP_pca$varprop[1]*100,0),"%)"),ylab=paste("PC2 (",round(KKLBNP_pca$varprop[2]*100,0),"%)"))
legend("bottomleft",fill=unique(colour_pops)[c(2,3,1,4,6,5)],legend=c("Kosi","Krom","Best","Lang","Pais","Noss"))
#text(KKLBNP_pca$eigenvect[37,1]~KKLBNP_pca$eigenvect[37,2], labels=KKLBNP_pca$sample[37], cex= 0.7)
dev.off()

pdf("PCA_KKNPLB.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,6))
par(mgp=c(3.5,0.75,0))
plot(KKLBNP_pca$eigenvect[,1]~KKLBNP_pca$eigenvect[,2], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("topright",fill=unique(colour_pops)[c(3,2,NA,5,6,NA,4,1)],legend=c("Krom","Kosi",NA,"Noss","Pais",NA,"Lang","Best"),cex=1.3,border=c("black","black",NA),inset=c(-0.3,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(KKLBNP_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC2 (",round(KKLBNP_pca$varprop[2]*100,0),"%)"),cex=1.5)

dev.off()


pdf("Explained_var_PCA_GS.pdf",width=8,height=5,paper="special",pointsize=12)
Var <- KKLBNP_pca$varprop[!is.na(KKLBNP_pca$varprop)]*100
bp<-barplot (Var, ylab ="Variance explained [%]", col="black", las=2,las=1,cex.lab=1.3)
barplot (Var, ylab ="Variance explained [%]", col="black", las=2,las=1,cex.lab=1.3)
abline (h=mean(Var,na.rm=T), col="red")
text(cex=1, x=bp-0.75, y=-1, paste0("PC",1:32), xpd=TRUE, srt=45)
dev.off()


pdf("PCA_KKNPLB_1_3.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,6))
par(mgp=c(3.5,0.75,0))
plot(KKLBNP_pca$eigenvect[,1]~KKLBNP_pca$eigenvect[,3], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("topright",fill=unique(colour_pops)[c(3,2,NA,5,6,NA,4,1)],legend=c("Krom","Kosi",NA,"Noss","Pais",NA,"Lang","Best"),cex=1.3,border=c("black","black",NA),inset=c(-0.3,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(KKLBNP_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC3 (",round(KKLBNP_pca$varprop[3]*100,0),"%)"),cex=1.5)

dev.off()

pdf("PCA_KKNPLB_1_4.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,6))
par(mgp=c(3.5,0.75,0))
plot(KKLBNP_pca$eigenvect[,1]~KKLBNP_pca$eigenvect[,4], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("topright",fill=unique(colour_pops)[c(3,2,NA,5,6,NA,4,1)],legend=c("Krom","Kosi",NA,"Noss","Pais",NA,"Lang","Best"),cex=1.3,border=c("black","black",NA),inset=c(-0.3,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(KKLBNP_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC4 (",round(KKLBNP_pca$varprop[4]*100,0),"%)"),cex=1.5)

dev.off()

pdf("PCA_KKNPLB_1_5.pdf",width=8,height=6,paper="special",pointsize=16)
par(mar=c(5,5,1,6))
par(mgp=c(3.5,0.75,0))
plot(KKLBNP_pca$eigenvect[,1]~KKLBNP_pca$eigenvect[,5], pch=17,col=colour_pops,xlab="",ylab="",cex=1.3,cex.lab=1.5,las=1,cex.axis=1.3)
legend("topright",fill=unique(colour_pops)[c(3,2,NA,5,6,NA,4,1)],legend=c("Krom","Kosi",NA,"Noss","Pais",NA,"Lang","Best"),cex=1.3,border=c("black","black",NA),inset=c(-0.3,0),xpd=T,bty="n")
mtext(side=1,line=2.5,paste0("PC1 (",round(KKLBNP_pca$varprop[1]*100,0),"%)"),cex=1.5)
mtext(side=2,line=3.5,paste0("PC5 (",round(KKLBNP_pca$varprop[5]*100,0),"%)"),cex=1.5)

dev.off()


