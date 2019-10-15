options(java.parameters = "-Xmx12g")
require(xlsx)

####################################
#Folded AF plots with Krom and Kosi#
####################################
AFdata<-read.table("HalleriKromKosiGS_new.csv",header=TRUE,sep="\t")
test<-as.character(AFdata$CHROM)
AFdata[,1]<-test
AFdata$AF<-AFdata$AC/AFdata$AN
AFdata$AF.1<-AFdata$AC.1/AFdata$AN.1
AFdata2<-data.frame(AFdata[,c(1,2,9)],AFdata[,10])
names(AFdata2)<-c("Chrom","POS","AF_Krom","AF_Kosi")
AFdata<-AFdata2
##########################
#Candidate genes Fst 0.1%#
##########################

RNASeqR<-read.xlsx2("E://Cu_project//RNA-Seq//Sample_counts//Deseq2_medofratios_stricter_Mapmanadded.xlsx",6,header=T)
RNASeqS<-read.xlsx2("E://Cu_project//RNA-Seq//Sample_counts//Deseq2_medofratios_stricter_Mapmanadded.xlsx",7,header=T)

Candidates2<-read.xlsx2("Genes_KromKosi_dedup.xlsx",1,header=TRUE)

Candidates5<-Candidates2[!duplicated(Candidates2$Ath_ID),]
Candidates3<-Candidates5[!is.na(Candidates5$Ath_ID),]
Candidates4<-Candidates3[!Candidates3$Ath_ID=="",]
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0006377",]
Candidates$gene_end<-1117289#two transcripts

Candidates<-Candidates4[Candidates4$Ah_ID%in%RNASeqR$Ah_ID|Candidates4$Ah_ID%in%RNASeqS$Ah_ID,]
#31
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0004954"|Candidates4$Ah_ID=="AHAL_G0003865"|Candidates4$Ah_ID=="AHAL_G0007270"|Candidates4$Ah_ID=="AHAL_G0006303"|Candidates4$Ah_ID=="AHAL_G0006286"|Candidates4$Ah_ID=="AHAL_G0010056"|Candidates4$Ah_ID=="AHAL_G0025694"|Candidates4$Ah_ID=="AHAL_G0030993"|Candidates4$Ah_ID=="AHAL_G0007270"|Candidates4$Ah_ID=="AHAL_G0025599"|Candidates4$Ah_ID=="AHAL_G0012884"|Candidates4$Ah_ID=="AHAL_G0001510"|Candidates4$Ah_ID=="AHAL_G0014338"|Candidates4$Ah_ID=="AHAL_G0014700"|Candidates4$Ah_ID=="AHAL_G0015690",]
test2<-as.character(Candidates$Ath_ID)
Candidates$Ath_ID<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0004954"|Candidates4$Ah_ID=="AHAL_G0003865"|Candidates4$Ah_ID=="AHAL_G0007270"|Candidates4$Ah_ID=="AHAL_G0006303"|Candidates4$Ah_ID=="AHAL_G0006286"|Candidates4$Ah_ID=="AHAL_G0010056"|Candidates4$Ah_ID=="AHAL_G0025694"|Candidates4$Ah_ID=="AHAL_G0030993"|Candidates4$Ah_ID=="AHAL_G0007270"|Candidates4$Ah_ID=="AHAL_G0025599"|Candidates4$Ah_ID=="AHAL_G0012884"|Candidates4$Ah_ID=="AHAL_G0001510"|Candidates4$Ah_ID=="AHAL_G0014338"|Candidates4$Ah_ID=="AHAL_G0014700"|Candidates4$Ah_ID=="AHAL_G0015690"|Candidates4$Ah_ID=="AHAL_G0002108"
|Candidates4$Ah_ID=="AHAL_G0024590"|Candidates4$Ah_ID=="AHAL_G0008885"|Candidates4$Ah_ID=="AHAL_G0024044"|Candidates4$Ah_ID=="AHAL_G0015076"|Candidates4$Ah_ID=="AHAL_G0012967"|Candidates4$Ah_ID=="AHAL_G0012903"|Candidates4$Ah_ID=="AHAL_G0012433"|Candidates4$Ah_ID=="AHAL_G0012418"|Candidates4$Ah_ID=="AHAL_G0011147"|Candidates4$Ah_ID=="AHAL_G0008318"|Candidates4$Ah_ID=="AHAL_G0017626"|Candidates4$Ah_ID=="AHAL_G0017746"|Candidates4$Ah_ID=="AHAL_G0024590"|Candidates4$Ah_ID=="AHAL_G0028807"|Candidates4$Ah_ID=="AHAL_G0021363"|Candidates4$Ah_ID=="AHAL_G0030810"|Candidates4$Ah_ID=="AHAL_G0030485"|Candidates4$Ah_ID=="AHAL_G0030123"|Candidates4$Ah_ID=="AHAL_G0029902"|Candidates4$Ah_ID=="AHAL_G0006377",]
test2<-as.character(Candidates$Ath_ID)
Candidates$Ath_ID<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

#convergent candidates
Candidates<-Candidates4[Candidates4$Ah_ID=="AHAL_G0012418"|Candidates4$Ah_ID=="AHAL_G0015692"|Candidates4$Ah_ID=="AHAL_G0021044"|Candidates4$Ah_ID=="AHAL_G0025603"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0030401"|Candidates4$Ah_ID=="AHAL_G0009828"|Candidates4$Ah_ID=="AHAL_G0027258"|Candidates4$Ah_ID=="AHAL_G0010679"|Candidates4$Ah_ID=="AHAL_G0013166"|Candidates4$Ah_ID=="AHAL_G0024026",]
test2<-as.character(Candidates$Ath_ID)
Candidates$Ath_ID<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)


test2<-as.character(Candidates$Ath_ID)
Candidates$Ath_ID<-test2
Candidates$gene_start<-as.numeric(as.character(Candidates$gene_start))
Candidates$gene_end<-as.numeric(as.character(Candidates$gene_end))
Candidates$gene_size<-as.numeric(as.character(Candidates$gene_size))
Candidates$Chr<-as.character(Candidates$Chr)

genemiddle=list((Candidates$gene_start+Candidates$gene_end)/2)
windowsize<-Candidates$gene_size*1.5 #bp genemiddle +/- windowsize will be displayed

for (j in 1:nrow(Candidates))
	{nam <- paste("AFofinterest", j, sep ="")
 	AF1<-AFdata[AFdata$Chrom==Candidates$Chr[j],]
	AFofinterest<-AF1[AF1$POS>=(genemiddle[[1]][j]-windowsize[j])&AF1$POS<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, AFofinterest)
	}
rm(AFofinterest)
list<-ls(pattern="^AFofinterest")


genelist<-read.table("AHAL_geneid.gff3",sep="\t",header=T)
colnames(genelist)<-c("Scaffold","Source","Type","Start_pos","End_pos","Score","Strand","Phase","Attributes","Ah_ID")
genelist<-genelist[genelist$Type=="gene",]
genelist<-droplevels(genelist)
genelist$Scaffold<-as.character(genelist$Scaffold)
genelist<-genelist[!duplicated(genelist[,-9]),]

for (j in 1:nrow(Candidates))
	{nam <- paste("Ath_IDsofinterest", j, sep ="")
 	Ath_IDs1<-genelist[genelist$Scaffold==Candidates$Chr[j],]
	Ath_IDsofinterest<-Ath_IDs1[Ath_IDs1$Start_pos>=(genemiddle[[1]][j]-windowsize[j])&Ath_IDs1$Start_pos<=(genemiddle[[1]][j]+windowsize[j])|Ath_IDs1$End_pos>=(genemiddle[[1]][j]-windowsize[j])&Ath_IDs1$End_pos<=(genemiddle[[1]][j]+windowsize[j]),]
	assign(nam, Ath_IDsofinterest)
	}
rm(Ath_IDsofinterest)

arrowdir2<-c(rep(1,nrow(Candidates)))
arrowdir2<-ifelse(Candidates$gene_strand=="+",2,1)

test<-as.factor(AFdata$Chrom)
AFdata[,1]<-test

winsize=10
AFmatrix=data.frame()
for (i in levels(AFdata[,1]))
        {datanow<-AFdata[AFdata[,1]==i,]
        nwins=ceiling(nrow(datanow)/winsize)
        wmatrix=matrix(nrow=nwins,ncol=5,data=0)
        wdt=1
        wend=winsize
        scaff<-i
        for(i in 1:nrow(wmatrix)) 
                {
                twin=datanow[wdt:wend,]
                wmatrix[i,1]=min(twin[,2])
                wmatrix[i,2]=max(twin[,2])
		    wmatrix[i,3]=wmatrix[i,1]+((wmatrix[i,2]-wmatrix[i,1])/2)
                wmatrix[i,4]=mean(twin[,3])
                wmatrix[i,5]=mean(twin[,4])
                wdt=wend+1
                wend=wend+winsize
                }
        wmatrix=cbind(scaff,wmatrix)
        AFmatrix=rbind(AFmatrix,wmatrix)
        }
names(AFmatrix)<-c("Chrom","Window_start","Window_end","Window_middle","AF_Krom","AF_Kosi")


test<-as.character(AFdata$Chrom)
AFdata[,1]<-test
test<-as.character(AFmatrix$Chrom)
AFmatrix[,1]<-test

test2<-AFmatrix
for (i in 2:6)
	{test2[,i]<-as.numeric(as.character(AFmatrix[,i]))
	}
AFmatrix<-test2


for (j in 1:nrow(Candidates))
	{nam <- paste("AFfreqofinterest", j, sep ="")
 	AF1<-AFmatrix[AFmatrix$Chrom==Candidates$Chr[j],]
	AFfreqofinterest<-AF1[AF1$Window_start>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_start<=(genemiddle[[1]][j]+windowsize[j]+100000)|AF1$Window_end>=(genemiddle[[1]][j]-windowsize[j]-100000)&AF1$Window_end<=(genemiddle[[1]][j]+windowsize[j]+100000),]
	assign(nam, AFfreqofinterest)
	}
rm(AFfreqofinterest)

SNPeffdf<-read.table("SNPeff_KromKosi.table",sep="\t",header=T)

AFdatasnp<-cbind(AFdata[,1:3],AFdata[,4])
names(AFdatasnp)<-c("CHROM","POS","AC","AC.1")
SNPeffdfAF<-merge(SNPeffdf,AFdatasnp)
SNPeffdfAFo<-SNPeffdfAF[SNPeffdfAF$Effect=="MODERATE",]
SNPeffdfAFp<-SNPeffdfAF[SNPeffdfAF$Effect=="MODIFIER",]
SNPeffdfAFr<-SNPeffdfAF[SNPeffdfAF$Effect=="HIGH",]

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
 	Exons1<-na.exclude(exonlist[exonlist$Scaffold==Candidates$Chr[j],])
	Exonsofinterestcand<-Exons1[Exons1$ID==Candidates$Ah_ID[j],]
	assign(nam, Exonsofinterestcand)
	}
rm(Exonsofinterestcand)

for (j in 1:nrow(Candidates))
	{
	arrowdir<-ifelse(get(paste("Ath_IDsofinterest",j,sep=""))$Strand=="+",2,1)
	nam <- paste("GS_AFKromKosi_absoluteAFplots_candidate_genes_",j,"_",Candidates$Ah_ID[j], sep ="")
	pdf(paste(nam, '.pdf', sep = ''), width=7, height=6,paper="special")
	par(mar=c(5,5,1,1)+0.1)
	par(oma=c(0,0,1,10)+0.1)
	par(mgp=c(3,0.7,0))
	plot(abs(get(paste("AFofinterest",j,sep=""))$AF_Krom-get(paste("AFofinterest",j,sep=""))$AF_Kosi)~get(paste("AFofinterest",j,sep=""))$POS,ylab="",ylim=c(0,1.2),cex.lab=1,xlab="",xlim=c(genemiddle[[1]][j]-windowsize[j],genemiddle[[1]][j]+windowsize[j]),xaxt='n',yaxt='n')
	rect(Candidates$gene_start[j],-2, Candidates$gene_end[j], 2,col="grey",border = NA)	
	points(abs(get(paste("AFofinterest",j,sep=""))$AF_Krom-get(paste("AFofinterest",j,sep=""))$AF_Kosi)~get(paste("AFofinterest",j,sep=""))$POS)
	text(x=genemiddle[[1]][j],y=1.15,labels=Candidates$Ah_ID[j],col="red")
	for(i in 1:nrow(get(paste("Ath_IDsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Ath_IDsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Ath_IDsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir[i],length=0.1,col="grey",lwd=1)
		}
	arrows(x0=Candidates$gene_start[j],y0=1.1,x1=Candidates$gene_end[j],y1=1.1,code=arrowdir2[j],length=0.1,col="red",lwd=1)
	for(i in 1:nrow(get(paste("Exonsofinterest",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterest",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterest",j,sep=""))$End_pos[i],y1=1.1,code=1,length=0,col="grey",lwd=3)
		}
	arrowdir3<-ifelse(get(paste("Exonsofinterestcand",j,sep=""))$Strand=="+",2,1)
	for(i in 1:nrow(get(paste("Exonsofinterestcand",j,sep="")))) 
		{arrows(x0=get(paste("Exonsofinterestcand",j,sep=""))$Start_pos[i],y0=1.1,x1=get(paste("Exonsofinterestcand",j,sep=""))$End_pos[i],y1=1.1,code=arrowdir3[i],length=0,col="red",lwd=3)
		}
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Krom[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="darkblue")
	lines(get(paste("AFfreqofinterest",j,sep=""))$AF_Kosi[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)]~
	get(paste("AFfreqofinterest",j,sep=""))$Window_middle[order(get(paste("AFfreqofinterest",j,sep=""))$Window_middle)],lwd=2,col="lightblue")
	legend("topright",legend=c(expression(10~SNP-window~bold(CS[Kr])),expression(10~SNP-window~bold(RS[Ko]))),pch=15,col=c("darkblue","lightblue"),bg="white",bty="n",xpd=NA,inset=c(-0.59,0.02),title="Average allele frequency")
	points(abs(SNPeffdfAFp$AC[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Ah_ID[j])]-SNPeffdfAFp$AC.1[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Ah_ID[j])])~SNPeffdfAFp$POS[as.character(SNPeffdfAFp$Gene)==as.character(Candidates$Ah_ID[j])],col="purple")
	points(abs(SNPeffdfAFo$AC[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Ah_ID[j])]-SNPeffdfAFo$AC.1[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Ah_ID[j])])~SNPeffdfAFo$POS[as.character(SNPeffdfAFo$Gene)==as.character(Candidates$Ah_ID[j])],col="orange")
	points(abs(SNPeffdfAFr$AC[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Ah_ID[j])]-SNPeffdfAFr$AC.1[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Ah_ID[j])])~SNPeffdfAFr$POS[as.character(SNPeffdfAFr$Gene)==as.character(Candidates$Ah_ID[j])],col="red")
	legend("topright",legend=c("low","moderate","modifier","high"),col=c("black","orange","purple","red"),pch=1,pt.lwd=2,title="Predicted effect",bty="n",xpd=NA,inset=c(-0.35,0.2),title.adj=2)
	axis(side=1,cex=1,at=axTicks(side=1),labels=formatC(axTicks(side=1)/1000,big.mark=",",format="d"),cex.axis=0.8,adj=0.1)
	axis(side=2,cex=1,at=axTicks(side=2),labels=axTicks(side=2),las=2,cex.axis=0.8)
	box()
	mtext(text="Scaffold position (kb)",side=1,line=2,cex=1)
	mtext(text=expression(Absolute~allele~frequency~difference~bold(CS[Kr])~-~bold(RS[Ko])),side=2,line=2,cex=1)
	dev.off()
	}



