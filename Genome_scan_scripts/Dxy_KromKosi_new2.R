################
#Dxy for Halleri#
################
data<-read.table("HalleriKromKosiGS_new.csv",sep="\t",header=TRUE)
names(data)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1")
data<-data.frame(data$CHROM,data$POS,data$AC,data$AC.1,data$AN,data$AN.1)

n1=data$data.AN #Krom
n2=data$data.AN.1 #Kosi
nt=n1+n2
winsize=10
d<-data
names(d)<-c("CHROM","POS","AC","AC.1","AN","AN.1")
d=d[d$AC<max(d$AC) | d$AC.1<max(d$AC.1),]
#whole data set

#calculation of dxy
AF1=data[,3]/n1
AF2=data[,4]/n2
dxy=AF1*(1-AF2)+(AF2*(1-AF1))
dxy[1:3]
#windows
Dxy_data<-data.frame(data[,1:4],dxy)
names(Dxy_data)<-c("CHROM","POS","AC","AC.1","Dxy")
Dxy_data[1:3,]

Dxydata=data.frame()
for (i in levels(Dxy_data[,1]))
	{datanow<-Dxy_data[Dxy_data[,1]==i,]
	nwins=ceiling(nrow(datanow)/winsize)
	wmatrix=matrix(nrow=nwins,ncol=3,data=0)
	wdt=1
	wend=winsize
	scaff<-i
	for(i in 1:nrow(wmatrix)) 
		{
		twin=datanow[wdt:wend,]
		wmatrix[i,1]=min(twin[,2])
		wmatrix[i,2]=max(twin[,2])
		wmatrix[i,3]=sum(twin[,5])/winsize
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	Dxydata<-rbind(Dxydata,wmatrix)
	}
write.table(Dxydata,"DxyKromKosi_new2.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
