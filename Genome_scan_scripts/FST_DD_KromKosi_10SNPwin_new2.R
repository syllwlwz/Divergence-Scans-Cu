################
#FST for Halleri#
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
d=d[(d$AC/d$AN)<max(d$AC/d$AN) | (d$AC.1/d$AN.1)<max(d$AC.1/d$AN.1),]
n1=d$AN #Krom
n2=d$AN.1 #Kosi
nt=n1+n2
winsize=10

#whole data set
d$atot=as.numeric(d[,3])+as.numeric(d[,4])
d$ptot=d$atot/(nt)
d=subset(d,d$ptot<1)
d$ht=2*d$ptot*(1-d$ptot)
d[,3]=as.numeric(d[,3])/n1
d[,4]=as.numeric(d[,4])/n2
d$h1=2*d[,3]*(1-d[,3])
d$h2=2*d[,4]*(1-d[,4])
d$h12=((d$h1*n1)+(d$h2*n2))/nt
FSTdata<-data.frame()
scaff<-data.frame()
#windows
for (i in levels(d[,1]))
	{datanow<-d[d[,1]==i,]
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
		wmatrix[i,3]=mean(abs(twin$ht-twin$h12)/twin$ht)
		wdt=wend+1
		wend=wend+winsize
		}
	wmatrix=cbind(scaff,wmatrix)
	FSTdata<-rbind(FSTdata,wmatrix)
	}
write.table(FSTdata,"FstKromKosiHalleri_10SNPwin_new2.csv",col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)


###############
#DD for Halleri#
###############

d<-data
names(d)<-c("CHROM","POS","AC","AC.1","AN","AN.1")
d=d[(d$AC/d$AN)<max(d$AC/d$AN) | (d$AC.1/d$AN.1)<max(d$AC.1/d$AN.1),]

wsize=10
cohort1_n=d$AN
cohort2_n=d$AN.1
d$AC[is.na(d$AC)]=0
d$AC.1[is.na(d$AC.1)]=0
d$cohort1_f=d$AC/cohort1_n
d$cohort2_f=d$AC.1/cohort2_n
d$diff=abs(d$cohort1_f-d$cohort2_f)
d$diffRaw=d$cohort1_f-d$cohort2_f
d$cohort1_pi=(2*d$cohort1_f*(1-d$cohort1_f))*(cohort1_n/(cohort1_n-1))
d$cohort2_pi=(2*d$cohort2_f*(1-d$cohort2_f))*(cohort2_n/(cohort2_n-1))

#DD from absolute values
scaff=data.frame()
DDdata<-data.frame()
for (i in levels(d[,1]))
        {datanow<-d[d[,1]==i,]
        w=matrix(data=0,ceiling(nrow(datanow)/wsize),9)
        wstart=1
        wend=wsize
        scaff<-i
        for(i in 1:nrow(w))
                {
                td=datanow[wstart:wend,]
                w[i,1]=min(td$POS)
                w[i,2]=max(td$POS)
                w[i,3]=w[i,1]+((w[i,2]-w[i,1])/2)
                w[i,4]=mean(td$cohort1_f,na.rm=T)
                w[i,5]=mean(td$cohort2_f,na.rm=T)
                w[i,6]=mean(td$diff,na.rm=T)
                w[i,7]=mean(td$cohort1_pi,na.rm=T)
                w[i,8]=mean(td$cohort2_pi,na.rm=T)
                w[i,9]=mean(td$diffRaw,na.rm=T)
                wstart=wend+1
                wend=wend+wsize
                }
        w=cbind(scaff,w)
        DDdata<-rbind(DDdata,w)
        }

DDdata<-na.exclude(DDdata)
DDdata[,8]<-as.numeric(as.character(DDdata[,8]))
DDdata[,7]<-as.numeric(as.character(DDdata[,7]))
lmdiversitytomeandiff<-lm(DDdata[,8]~DDdata[,7])
Intercept=coef(lmdiversitytomeandiff)[1]
Slope=coef(lmdiversitytomeandiff)[2]
b=DDdata[,7]-((DDdata[,8]-Intercept)/Slope)
a=Slope*DDdata[,7]+Intercept-DDdata[,8]
c=sqrt((a^2)+(b^2))
h=a*b/c
h_raw<-ifelse(DDdata[,7]<((DDdata[,8]-Intercept)/Slope)|DDdata[,8]>(Slope*DDdata[,7]+Intercept),h,-h)
dataDD<-cbind(DDdata,h_raw)
names(dataDD)<-c("Scaffold","Start_pos","End_pos","Mean_pos","Mean_allele_frequency_cohort1","Mean_allele_frequency_cohort2","Mean_abs_diff","Mean_pi_cohort1","Mean_pi_cohort2","Mean_raw_diff","DD")
write.table(dataDD,"DDKromKosiHalleri_10SNPwin_new2.csv",sep="\t",quote=F,row.names=F,col.names=F)


