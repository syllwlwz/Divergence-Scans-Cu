###########################
#Nielsen 2dSFS for halleri#
###########################
require("doParallel")
registerDoParallel(cores=12)
require("foreach")
dataA<-read.table("HalleriKromKosiGS_new.csv",sep="\t",header=TRUE)
names(dataA)<-c("CHROM","POS","AC","AN","CHROM.1","POS.1","AC.1","AN.1")
require("numbers")
AF_test1<-dataA$AC/dataA$AN
AF_test2<-dataA$AC.1/dataA$AN.1
AF_folded1<-ifelse(AF_test1>0.5,1-AF_test1,AF_test1)
AF_folded2<-ifelse(AF_test2>0.5,1-AF_test2,AF_test2)


AF1<-as.integer(AF_folded1*mLCM(unique(c(dataA$AN,dataA$AN.1))))
AF2<-as.integer(AF_folded2*mLCM(unique(c(dataA$AN,dataA$AN.1))))
data<-data.frame(dataA[,c(1:2)],AF1,dataA[,c(1:2)],AF2)
names(data)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")

pos<-2
f1<-3
f2<-6
##your window size!
wsize=10
d<-data
names(d)<-c("CHROM","POS","AC","CHROM.1","POS.1","AC.1")

#count occurrences of chromosomes
Chromcount<-as.data.frame(table(d[,1]))
Exclude<-Chromcount[Chromcount$Freq<10,]
Excludeonly<-droplevels(Exclude)
Excl<-1:nrow(Excludeonly)
for (i in 1:nrow(Excludeonly))
	{Excl[i]<-as.character(Excludeonly$Var1[i])
	}

Chrom<-1:nrow(d)
for (i in 1:nrow(d))
	{Chrom[i]<-as.character(d$CHROM[i])
	}
d$CHROM<-Chrom

for (i in 1:nrow(d))
	{if (d[i,1] %in% Excl)
		{d[i,1]<-"remove"
		}
	}

d<-d[d[,1]!="remove",]
d[,1]<-as.factor(d[,1])
d<-droplevels(d)



d$AC[is.na(d[,3])]=0
d$AC.1[is.na(d[,6])]=0
d$AC<-as.integer(d$AC)
d$AC.1<-as.integer(d$AC.1)

Nielsendata<-data.frame()
count=0
#windows
result<-foreach (scaff=levels(d[,1]),.combine=rbind) %dopar%
	{
	count=count+1
	datanow<-d[d[,1]==scaff,]
	##initialize your matrix for 2 dimensional sfs
	m=matrix(nrow=length(unique(datanow[,f1])),ncol=length(unique(datanow[,f2])),data=0)
	##enter values for 2dSFS - rows are pop 1 sfs,  columns are pop 2 sfs
	for(i in sort(unique(datanow[,f1]))) 
		{
		if (min(datanow[,f1])=="0")
			{
			k=which(sort(unique(datanow[,f1]))==i)
			for(j in sort(unique(datanow[,f2]))) 
				{
				if (min(datanow[,f2])=="0")
					{l=which(sort(unique(datanow[,f2]))==j)
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				else
					{l=which(sort(unique(datanow[,f2]))==j)
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				}
			}
		else
			{
			k=which(sort(unique(datanow[,f1]))==i)
			for(j in sort(unique(datanow[,f2]))) 
				{
				if (min(datanow[,f2])=="0")
					{l=which(sort(unique(datanow[,f2]))==j)
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				else
					{l=which(sort(unique(datanow[,f2]))==j)
					m[k,l]=m[k,l]+nrow(datanow[datanow[,f1]==i & datanow[,f2]==j,])
					}
				}
			}

		}
	m[nrow(m),ncol(m)]=0 ##presumably you wont consider sites fixed in both pops?
	mn=m
	## convert number of sites to proportions (= genome-wide probabilities)
	m=mn/sum(mn)

	##calculate number of windows in your dataframe
	wcl=floor(length(datanow[,1])/wsize)
	s=1
	e=wsize
	##initialize matrix to record window info
	stat=matrix(0,wcl,5)
	for(w in 1:nrow(stat)) 
		{
		tw=datanow[s:e,]
		mw=matrix(nrow=length(unique(datanow[,f1])),ncol=length(unique(datanow[,f2])),data=0) ##Window-Specific Matrix!
		stat[w,1]=1
		stat[w,2]=mean(tw[,pos])
		stat[w,3]=min(tw[,pos])
		stat[w,4]=max(tw[,pos])
		LnL0=0
		LnL1=0
		if (min(datanow[,f1])=="0")
			{for(i in sort(unique(datanow[,f1])))  ##Get 2dSFS for this window
				{
				k=which(sort(unique(datanow[,f1]))==i)
				if (min(datanow[,f2])=="0")
					{for(j in sort(unique(datanow[,f2]))) 
						{
						l=which(sort(unique(datanow[,f2]))==j)
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				else
					{for(j in sort(unique(datanow[,f2]))) 
						{
						l=which(sort(unique(datanow[,f2]))==j)
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				}
			}
		else
			{
			for(i in sort(unique(datanow[,f1])))
				{
				k=which(sort(unique(datanow[,f1]))==i)
				if (min(datanow[,f2])=="0")
					{for(j in sort(unique(datanow[,f2]))) 
						{
						l=which(sort(unique(datanow[,f2]))==j)
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				else
					{for(j in sort(unique(datanow[,f2]))) 
						{
						l=which(sort(unique(datanow[,f2]))==j)
						mw[k,l]=length(subset(tw[,f1],tw[,f1]==i & tw[,f2]==j))
						}
					}
				}
			}
		mw[nrow(mw),ncol(mw)]=0
		pw=mw/sum(mw)
		LnL1=log(prod(pw^mw))
		if(LnL1==(-Inf)) 
			{
			LnL1=3^-600
			}
		LnL0=log(prod(m^mw))
		if(LnL0==(-Inf)) 
			{
			LnL0=3^-600
			}
		stat[w,5]=2*(LnL1-LnL0) ##your G-test statistic
#		print(stat[w,5])
		s=e+1
		e=e+wsize
		}
	stat=cbind(scaff,stat)
	Nielsendata<-rbind(Nielsendata,stat)
	Nielsendata
	}

write.table(result,"NielsenKromKosi_folded_parallel2.csv",quote=F,row.names=F,col.names=F)

