#take numbers from Venn diagram

#phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
#where q=size of overlap-1; m=number of upregulated genes in experiment #1; 
#n=(total number of genes on platform-m); k=number of upregulated genes in experiment #2.


#‘pick m balls from a jar of n balls, then repeat by picking k balls: what is the significance of exactly q overlap. 
#genome scans and high effect var
q=4+1-1
m=511+4+4+1+45
k=29+4+5+1
n=32739-(511+4+4+1+45)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0005332951

#genome scans and CNVs
q=4-1
m=511+4+4+1+45
k=317+58+4
n=32739-(511+4+4+1+45)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.8943541

#genome scans and RNA-Seq
q=45+1-1
m=511+4+4+1+45
k=1675+58+5+1+45
n=32739-(511+4+4+1+45)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.004588276

#high effect var and RNA-Seq
q=5+1-1
m=29+4+5+1
k=1675+58+5+1+45
n=32739-(1675+58+5+1+45)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0230071

#high effect var and CNV
#0

#CNV and RNA-Seq
q=58-1
m=317+58+4
k=1675+58+5+1+45
n=32739-(1675+58+5+1+45)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#6.771165e-12


#############################################
#all expressed genes#
#############################################
All<-read.table("F://Cu_project//RNA-Seq//Mapman_RBH//Log_cpm.table",sep="\t",header=T)
#32739
All_LCu<-All[apply(All[,c(4:6,10:12,16:18,22:24)],1,function(x) any(x>1.2)),]
All_n<-nrow(All_LCu)
All_n
#21329

#genome scans and high effect var
q=2+1-1
m=353+2+1+1+41
k=14+2+1+5
n=21329-(353+2+1+1+41)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.0076

#genome scans and CNVs
#only 1

#genome scans and RNA-Seq
q=41+1-1
m=355+2+1+1+42
k=1472+47+5+1+41
n=21329-(353+2+1+1+41)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.013

#high effect var and RNA-Seq
q=5+1-1
m=14+2+1+5
k=1472+47+5+1+41
n=21329-(1472+47+5+1+41)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#0.006

#high effect var and CNV
#0

#CNV and RNA-Seq
q=47-1
m=73+47+1
k=1472+47+5+1+41
n=21329-(1472+47+5+1+41)

phyper(q, m, n, k, lower.tail = F, log.p = FALSE)
#2.1e-21







