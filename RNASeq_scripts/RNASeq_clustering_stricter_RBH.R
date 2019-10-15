options(java.parameters = "-Xmx12000m")
require(xlsx)

Sample10<-read.table("Sample10B_counts.txt",sep="\t")
Sample11<-read.table("Sample11_counts.txt",sep="\t")
Sample1<-read.table("Sample1_counts.txt",sep="\t")
Sample2<-read.table("Sample2B_counts.txt",sep="\t")
Sample3<-read.table("Sample3_counts.txt",sep="\t")
Sample4<-read.table("Sample4_counts.txt",sep="\t")
Sample5<-read.table("Sample5_counts.txt",sep="\t")
Sample6<-read.table("Sample6_counts.txt",sep="\t")
Sample7<-read.table("Sample7_counts.txt",sep="\t")
Sample8<-read.table("Sample8_counts.txt",sep="\t")
Sample9<-read.table("Sample9B_counts.txt",sep="\t")
Sample12<-read.table("Sample12_counts.txt",sep="\t")
Sample13<-read.table("Sample13_counts.txt",sep="\t")
Sample14<-read.table("Sample14B_counts.txt",sep="\t")
Sample15<-read.table("Sample15B_counts.txt",sep="\t")
Sample16<-read.table("Sample16_counts.txt",sep="\t")
Sample17<-read.table("Sample17B_counts.txt",sep="\t")
Sample18<-read.table("Sample18B_counts.txt",sep="\t")
Sample19<-read.table("Sample19_counts.txt",sep="\t")
Sample20<-read.table("Sample20_counts.txt",sep="\t")
Sample21<-read.table("Sample21_counts.txt",sep="\t")
Sample22<-read.table("Sample22_counts.txt",sep="\t")
Sample23<-read.table("Sample23_counts.txt",sep="\t")
Sample24<-read.table("Sample24_counts.txt",sep="\t")

All_samples<-data.frame(Sample1[,1],Sample1[,2],Sample2[,2],Sample3[,2],Sample4[,2],Sample5[,2],Sample6[,2],Sample7[,2],Sample8[,2],Sample9[,2],Sample10[,2],Sample11[,2],Sample12[,2],Sample13[,2],Sample14[,2],Sample15[,2],Sample16[,2],Sample17[,2],Sample18[,2],Sample19[,2],Sample20[,2],Sample21[,2],Sample22[,2],Sample23[,2],Sample24[,2])
colnames(All_samples)<-c("Gene","Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8","Sample9","Sample10","Sample11","Sample12","Sample13","Sample14","Sample15","Sample16","Sample17","Sample18","Sample19","Sample20","Sample21","Sample22","Sample23","Sample24")

filelist = list.files(pattern = "*.txt") # Make a file list from all count text files
require(gtools)
filelist <- mixedsort(filelist) # Sort the numbers: 1 to 72

datafr = do.call(cbind,lapply(filelist,function(fn)read.table(fn,header=FALSE, sep="\t")[,2])) #merge all count files
genes <- read.delim("Sample1_counts.txt", head = FALSE)[,1] # for the first column
genes2 <- as.character(genes)

count0 = cbind(genes2, datafr) #merge with gene name
countdata <- count0[,-1]
rownames(countdata) <- count0[,1] 
colnames(countdata) <- c(1:24) # add column names
colnames(countdata)[1:24] <- paste("S", colnames(countdata[,c(1:24)]), sep = "") # add a letter P to the column names just to avoid confusion
mode(countdata)
mode(countdata) = "numeric" #convert the class of a matrix from character to numeric
class(countdata)

statsPerSample <- data.frame(t(apply(countdata, 2, summary)))
head(statsPerSample)

setwd("E:/RNASeq_Cu")
# Load the pre-made coldata to represent experimental designs
#coldata <- read.delim("RNA_samples_names.csv", header = TRUE,sep=";")
#rownames(coldata) <- coldata$ID

colnames(countdata) <- coldata$ID
countdata_df<-data.frame(row.names(countdata),as.data.frame(countdata))
colnames(countdata_df)<-c("ID",paste("S",1:24, sep = ""))
## Calculation of FPKM

# 1) Load gff and Merge
require(Biostrings)
CDS=readDNAStringSet("AHAL.cds.longest.renamed.fasta")
seq_name = names(CDS)
sequence = paste(CDS)
CDS_df <- data.frame(seq_name, sequence)
CDS_df$Tlen<-nchar(as.character(CDS_df$sequence))
PS_Hyp_genes<-read.table("Pseudogenes_hypothetical_Tlen.table",sep="\t",header=T) 
        
FPKM_df1 = merge(CDS_df[,c(1,3)],countdata_df,by.x="seq_name",by.y="ID",all.y=T)
FPKM_df1$Tlen[is.na(FPKM_df1$Tlen)]<-PS_Hyp_genes$Tlen[match(FPKM_df1$seq_name[is.na(FPKM_df1$Tlen)],PS_Hyp_genes$V10)]
FPKM_df<-FPKM_df1
head(FPKM_df)


# 2) Calculate
  # FPKM=(Counts/(gene size/1000))/(sum mapped fragments/1 million)
for (i in 1:24)
	{Si<-data.frame((FPKM_df[,2+i]/(FPKM_df$Tlen/1000))/(sum(FPKM_df[,2+i])/1000000))
	colnames(Si)<-paste("FPKM_S",i,sep="")
	FPKM_df<-data.frame(FPKM_df,Si)
	}
write.table(as.data.frame(FPKM_df),"FPKM_df.txt",sep="\t",row.names=FALSE)

# 3) Plot FPKM (vs counts)
plot(log2(FPKM_df$S1+1), log2(FPKM_df$FPKM_S1+1), xlab = "log2(FPKM+1) ", ylab = "log2(count+1)", cex = .7, main = "S1",xlim=c(0,20),ylim=c(0,20))
plot(log2(FPKM_df$S2+1), log2(FPKM_df$FPKM_S2+1), xlab = "log2(FPKM+1) ", ylab = "log2(count+1)", cex = .7, main = "S2",xlim=c(0,20),ylim=c(0,20))


# 4) Explore
# Create a histogram of expression (FPKM) using a log-transformation.
D_S1 = density((log10(FPKM_df$FPKM_S1+1)))
plot(D_S1, main = "FPKM_S1", xlab = "log2(FPKM + 1)")

D_S2 = density((log10(FPKM_df$FPKM_S2+1)))
plot(D_S2, main = "FPKM_S2", xlab = "log2(FPKM + 1)")

hist(log2(FPKM_df$FPKM_S1+1), main = "FPKM_S1", xlab = "log2(FPKM + 1)")

hist(FPKM_df$Tlen, breaks=50)

plot(FPKM_df$Tlen,log2(FPKM_df$FPKM_S2+1), xlab ="Transcript length [bp]", ylab = "log2(FPKM+1)", cex = .7, main = "S2")
plot(FPKM_df$Tlen,log2(FPKM_df$FPKM_S1+1), xlab ="Transcript length [bp]", ylab = "log2(FPKM+1)", cex = .7, main = "S1")

#hierarchical clustering
require("vegan")
FPKM_stand<-decostand(log2(na.omit(FPKM_df[,27:50])+1),"standardize")
colnames(FPKM_stand)<-coldata$Sample.Name
FPKM_dist<-dist(t(FPKM_stand))

hcluster<-hclust(FPKM_dist,method="ward.D")
pdf("Hierarchical_clustering_Cu_stand_FPKM.pdf",width=4.72441,height=4.72441,paper="special")
plot(hcluster)
dev.off()


#Make a data frame form a matrix
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
require(DESeq2)
Data_info<-strsplit(as.character(coldata$Sample.Name)," ")
Cu_condition1<-lapply(Data_info,"[",2)
Cu_condition<-data.frame(matrix(unlist(Cu_condition1),nrow=length(Cu_condition1),byrow=T))
Cu_Population1<-lapply(Data_info,"[",1)
Cu_Population<-data.frame(matrix(unlist(Cu_Population1),nrow=length(Cu_Population1),byrow=T))
levels(Cu_condition[,1])[levels(Cu_condition[,1])=="+Cu"]<-"high_Cu"
levels(Cu_condition[,1])[levels(Cu_condition[,1])=="-Cu"]<-"normal_Cu"
Cu_Tissue1<-lapply(Data_info,"[",3)
Cu_Tissue2<-data.frame(matrix(unlist(Cu_Tissue1),nrow=length(Cu_Tissue1),byrow=T))
Cu_Tissue<-substr(as.character(Cu_Tissue2[,1]),1,1)

coldata$Cu_condition<-Cu_condition[,1]
coldata$Population<-Cu_Population[,1]
coldata$Tissue<-as.factor(Cu_Tissue)


dds <- DESeqDataSetFromMatrix(countData=countdata,colData=coldata,design=~Tissue+Population+Cu_condition)

## Check the general properties of the DESeq dataset
print(dds)
dds_root<-dds[,dds$Tissue=="R"]
dds_shoot<-dds[,dds$Tissue=="S"]
dds_root$Tissue<-droplevels(dds_root$Tissue)
dds_shoot$Tissue<-droplevels(dds_shoot$Tissue)
colData(dds_root)<-colData(dds_root)[,c(1:4)]
colData(dds_shoot)<-colData(dds_shoot)[,c(1:4)]
design(dds_root)<-~Population+Cu_condition+Population:Cu_condition
design(dds_shoot)<-~Population+Cu_condition+Population:Cu_condition

#### Normalization
# Normalizing for different numbers of aligned reads per library  
dds_root.norm <-  estimateSizeFactors(dds_root)
sizeFactors(dds_root.norm)
#1.1166632 1.0562173 0.8989707 1.0453721 0.9306971 0.8460765 0.8579300
#0.8134182 0.9880103 1.3013973 1.1267492 1.1794839
data.frame(dds_root$Sample.Name,sizeFactors(dds_root.norm))


dds_shoot.norm <-  estimateSizeFactors(dds_shoot)
sizeFactors(dds_shoot.norm)
#1.0188702 0.9085191 0.9786146 1.1507284 0.8299584 1.0014517 1.0884745
#0.9453697 1.0156135 1.1536615 0.9359362 1.0737975
data.frame(dds_shoot$Sample.Name,sizeFactors(dds_shoot.norm))

require(affy)
# Checking the normalization
jpeg("Distributions_raw_counts_and_normalized_counts_root.jpeg",width =8.447917,height=6.072917,units="in",res=1000)
par(mfrow=c(2,2),mar=c(2,7,1,0.5),cex.lab=0.7)#set parameters for the plotting window
epsilon <- 1 # pseudo-count to avoid problems with log(0)
boxplot(log2(counts(dds_root.norm)+epsilon),col=c("red","red","red","grey","grey","grey"), cex.axis=0.7, 
        las=1, xlab="log2(counts+1)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds_root.norm, normalized=TRUE)+epsilon), col=c("red","red","red","grey","grey","grey"),cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds_root.norm)+epsilon), 
            xlab="log2(counts+1)", col=c("red","red","red","grey","grey","grey"),cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds_root.norm, normalized=TRUE)+epsilon), 
            xlab="log2(normalized counts)",col=c("red","red","red","grey","grey","grey"), cex.lab=0.7, panel.first=grid()) 
dev.off()
jpeg("Distributions_raw_counts_and_normalized_counts_shoot.jpeg",width =8.447917,height=6.072917,units="in",res=1000)
par(mfrow=c(2,2),mar=c(2,7,1,0.5),cex.lab=0.7)#set parameters for the plotting window
epsilon <- 1 # pseudo-count to avoid problems with log(0)
boxplot(log2(counts(dds_shoot.norm)+epsilon),col=c("red","red","red","grey","grey","grey"), cex.axis=0.7, 
        las=1, xlab="log2(counts+1)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds_shoot.norm, normalized=TRUE)+epsilon), col=c("red","red","red","grey","grey","grey"),cex.axis=0.7, 
        las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts") 
plotDensity(log2(counts(dds_shoot.norm)+epsilon), 
            xlab="log2(counts+1)", col=c("red","red","red","grey","grey","grey"),cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds_shoot.norm, normalized=TRUE)+epsilon), 
            xlab="log2(normalized counts)",col=c("red","red","red","grey","grey","grey"), cex.lab=0.7, panel.first=grid()) 
dev.off()

# Restore default parameters
par(mfrow=c(1,1), cex.lab=1,mar=c(5.1, 4.1, 4.1, 2.1))

# Performing estimation of dispersion parameter to account for different variances within populations for each gene
dds_root.disp <- estimateDispersions(dds_root.norm, fitType='local')
jpeg("Dispest_root_local.jpeg",width =8,height=8,units="in",res=1000)
plotDispEsts(dds_root.disp)
dev.off()

dds_root.disp <- estimateDispersions(dds_root.norm, fitType='parametric')
jpeg("Dispest_root_parametric.jpeg",width =8,height=8,units="in",res=1000)
plotDispEsts(dds_root.disp)
dev.off()

dds_root.disp <- estimateDispersions(dds_root.norm, fitType='mean')
jpeg("Dispest_root_mean.jpeg",width =8,height=8,units="in",res=1000)
plotDispEsts(dds_root.disp)
dev.off()

#go with default parametric, local automatically substituted if parametric does not converge,local gives more significant genes
#### Calculate Differential Expression
alpha <- 0.05 #significance level for adjusted p-value, default is 0.1
waldTestResult <- nbinomWaldTest(dds_root.disp)
resultDESeq2 <- results(waldTestResult, alpha=alpha, pAdjustMethod="BH")

head(resultDESeq2)

# Instead of normalizing and calculating differential expression stepwise, we can also run DESeq() on our "dds" data set. 
# Run the DESeq calculation (it might take a couple of minutes)
#dds_root$Population <- relevel(dds_root$Population, ref = "Krom")
#dds_root$Cu_condition <- relevel(dds_root$Cu_condition, ref = "normal_Cu")

design(dds_root)<- ~ Population + Cu_condition + Population:Cu_condition
design(dds_shoot)<- ~ Population + Cu_condition + Population:Cu_condition

dds_root_calc = DESeq(dds_root,test="Wald",fitType="local")
dds_shoot_calc = DESeq(dds_shoot,test="Wald",fitType="local")
resultsNames(dds_root_calc)
All<-results(dds_root_calc)
Kosi_Cu_treatment_R<-results(dds_root_calc,contrast=c("Cu_condition","high_Cu","normal_Cu"))
Krom_Cu_treatment_R<-results(dds_root_calc,list(c("Cu_condition_high_Cu_vs_normal_Cu","PopulationKrom.Cu_conditionhigh_Cu")))
IA_R<-results(dds_root_calc,alpha=0.05)
Krom_Kosi_R<-results(dds_root_calc,contrast=list("Population_Krom_vs_Kosi"))
Krom_Kosi_S<-results(dds_shoot_calc,contrast=list("Population_Krom_vs_Kosi"))

#qvalues_IA_R<-qvalue(na.omit(IA_R$pvalue))
#same result: no significance 

Kosi_Cu_treatment_R_sig<-as.data.frame(Kosi_Cu_treatment_R[which(Kosi_Cu_treatment_R$padj<0.05),])
Krom_Cu_treatment_R_sig<-as.data.frame(Krom_Cu_treatment_R[which(Krom_Cu_treatment_R$padj<0.05),])
IA_R_sig<-as.data.frame(IA_R[which(IA_R$padj<0.05),])
Krom_Kosi_R_sig<-as.data.frame(Krom_Kosi_R[which((Krom_Kosi_R$padj<0.05)&(Krom_Kosi_R$log2FoldChange>=1|Krom_Kosi_R$log2FoldChange<=(-1))),])
Krom_Kosi_S_sig<-as.data.frame(Krom_Kosi_S[which((Krom_Kosi_S$padj<0.05)&(Krom_Kosi_S$log2FoldChange>=1|Krom_Kosi_S$log2FoldChange<=(-1))),])
#IA not significant

resultsNames(dds_shoot_calc)
Kosi_Cu_treatment_S<-results(dds_shoot_calc,contrast=c("Cu_condition","high_Cu","normal_Cu"))
Krom_Cu_treatment_S<-results(dds_shoot_calc,list(c("Cu_condition_high_Cu_vs_normal_Cu","PopulationKrom.Cu_conditionhigh_Cu")))
IA_S<-results(dds_shoot_calc,name="PopulationKrom.Cu_conditionhigh_Cu",alpha=0.05)
Krom_Kosi_S<-results(dds_shoot_calc,contrast=list("Population_Krom_vs_Kosi"))

Kosi_Cu_treatment_S_sig<-as.data.frame(Kosi_Cu_treatment_S[which(Kosi_Cu_treatment_S$padj<0.05),])
Krom_Cu_treatment_S_sig<-as.data.frame(Krom_Cu_treatment_S[which(Krom_Cu_treatment_S$padj<0.05),])
IA_S_sig<-as.data.frame(IA_S[which(IA_S$padj<0.1),])
Krom_Kosi_S_sig<-as.data.frame(Krom_Kosi_S[which(Krom_Kosi_S$padj<0.05),])
#IA not significant
IA_S_sorted<-IA_S_sig[order(IA_S_sig$padj),]
IA_R_sorted<-IA_R_sig[order(IA_R_sig$padj),]


#different design with same meaning
dds_root$IA<-factor(paste0(dds_root$Population,dds_root$Cu_condition))
dds_shoot$IA<-factor(paste0(dds_shoot$Population,dds_shoot$Cu_condition))
design(dds_root)<-~IA
design(dds_shoot)<-~IA
dds_root_calc = DESeq(dds_root,test="Wald",fitType="local")
dds_shoot_calc = DESeq(dds_shoot,test="Wald",fitType="local")
resultsNames(dds_root_calc)
resultsNames(dds_shoot_calc)
IA_R_Krom<-results(dds_root_calc,contrast=c("IA","Kromhigh_Cu","Kromnormal_Cu"),alpha=0.05)
IA_R_Kosi<-results(dds_root_calc,contrast=c("IA","Kosihigh_Cu","Kosinormal_Cu"),alpha=0.05)
IA_R_Krom_Kosi_highCu<-results(dds_root_calc,contrast=c("IA","Kromhigh_Cu","Kosihigh_Cu"),alpha=0.05)
IA_R_Krom_Kosi_normalCu<-results(dds_root_calc,contrast=c("IA","Kromnormal_Cu","Kosinormal_Cu"),alpha=0.05)

IA_R_Krom_sig<-as.data.frame(IA_R_Krom[which(IA_R_Krom$padj<=0.05&(IA_R_Krom$log2FoldChange>=1|IA_R_Krom$log2FoldChange<=(-1))),])
IA_R_Kosi_sig<-as.data.frame(IA_R_Kosi[which(IA_R_Kosi$padj<=0.05&(IA_R_Kosi$log2FoldChange>=1|IA_R_Kosi$log2FoldChange<=(-1))),])
IA_R_Krom_Kosi_highCu_sig<-as.data.frame(IA_R_Krom_Kosi_highCu[which(IA_R_Krom_Kosi_highCu$padj<=0.05&(IA_R_Krom_Kosi_highCu$log2FoldChange>=1|IA_R_Krom_Kosi_highCu$log2FoldChange<=(-1))),])
IA_R_Krom_Kosi_normalCu_sig<-as.data.frame(IA_R_Krom_Kosi_normalCu[which(IA_R_Krom_Kosi_normalCu$padj<=0.05&(IA_R_Krom_Kosi_normalCu$log2FoldChange>=1|IA_R_Krom_Kosi_normalCu$log2FoldChange<=(-1))),])


IA_S_Krom<-results(dds_shoot_calc,contrast=c("IA","Kromhigh_Cu","Kromnormal_Cu"),alpha=0.05)
IA_S_Kosi<-results(dds_shoot_calc,contrast=c("IA","Kosihigh_Cu","Kosinormal_Cu"),alpha=0.05)
IA_S_Krom_Kosi_highCu<-results(dds_shoot_calc,contrast=c("IA","Kromhigh_Cu","Kosihigh_Cu"),alpha=0.05)
IA_S_Krom_Kosi_normalCu<-results(dds_shoot_calc,contrast=c("IA","Kromnormal_Cu","Kosinormal_Cu"),alpha=0.05)

IA_S_Krom_sig<-as.data.frame(IA_S_Krom[which(IA_S_Krom$padj<=0.05&(IA_S_Krom$log2FoldChange>=1|IA_S_Krom$log2FoldChange<=(-1))),])
IA_S_Kosi_sig<-as.data.frame(IA_S_Kosi[which(IA_S_Kosi$padj<=0.05&(IA_S_Kosi$log2FoldChange>=1|IA_S_Kosi$log2FoldChange<=(-1))),])
IA_S_Krom_Kosi_highCu_sig<-as.data.frame(IA_S_Krom_Kosi_highCu[which(IA_S_Krom_Kosi_highCu$padj<=0.05&(IA_S_Krom_Kosi_highCu$log2FoldChange>=1|IA_S_Krom_Kosi_highCu$log2FoldChange<=(-1))),])
IA_S_Krom_Kosi_normalCu_sig<-as.data.frame(IA_S_Krom_Kosi_normalCu[which(IA_S_Krom_Kosi_normalCu$padj<=0.05&(IA_S_Krom_Kosi_normalCu$log2FoldChange>=1|IA_S_Krom_Kosi_normalCu$log2FoldChange<=(-1))),])

IA_S_Krom_sig_nolog2foldchange<-as.data.frame(IA_S_Krom[which(IA_S_Krom$padj<=0.05),])
IA_S_Kosi_sig_nolog2foldchange<-as.data.frame(IA_S_Kosi[which(IA_S_Kosi$padj<=0.05),])
IA_R_Krom_sig_nolog2foldchange<-as.data.frame(IA_R_Krom[which(IA_R_Krom$padj<=0.05),])
IA_R_Kosi_sig_nolog2foldchange<-as.data.frame(IA_R_Kosi[which(IA_R_Kosi$padj<=0.05),])



#NormalizedCount_all <- counts(dds.SF, normalized = T) #Extract the normalized counts
#head(NormalizedCount_all)
#write.table(NormalizedCount_all, file="Normalized_counts_woGTF.txt", sep="\t", quote=F, col.names=NA)

pdf("Cooks_distance_root_for_outliers.pdf",width=4.72441,height=4.72441,paper="special")
boxplot(log10(assays(dds_root_calc)[["cooks"]]), range=0, las=2, main="Cook's distances") # Check any outliers by Cook's D 
dev.off()
pdf("Cooks_distance_shoot_for_outliers.pdf",width=4.72441,height=4.72441,paper="special")
boxplot(log10(assays(dds_shoot_calc)[["cooks"]]), range=0, las=2, main="Cook's distances") # Check any outliers by Cook's D 
dev.off()

#merge by rownames with orthologues and thaliana gene description
OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

FPKM2<-FPKM_df[,c(27:50)]
colnames(FPKM2)<-coldata$Sample.Name
rownames(FPKM2)<-FPKM_df[,1]
FPKM2<-FPKM2[,c(19,11,3,20,12,8,24,16,4,21,13,5,17,9,1,18,10,2,23,15,7,22,14,6)]
Krom_R_Cu_Ah<-merge(IA_R_Krom_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Krom_R_Cu_Ah)[1]<-"Ah_gene"
Krom_R_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_R_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Kosi_R_Cu_Ah<-merge(IA_R_Kosi_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Kosi_R_Cu_Ah)[1]<-"Ah_gene"
Kosi_R_Cu_Ath_desc<-merge(OG_Ahalleri,Kosi_R_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_R_high_Cu_Ah<-merge(IA_R_Krom_Kosi_highCu_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Krom_Kosi_R_high_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_R_high_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_R_high_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_R_normal_Cu_Ah<-merge(IA_R_Krom_Kosi_normalCu_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Krom_Kosi_R_normal_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_R_normal_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_R_normal_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_S_Cu_Ah<-merge(IA_S_Krom_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Krom_S_Cu_Ah)[1]<-"Ah_gene"
Krom_S_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_S_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Kosi_S_Cu_Ah<-merge(IA_S_Kosi_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Kosi_S_Cu_Ah)[1]<-"Ah_gene"
Kosi_S_Cu_Ath_desc<-merge(OG_Ahalleri,Kosi_S_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_S_high_Cu_Ah<-merge(IA_S_Krom_Kosi_highCu_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Krom_Kosi_S_high_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_S_high_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_S_high_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_S_normal_Cu_Ah<-merge(IA_S_Krom_Kosi_normalCu_sig,FPKM2,by="row.names",all.x=TRUE)
colnames(Krom_Kosi_S_normal_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_S_normal_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_S_normal_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

#sort by fold change
Krom_R_Cu_Ath_desc_sorted<-Krom_R_Cu_Ath_desc[order(-abs(Krom_R_Cu_Ath_desc$log2FoldChange)),]
Kosi_R_Cu_Ath_desc_sorted<-Kosi_R_Cu_Ath_desc[order(-abs(Kosi_R_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_R_high_Cu_Ath_desc_sorted<-Krom_Kosi_R_high_Cu_Ath_desc[order(-abs(Krom_Kosi_R_high_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_R_normal_Cu_Ath_desc_sorted<-Krom_Kosi_R_normal_Cu_Ath_desc[order(-abs(Krom_Kosi_R_normal_Cu_Ath_desc$log2FoldChange)),]
Krom_S_Cu_Ath_desc_sorted<-Krom_S_Cu_Ath_desc[order(-abs(Krom_S_Cu_Ath_desc$log2FoldChange)),]
Kosi_S_Cu_Ath_desc_sorted<-Kosi_S_Cu_Ath_desc[order(-abs(Kosi_S_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_S_high_Cu_Ath_desc_sorted<-Krom_Kosi_S_high_Cu_Ath_desc[order(-abs(Krom_Kosi_S_high_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_S_normal_Cu_Ath_desc_sorted<-Krom_Kosi_S_normal_Cu_Ath_desc[order(-abs(Krom_Kosi_S_normal_Cu_Ath_desc$log2FoldChange)),]

#add gff description and info
gff<-read.table("AHAL_sorted.gff3",sep="\t",quote="")
gff_gene<-gff[gff$V3=="gene",]
Anno_genea<-strsplit(as.character(gff_gene$V9),"ID=",fixed=T)
Anno_geneb<-lapply(Anno_genea,"[",2)
Anno_genec<-data.frame(matrix(unlist(Anno_geneb),nrow=length(Anno_geneb)))
Anno_genee<-strsplit(as.character(Anno_genec[,1]),";",fixed=T)
Anno_genef<-lapply(Anno_genee,"[",1)
Anno_gene<-data.frame(matrix(unlist(Anno_genef),nrow=length(Anno_genef)))
gff_gene$Gene<-Anno_gene[,1]

Krom_R_Cu_Ath_desc_sorted_gff<-merge(Krom_R_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Kosi_R_Cu_Ath_desc_sorted_gff<-merge(Kosi_R_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_R_high_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_R_normal_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_S_Cu_Ath_desc_sorted_gff<-merge(Krom_S_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Kosi_S_Cu_Ath_desc_sorted_gff<-merge(Kosi_S_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_S_high_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_S_normal_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)

#write.table(Krom_R_Cu_Ath_desc_sorted_gff,"Krom_root_Cu_reaction.table",sep="\t",row.names=F)
#write.table(Kosi_R_Cu_Ath_desc_sorted_gff,"Kosi_root_Cu_reaction.table",sep="\t",row.names=F)
#write.table(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff,"Krom_Kosi_root_high_Cu.table",sep="\t",row.names=F)
#write.table(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff,"Krom_Kosi_root_normal_Cu.table",sep="\t",row.names=F)
#write.table(Krom_S_Cu_Ath_desc_sorted_gff,"Krom_shoot_Cu_reaction.table",sep="\t",row.names=F)
#write.table(Kosi_S_Cu_Ath_desc_sorted_gff,"Kosi_shoot_Cu_reaction.table",sep="\t",row.names=F)
#write.table(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff,"Krom_Kosi_shoot_high_Cu.table",sep="\t",row.names=F)
#write.table(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff,"Krom_Kosi_shoot_normal_Cu.table",sep="\t",row.names=F)

#Unknowns<-c(Krom_R_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Krom_R_Cu_Ath_desc_sorted_gff$Ath)],Kosi_R_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Kosi_R_Cu_Ath_desc_sorted_gff$Ath)],Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff$Ath)],Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff$Ath)],
#Krom_S_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Krom_S_Cu_Ath_desc_sorted_gff$Ath)],Kosi_S_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Kosi_S_Cu_Ath_desc_sorted_gff$Ath)],Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff$Ath)],Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff$Ah_gene[is.na(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff$Ath)])
#Unknowns2<-unique(Unknowns)
#write.table(Unknowns2,"Unknowns_genes_Cu.table",row.names=F,col.names=F,quote=F)

colnames(Krom_R_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Kosi_R_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Kosi_S_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")

#change order
Krom_R_Cu_Ath_desc_sorted_gff<-Krom_R_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Kosi_R_Cu_Ath_desc_sorted_gff<-Kosi_R_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff<-Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff<-Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Kosi_S_Cu_Ath_desc_sorted_gff<-Kosi_S_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff<-Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff<-Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]

#add lists
HM_list<-read.xlsx2("metal_homeostasis_Ute_2018_04_23_later.xlsx",1,header=T)
HM_list$AGI_Number<-toupper(HM_list$AGI_Number)
Krom_R_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_R_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
Kosi_R_Cu_Ath_desc_sorted_gff_lists<-merge(Kosi_R_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
Kosi_S_Cu_Ath_desc_sorted_gff_lists<-merge(Kosi_S_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff,HM_list,all.x=T,by.x="Ath_ID",by.y="AGI_Number")

write.table(Krom_R_Cu_Ath_desc_sorted_gff_lists,"Krom_root_Cu_reaction_withBlast_stricter.table",sep="\t",row.names=F)
write.table(Kosi_R_Cu_Ath_desc_sorted_gff_lists,"Kosi_root_Cu_reaction_withBlast_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_root_high_Cu_withBlast_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_root_normal_Cu_withBlast_stricter.table",sep="\t",row.names=F)
#write.table(Krom_S_Cu_Ath_desc_sorted_gff_lists,"Krom_shoot_Cu_reaction_withBlast_stricter.table",sep="\t",row.names=F)
write.table(Kosi_S_Cu_Ath_desc_sorted_gff_lists,"Kosi_shoot_Cu_reaction_withBlast_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_shoot_high_Cu_withBlast_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_shoot_normal_Cu_withBlast_stricter.table",sep="\t",row.names=F)



###################################################
#Datasets with median of ratios/ transcript length#
###################################################
#separate root and shoot
dds_root <- estimateSizeFactors(dds_root)
sizeFactors(dds_root)
normalized_counts_R <-data.frame(counts(dds_root, normalized=TRUE))
colnames(normalized_counts_R)<-dds_root$Sample.Name

dds_shoot <- estimateSizeFactors(dds_shoot)
sizeFactors(dds_shoot)
normalized_counts_S <-data.frame(counts(dds_shoot, normalized=TRUE))
colnames(normalized_counts_S)<-dds_shoot$Sample.Name

countdata_df<-cbind(normalized_counts_R,normalized_counts_S)
countdata_rlog<-countdata_df[,c(9,5,1,10,6,4,12,8,2,11,7,3,21,17,13,22,18,14,24,20,16,23,19,15)]

countdata_rlog2<-merge(CDS_df[,c(1,3)],countdata_rlog,by.x="seq_name",by.y="row.names",all.y=T)
countdata_rlog2$Tlen[is.na(countdata_rlog2$Tlen)]<-PS_Hyp_genes$Tlen[match(countdata_rlog2$seq_name[is.na(countdata_rlog2$Tlen)],PS_Hyp_genes$V10)]

countdata_rlog2[,3:26]<-countdata_rlog2[,3:26]*1000/countdata_rlog2$Tlen
countdata_rlog_df<-countdata_rlog2[,-2]
colnames(countdata_rlog_df)[2:25]<-c("Krom_Cu_R1","Krom_Cu_R2","Krom_Cu_R3","Krom_Ctrl_R1","Krom_Ctrl_R2","Krom_Ctrl_R3","Kosi_Cu_R1","Kosi_Cu_R2","Kosi_Cu_R3","Kosi_Ctrl_R1","Kosi_Ctrl_R2","Kosi_Ctrl_R3","Krom_Cu_S1","Krom_Cu_S2","Krom_Cu_S3","Krom_Ctrl_S1","Krom_Ctrl_S2","Krom_Ctrl_S3","Kosi_Cu_S1","Kosi_Cu_S2","Kosi_Cu_S3","Kosi_Ctrl_S1","Kosi_Ctrl_S2","Kosi_Ctrl_S3")

countdata_rlog_df_Tair<-merge(countdata_rlog_df,OG_Ahalleri,by.x="seq_name",by.y="Ah_ID",all.x=TRUE)

write.table(countdata_rlog_df_Tair,"Countdata_medofratios.table",sep="\t",row.names=F)






Krom_R_Cu_Ah<-merge(IA_R_Krom_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Krom_R_Cu_Ah)[1]<-"Ah_gene"
Krom_R_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_R_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Kosi_R_Cu_Ah<-merge(IA_R_Kosi_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Kosi_R_Cu_Ah)[1]<-"Ah_gene"
Kosi_R_Cu_Ath_desc<-merge(OG_Ahalleri,Kosi_R_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_R_high_Cu_Ah<-merge(IA_R_Krom_Kosi_highCu_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Krom_Kosi_R_high_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_R_high_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_R_high_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_R_normal_Cu_Ah<-merge(IA_R_Krom_Kosi_normalCu_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Krom_Kosi_R_normal_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_R_normal_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_R_normal_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_S_Cu_Ah<-merge(IA_S_Krom_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Krom_S_Cu_Ah)[1]<-"Ah_gene"
Krom_S_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_S_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Kosi_S_Cu_Ah<-merge(IA_S_Kosi_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Kosi_S_Cu_Ah)[1]<-"Ah_gene"
Kosi_S_Cu_Ath_desc<-merge(OG_Ahalleri,Kosi_S_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_S_high_Cu_Ah<-merge(IA_S_Krom_Kosi_highCu_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Krom_Kosi_S_high_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_S_high_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_S_high_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

Krom_Kosi_S_normal_Cu_Ah<-merge(IA_S_Krom_Kosi_normalCu_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
colnames(Krom_Kosi_S_normal_Cu_Ah)[1]<-"Ah_gene"
Krom_Kosi_S_normal_Cu_Ath_desc<-merge(OG_Ahalleri,Krom_Kosi_S_normal_Cu_Ah,by.y="Ah_gene",by.x="Ah_ID",all.y=TRUE)

#sort by fold change
Krom_R_Cu_Ath_desc_sorted<-Krom_R_Cu_Ath_desc[order(-abs(Krom_R_Cu_Ath_desc$log2FoldChange)),]
Kosi_R_Cu_Ath_desc_sorted<-Kosi_R_Cu_Ath_desc[order(-abs(Kosi_R_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_R_high_Cu_Ath_desc_sorted<-Krom_Kosi_R_high_Cu_Ath_desc[order(-abs(Krom_Kosi_R_high_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_R_normal_Cu_Ath_desc_sorted<-Krom_Kosi_R_normal_Cu_Ath_desc[order(-abs(Krom_Kosi_R_normal_Cu_Ath_desc$log2FoldChange)),]
Krom_S_Cu_Ath_desc_sorted<-Krom_S_Cu_Ath_desc[order(-abs(Krom_S_Cu_Ath_desc$log2FoldChange)),]
Kosi_S_Cu_Ath_desc_sorted<-Kosi_S_Cu_Ath_desc[order(-abs(Kosi_S_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_S_high_Cu_Ath_desc_sorted<-Krom_Kosi_S_high_Cu_Ath_desc[order(-abs(Krom_Kosi_S_high_Cu_Ath_desc$log2FoldChange)),]
Krom_Kosi_S_normal_Cu_Ath_desc_sorted<-Krom_Kosi_S_normal_Cu_Ath_desc[order(-abs(Krom_Kosi_S_normal_Cu_Ath_desc$log2FoldChange)),]

#add gff description and info
gff<-read.table("AHAL_sorted2.gff3",sep="\t",quote="")
gff_gene<-gff[gff$V3=="gene",]
Anno_genea<-strsplit(as.character(gff_gene$V9),"ID=",fixed=T)
Anno_geneb<-lapply(Anno_genea,"[",2)
Anno_genec<-data.frame(matrix(unlist(Anno_geneb),nrow=length(Anno_geneb)))
Anno_genee<-strsplit(as.character(Anno_genec[,1]),";",fixed=T)
Anno_genef<-lapply(Anno_genee,"[",1)
Anno_gene<-data.frame(matrix(unlist(Anno_genef),nrow=length(Anno_genef)))
gff_gene$Gene<-Anno_gene[,1]

Krom_R_Cu_Ath_desc_sorted_gff<-merge(Krom_R_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Kosi_R_Cu_Ath_desc_sorted_gff<-merge(Kosi_R_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_R_high_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_R_normal_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_S_Cu_Ath_desc_sorted_gff<-merge(Krom_S_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Kosi_S_Cu_Ath_desc_sorted_gff<-merge(Kosi_S_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_S_high_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff<-merge(Krom_Kosi_S_normal_Cu_Ath_desc_sorted,gff_gene,by.x="Ah_ID",by.y="Gene",all.x=T)

colnames(Krom_R_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Kosi_R_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Kosi_S_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")
colnames(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff)[40:48]<-c("Scaffold","Source","Type","Start","End","Score","Strand","Phase","Attributes")

#change order
Krom_R_Cu_Ath_desc_sorted_gff<-Krom_R_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Kosi_R_Cu_Ath_desc_sorted_gff<-Kosi_R_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff<-Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff<-Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Kosi_S_Cu_Ath_desc_sorted_gff<-Kosi_S_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff<-Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff<-Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff[,c(1,40:48,2:39)]



Mapman<-read.table("F:/Cu_project/20190828_MAPMAN_Ath_AGI_LOCUS_TAIR10_UK_LabEdit_strict_LS.txt",header=T,sep="\t",fill=T)
Mapman$Category<-paste(Mapman$BINCODE, Mapman$NAME,sep=":")
Mapman$Main<-as.numeric(sapply(strsplit(as.character(Mapman$BINCODE),split=".",fixed=T),'[',1))
Mapman$Main_Name<-sapply(strsplit(as.character(Mapman$NAME),split=".",fixed=T),'[',1)
Mapman$Main_category<-paste(Mapman$Main, Mapman$Main_Name,sep=":")
Mapman<-Mapman[!(is.na(Mapman$IDENTIFIER)|Mapman$IDENTIFIER==""),]
write.table(Mapman,"Mapman_for_merge.table",sep="\t",row.names=F,quote=F)

Mapman<-read.table("Mapman_for_merge.table",header=T,sep="\t",fill=T)
Mapman<-unique(Mapman)
Mapman$IDENTIFIER<-toupper(Mapman$IDENTIFIER)
#460 no AGI code but name, with M for metabolite against T for transcript
#str(Mapman[Mapman$TYPE=="M",])
Mapman<-Mapman[Mapman$TYPE=="T",]
Mapman_df<-unique(Mapman[,c(3,9)])
require(data.table)
DT<-data.table(Mapman_df)
setkey(DT,IDENTIFIER)
Mapman_dt<-DT[,paste(Main_category,collapse="/"),by=IDENTIFIER]

Krom_R_Cu_Ath_desc_sorted_gff_M<-merge(Krom_R_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Krom_S_Cu_Ath_desc_sorted_gff_M<-merge(Krom_S_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Kosi_R_Cu_Ath_desc_sorted_gff_M<-merge(Kosi_R_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Kosi_S_Cu_Ath_desc_sorted_gff_M<-merge(Kosi_S_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_M<-merge(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_M<-merge(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_M<-merge(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_M<-merge(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff,Mapman_dt,by.x="Ath_ID",by.y="IDENTIFIER",all.x=T)

colnames(Krom_R_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Krom_S_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Kosi_R_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Kosi_S_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"
colnames(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_M)[49]<-"Mapman_category"

write.table(Krom_R_Cu_Ath_desc_sorted_gff_M,"Krom_root_Cu_reaction_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Kosi_R_Cu_Ath_desc_sorted_gff_M,"Kosi_root_Cu_reaction_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_M,"Krom_Kosi_root_high_Cu_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_M,"Krom_Kosi_root_normal_Cu_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Krom_S_Cu_Ath_desc_sorted_gff_M,"Krom_shoot_Cu_reaction_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Kosi_S_Cu_Ath_desc_sorted_gff_M,"Kosi_shoot_Cu_reaction_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_M,"Krom_Kosi_shoot_high_Cu_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)
write.table(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_M,"Krom_Kosi_shoot_normal_Cu_RBH_medofratios_stricter_M.table",sep="\t",row.names=F)

#add lists
HM_list<-read.csv("F:/Cu_project/metal_homeostasis_2019_v4_stringent.csv",header=T)
HM_list$AGI_Number<-toupper(HM_list$AGI_Number)
HM_list<-HM_list[,c(2:4,6,8,9,11,14,15,16)]
colnames(HM_list)[1]<-"AGI_Code"
Cu_proteins<-read.xlsx2("F:/Cu_project/Cu_proteins_integrated_mod_06_2013.xlsx",3,header=T)
Cu_proteins$AGI_Id<-toupper(Cu_proteins$AGI_Id)
colnames(Cu_proteins)[1]<-"AGI_Code"
DNA_repair<-read.xlsx2("F:/Cu_project/At_DNA_dam_rep_genes_short.xlsx",2,header=T)
DNA_repair<-DNA_repair[,c(2:3,5,6,7,8)]
colnames(DNA_repair)<-c("AGI_Code","Name","Pathway","Function","Activity","Comment")
Zn_metalloproteins<-read.xlsx2("F:/Cu_project/Zn_Cu_Mn_Fe_heme_Fe_ion_metalloproteins_Table_5_Zhang_FPLS_LS.XLS",1,header=T)
colnames(Zn_metalloproteins)[1]<-"AGI_Code"
Mn_metalloproteins<-read.xlsx2("F:/Cu_project/Zn_Cu_Mn_Fe_heme_Fe_ion_metalloproteins_Table_5_Zhang_FPLS_LS.XLS",3,header=T)
colnames(Mn_metalloproteins)[1]<-"AGI_Code"
Heme_metalloproteins<-read.xlsx2("F:/Cu_project/Zn_Cu_Mn_Fe_heme_Fe_ion_metalloproteins_Table_5_Zhang_FPLS_LS.XLS",4,header=T)
colnames(Heme_metalloproteins)[1]<-"AGI_Code"
Fe_metalloproteins<-read.xlsx2("F:/Cu_project/Zn_Cu_Mn_Fe_heme_Fe_ion_metalloproteins_Table_5_Zhang_FPLS_LS.XLS",5,header=T)
colnames(Fe_metalloproteins)[1]<-"AGI_Code"
FeS_biogenesis<-read.xlsx2("F:/Cu_project/FeS_cluster_assembly_Table_2_Zhang_FPLS_LS.XLS",1,header=T)
colnames(FeS_biogenesis)[2]<-"AGI_Code"
FeS_biogenesis<-FeS_biogenesis[!(is.na(FeS_biogenesis$AGI_Code)|FeS_biogenesis$AGI_Code==""),]
FeS_use<-read.xlsx2("F:/Cu_project/FeS_cluster_use_Table_1_Zhang_FPLS_LS.xlsx",1,header=T)
colnames(FeS_use)[2]<-"AGI_Code"
FeS_use<-FeS_use[!(is.na(FeS_use$AGI_Code)|FeS_use$AGI_Code==""),]
Meiosis<-read.xlsx2("F:/Cu_project/meiosis_all_update_2019_09_17.xlsx",1,header=T)
colnames(Meiosis)<-c("AGI_Code","Gene.short.name","Name","Meiotic.phase","Description_custom")
Casp<-read.xlsx2("F:/Cu_project/Casparian_strip_Suberin_2019_Baohai_LS.xlsx",1,header=T)
Casp<-Casp[,c(1:4,6:7)]
colnames(Casp)[2]<-"AGI_Code"
Suberin<-read.xlsx2("F:/Cu_project/Casparian_strip_Suberin_2019_Baohai_LS.xlsx",2,header=T)
Suberin<-Suberin[,c(1:3,5)]
colnames(Suberin)[1]<-"AGI_Code"

All_lists<-Reduce(function(x,y) merge(x,y,all=T,by="AGI_Code"),list(HM_list,Cu_proteins,DNA_repair,Zn_metalloproteins,Mn_metalloproteins,Heme_metalloproteins,Fe_metalloproteins,FeS_biogenesis,FeS_use,Meiosis,Casp,Suberin))
write.table(All_lists,"All_lists_Cu_project.table",sep="\t",row.names=F,quote=F)

Krom_R_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_R_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Kosi_R_Cu_Ath_desc_sorted_gff_lists<-merge(Kosi_R_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Kosi_S_Cu_Ath_desc_sorted_gff_lists<-merge(Kosi_S_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists<-merge(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_M,All_lists,all.x=T,by.x="Ath_ID",by.y="AGI_Code")

Krom_R_Cu_Ath_desc_sorted_gff_lists<-Krom_R_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]
Kosi_R_Cu_Ath_desc_sorted_gff_lists<-Kosi_R_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists<-Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists<-Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]
Kosi_S_Cu_Ath_desc_sorted_gff_lists<-Kosi_S_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists<-Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists<-Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists[,c(2:11,1,12:92)]


write.xlsx2(unique(Krom_R_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Krom_R_Cu_reaction",row.names=F)
write.xlsx2(unique(Kosi_R_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Kosi_R_Cu_reaction",append=T,row.names=F)
write.xlsx2(unique(Kosi_S_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Kosi_S_Cu_reaction",append=T,row.names=F)
write.xlsx2(unique(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Krom_Kosi_R_high_Cu",append=T,row.names=F)
write.xlsx2(unique(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Krom_Kosi_S_high_Cu",append=T,row.names=F)
write.xlsx2(unique(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Krom_Kosi_R_normal_Cu",append=T,row.names=F)
write.xlsx2(unique(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists),"Deseq2_medofratios_stricter_Mapmanadded.xlsx",sheetName="Krom_Kosi_S_normal_Cu",append=T,row.names=F)

write.table(Krom_R_Cu_Ath_desc_sorted_gff_lists,"Krom_root_Cu_reaction_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Kosi_R_Cu_Ath_desc_sorted_gff_lists,"Kosi_root_Cu_reaction_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_root_high_Cu_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_root_normal_Cu_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Krom_S_Cu_Ath_desc_sorted_gff_lists,"Krom_shoot_Cu_reaction_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Kosi_S_Cu_Ath_desc_sorted_gff_lists,"Kosi_shoot_Cu_reaction_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_shoot_high_Cu_RBH_medofratios_stricter.table",sep="\t",row.names=F)
write.table(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff_lists,"Krom_Kosi_shoot_normal_Cu_RBH_medofratios_stricter.table",sep="\t",row.names=F)


pdf( "Candidates_KromKosiconst_final_divided.pdf",width=12,height=12,pointsize=13)

layout(matrix(c(1:20),4,5,byrow=TRUE))
Candidates_Kromconst_Kosireac<-countdata_rlog_df_Tair[countdata_rlog_df_Tair$seq_name=="AHAL_G0006247"|countdata_rlog_df_Tair$seq_name=="AHAL_G0027080"|countdata_rlog_df_Tair$seq_name=="AHAL_G0003607"|countdata_rlog_df_Tair$seq_name=="AHAL_G0006030"|countdata_rlog_df_Tair$seq_name=="AHAL_G0019842"|countdata_rlog_df_Tair$seq_name=="AHAL_G0006748"|countdata_rlog_df_Tair$seq_name=="AHAL_G0016520"|countdata_rlog_df_Tair$seq_name=="AHAL_G0006377"|countdata_rlog_df_Tair$seq_name=="AHAL_G0029705"
|countdata_rlog_df_Tair$seq_name=="AHAL_G0018941"|countdata_rlog_df_Tair$seq_name=="AHAL_G0010517"|countdata_rlog_df_Tair$seq_name=="AHAL_G0017109"|countdata_rlog_df_Tair$seq_name=="AHAL_G0024785"|countdata_rlog_df_Tair$seq_name=="AHAL_G0006054"|countdata_rlog_df_Tair$seq_name=="AHAL_G0010288"|countdata_rlog_df_Tair$seq_name=="AHAL_G0002779"|countdata_rlog_df_Tair$seq_name=="AHAL_G0011944"|countdata_rlog_df_Tair$seq_name=="AHAL_G0014700"|countdata_rlog_df_Tair$seq_name=="AHAL_G0007128"|
countdata_rlog_df_Tair$seq_name=="AHAL_G0017516"|countdata_rlog_df_Tair$seq_name=="AHAL_G0007270"|countdata_rlog_df_Tair$seq_name=="AHAL_G000637",]
Candidates_Kromconst_Kosireac2<-Candidates_Kromconst_Kosireac[!duplicated(Candidates_Kromconst_Kosireac$seq_name),]
Candidates_R<-Candidates_Kromconst_Kosireac2[c(20,9,7,5,2,18,13,8,3,16),]

par(mar=c(6,3,3,1))
par(oma=c(0,4,0,3))
Sub<-c("AT5G35525","AT1G71140","AT5G44430","AT1G15960","AT1G21850","AT3G56980","AT2G38460","AT1G72520","AT1G01480","AT3G13226")


Main<-c("PCR3","DTX14","PDF1.2c","NRAMP6","SKS8","bHLH39","IREG1","LOX4","ACS2","Regul. prot. RecX fam.")
for (i in 1:10)
{
colours<-c("black","red")
Cand<-cbind(as.numeric(t(Candidates_R[i,11:13])),as.numeric(t(Candidates_R[i,8:10])),as.numeric(t(Candidates_R[i,5:7])),as.numeric(t(Candidates_R[i,2:4])))
bpl_Cand<-boxplot(t(Cand)~c(1:4),names=NA,main="",las=2,boxcol=colours,whiskcol=colours,staplecol=colours,medcol=colours,cex.names=1.3)
Mean <- apply(Cand,2,mean)
bpl_Cand$stats[3,]<-Mean
bxp(bpl_Cand,medlty="dotted",names=NA,main="",las=2,boxcol=colours,whiskcol=colours,staplecol=colours,medcol=colours,add=T,cex.names=1.3)
boxplot(t(Cand)~c(1:4),names=NA,main="",las=2,add=T,boxcol=colours,whiskcol=colours,staplecol=colours,medcol=colours,cex.names=1.3)
text(c(expression(paste(bold(RS[KO])," ",bold(L),"Cu")),expression(paste(bold(RS[Ko])," ",bold(H),"Cu")),expression(paste(bold(CS[Kr])," ",bold(L),"Cu")),expression(paste(bold(CS[Kr])," ",bold(H),"Cu"))),cex=1.5,xpd=TRUE, srt=45,pos=1,x=c(0.5,1.5,2.5,3.5),y=par("usr")[3],adj=1,offset=2.2)
Lines <- list(Main[i],Sub[i])
mtext(do.call(expression, Lines),side=3,line=1:0,cex=1,padj=c(-0.5,0))
}
mtext(expression(Normalized~counts~(per~kb~transcript)),2,line=1,xpd=T,outer=T,cex=1.3)

Candidates_S<-Candidates_Kromconst_Kosireac2[c(19,2,7,10,15,21,17,12,6,1),]
Main_S<-c("HIPP05","SKS8","PDF1.2c","MOT1","Plantacyanin","Cupredoxin superfam.","LSU3","SULTR4;2","MHF2","similar to CLPC hom. 1")
Sub_S<-c("AT4G05030","AT1G21850","AT5G44430","AT2G25680","AT2G02850","AT4G39830","AT3G49570","AT3G13226","AT1G78790","similar to AT5G50920")
for (i in 1:10)
{
colours<-c("black","red")
Cand<-cbind(as.numeric(t(Candidates_S[i,23:25])),as.numeric(t(Candidates_S[i,20:22])),as.numeric(t(Candidates_S[i,17:19])),as.numeric(t(Candidates_S[i,14:16])))
bpl_Cand<-boxplot(t(Cand)~c(1:4),names=NA,main="",las=2,boxcol=colours,whiskcol=colours,staplecol=colours,medcol=colours,cex.names=1.3)
Mean <- apply(Cand,2,mean)
bpl_Cand$stats[3,]<-Mean
bxp(bpl_Cand,medlty="dotted",names=NA,main="",las=2,boxcol=colours,whiskcol=colours,staplecol=colours,medcol=colours,add=T,cex.names=1.3)
boxplot(t(Cand)~c(1:4),names=NA,main="",las=2,add=T,boxcol=colours,whiskcol=colours,staplecol=colours,medcol=colours,cex.names=1.3)
text(c(expression(paste(bold(RS[KO])," ",bold(L),"Cu")),expression(paste(bold(RS[Ko])," ",bold(H),"Cu")),expression(paste(bold(CS[Kr])," ",bold(L),"Cu")),expression(paste(bold(CS[Kr])," ",bold(H),"Cu"))),cex=1.5,xpd=TRUE, srt=45,pos=1,x=c(0.5,1.5,2.5,3.5),y=par("usr")[3],adj=1,offset=2.2)
Lines <- list(Main_S[i],Sub_S[i])
mtext(do.call(expression, Lines),side=3,line=1:0,cex=1,padj=c(-0.5,0))
}

dev.off()





Log <- rlog(dds, blind = FALSE)
LogS <- rlog(dds_shoot, blind = FALSE)
LogR <- rlog(dds_root, blind = FALSE)
pdf("PCA_all.pdf", width=8, height=8, paper="special")
par(oma=c(1,1,1,1))
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(1,1.35),heights=c(0.75,1))
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(Log, intgroup=c("Population", "Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=pcaData$Population,pch=c(16,17)[as.numeric(pcaData$Tissue)],xlab="",ylab="",main="Total",cex=1.5,cex.lab=1.2,axes=F)
legend("topright",inset=c(-0.16,0),legend=c("Krom","Kosi"),fill=c("red","black"),title="Population",bty="n",cex=1.2)
legend("topright",inset=c(-0.16,0.43),legend=c("Root","Shoot"),pch=c(16,17),title="Tissue",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,1))
pcaData <- plotPCA(LogR, intgroup=c("Population", "Cu_condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=pcaData$Population,pch=c(16,17)[as.numeric(pcaData$Cu_condition)],xlab="",ylab="",main="Root",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
par(mar=c(3,3,2,8),xpd=T)
pcaData <- plotPCA(LogS, intgroup=c("Population", "Cu_condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData$PC2~pcaData$PC1,col=pcaData$Population,pch=c(16,17)[as.numeric(pcaData$Cu_condition)],xlab="",ylab="",main="Shoot",cex=1.5,cex.lab=1.2,axes=F)
axis(1,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
axis(2,line=0,cex.axis=1,cex.lab=1,mgp=c(5,0.75,0),las=1)
mtext(side=1,line=2,paste0("PC1 (",percentVar[1],"%) "))
mtext(side=2,line=2,paste0("PC2 (",percentVar[2],"%) "))
box()
legend("topright",inset=c(-0.35,0),legend=c("Krom","Kosi"),fill=c("red","black"),title="Population",bty="n",cex=1.2)
legend("topright",inset=c(-0.35,0.3),legend=c(expression(paste(bold(L),"Cu")),expression(paste(bold(H),"Cu"))),pch=c(16,17),title="Treatment",bty="n",cex=1.2,pt.cex=1.5,title.adj=0)

dev.off()


#Calculate sample distance
sampleDists <- dist( t( assay(LogS) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(dds_shoot$Population, dds_shoot$Cu_condition, sep="-" )
colnames(sampleDistMatrix) <- NULL

attributes(sampleDists)$Labels<-paste(dds_shoot$Population, dds_shoot$Cu_condition, sep="-" )
hcluster<-hclust(sampleDists,method="ward.D") 
pdf("Hierarchical_clustering_Cu_rlog_shoot.pdf",width=4.72441,height=4.72441,paper="special")
plot(hcluster) 
dev.off()

require("RColorBrewer")
require("gplots")

pdf("Hierarchical_clustering_with_heatmap_Cu_rlog_shoot.pdf",width=4.72441,height=4.72441,paper="special")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(5,10))
dev.off()

sampleDists <- dist( t( assay(LogR) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(dds_root$Population, dds_root$Cu_condition, sep="-" )
colnames(sampleDistMatrix) <- NULL

attributes(sampleDists)$Labels<-paste(dds_root$Population, dds_root$Cu_condition, sep="-" )
hcluster<-hclust(sampleDists,method="ward.D") 
pdf("Hierarchical_clustering_Cu_rlog_root.pdf",width=4.72441,height=4.72441,paper="special")
plot(hcluster) 
dev.off()
pdf("Hierarchical_clustering_with_heatmap_Cu_rlog_root.pdf",width=4.72441,height=4.72441,paper="special")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(5,10))
dev.off()

sampleDists <- dist( t( assay(Log) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(dds$Population, dds$Cu_condition, dds$Tissue,sep="-" )
colnames(sampleDistMatrix) <- NULL

attributes(sampleDists)$Labels<-paste(dds$Population, dds$Cu_condition, dds$Tissue, sep="-" )
hcluster<-hclust(sampleDists,method="ward.D") 
pdf("Hierarchical_clustering_Cu_rlog_total.pdf",width=4.72441,height=4.72441,paper="special")
plot(hcluster) 
dev.off()
pdf("Hierarchical_clustering_with_heatmap_Cu_rlog_total.pdf",width=4.72441,height=4.72441,paper="special")
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins = c(5,10))
dev.off()



#Test PCA
plotPCA(LogR, intgroup=c("Population", "Cu_condition"))

pc<-prcomp(sampleDistMatrix,center=T,scale=T)
pc_val<-as.data.frame(predict(pc))
plot(pc_val$PC1,pc_val$PC2,col=colData(LogS)$Cu_condition,pch=as.numeric(colData(LogS)$Population),bg=colData(LogS)$Cu_condition)

text(pc_val[, 1:2],)
#Visualization of PCA
LogS <- rlog(dds_shoot, blind = FALSE)
LogR <- rlog(dds_root, blind = FALSE)

library("extrafont")
loadfonts(device="win")
windowsFonts(Arial=windowsFont("TT Arial")) 

pcaData <- plotPCA(LogR, intgroup=c("Population", "Cu_condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf('Root_PCA.pdf', width=10, height=10,paper="special")
RootPC <- ggplot(pcaData, aes(PC1, PC2, color=Population, shape=Cu_condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (",percentVar[1],"%) ")) +
  ylab(paste0("PC2 (",percentVar[2],"%) ")) + 
  coord_fixed() +
  theme_bw() +
  ggtitle("Root") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = pcaData$Population,
  values=c( "red", "blue")) +
  theme(plot.title = element_text(size = 14),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12)) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
  scale_shape_discrete(name="Treatment",labels=c("normal Cu","high Cu")) +
  theme(axis.title.x=element_text(element_text(size=12))) + 
  theme(plot.margin=unit(c(1,1,1,1),"cm")) +
  theme(axis.title.y=element_text(element_text(size=12)))
RootPC
dev.off()

pcaData <- plotPCA(LogS, intgroup=c("Population", "Cu_condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("PCA_shoot.png",type="cairo",units="in",width=8,height=8,pointsize=20,res=300)
ShootPC <- ggplot(pcaData, aes(PC1, PC2, color=Population, shape=Cu_condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (",percentVar[1],"%) ")) +
  ylab(paste0("PC2 (",percentVar[2],"%) ")) + 
  coord_fixed() +
  theme_bw() +
  ggtitle("Shoot") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = pcaData$Population,
  values=c( "red", "blue")) +
  theme(plot.title = element_text(size = 14),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12)) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
  scale_shape_discrete(name="Treatment",labels=c("normal Cu","high Cu")) +
  theme(axis.title.x=element_text(element_text(size=12))) +
  theme(axis.title.y=element_text(element_text(size=12)))
ShootPC
dev.off()

pcaData <- plotPCA(Log, intgroup=c("Population", "Cu_condition","Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("PCA_total.png",type="cairo",units="in",width=8,height=8,pointsize=20,res=300)
AllPC <- ggplot(pcaData, aes(PC1, PC2, color=Population, shape=Cu_condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1 (",percentVar[1],"%) ")) +
  ylab(paste0("PC2 (",percentVar[2],"%) ")) + 
  coord_fixed() +
  theme_bw() +
  ggtitle("All") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(breaks = pcaData$Population,
  values=c( "red", "blue")) +
  theme(plot.title = element_text(size = 14),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12)) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank()) +
  scale_shape_discrete(name="Treatment",labels=c("normal Cu","high Cu")) +
  theme(axis.title.x=element_text(element_text(size=12))) +
  theme(axis.title.y=element_text(element_text(size=12)))
AllPC
dev.off()

#hierarchical clustering with rlog
library("RColorBrewer")
library("pheatmap")
Log <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(Log)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- Log$Sample.Name
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Hierarchical_clustering_Cu_rlog.pdf",width=4.72441,height=4.72441,paper="special")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


#MapMan
# 2) Calculate cpm for cutoff for MapMan
  # cpm=(Counts/(sum mapped fragments/1 million)
CPM_df<-FPKM_df
for (i in 1:24)
	{Si<-data.frame(CPM_df[,2+i]/(sum(CPM_df[,2+i])/1000000))
	colnames(Si)<-paste("CPM_S",i,sep="")
	CPM_df<-data.frame(CPM_df,Si)
	}
Log_cpm<-log2(CPM_df[,51:74]+1)
colnames(Log_cpm)<-coldata$Sample.Name
rownames(Log_cpm)<-FPKM_df[,1]
Log_cpm<-Log_cpm[,c(19,11,3,20,12,8,24,16,4,21,13,5,17,9,1,18,10,2,23,15,7,22,14,6)]

Density_df<-apply(Log_cpm,2,density,na.rm=T)
plot(NA, xlim=range(sapply(Density_df, "[", "x")), ylim=range(sapply(Density_df, "[", "y")))
mapply(lines, Density_df, col=1:length(Density_df))
#or
library(affy)
plotDensity(Log_cpm, lty=1, lwd=2,  main="Log2-transformed counts per million (cpm)", xlab="log2(counts+1)", ylab="Number of genes",)
#legend("topright",pch=20, lty=1,lwd=0.75, bty="n", cex=0.75)
abline(v=1.2,col="red") 
#abline(v=1) 
#abline(v=1.5) 

#all with thaliana orthologue

Krom_R_Cu_Ath_desc_sorted_gff<-read.delim("Krom_root_Cu_reaction_RBH_medofratios_stricter.table",sep="\t")
Kosi_R_Cu_Ath_desc_sorted_gff<-read.delim("Kosi_root_Cu_reaction_RBH_medofratios_stricter.table",sep="\t")
Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff<-read.delim("Krom_Kosi_root_high_Cu_RBH_medofratios_stricter.table",sep="\t")
Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff<-read.delim("Krom_Kosi_root_normal_Cu_RBH_medofratios_stricter.table",sep="\t")

Kosi_S_Cu_Ath_desc_sorted_gff<-read.delim("Kosi_shoot_Cu_reaction_RBH_medofratios_stricter.table",sep="\t")
Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff<-read.delim("Krom_Kosi_shoot_high_Cu_RBH_medofratios_stricter.table",sep="\t")
Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff<-read.delim("Krom_Kosi_shoot_normal_Cu_RBH_medofratios_stricter.table",sep="\t")

names(Log_cpm)<-paste("CPM",1:24)
head(row.names(Log_cpm))

colnames(Krom_R_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Krom_R_Cu_cpm<-merge(Krom_R_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Krom_R_Cu_cpm_fil<-Krom_R_Cu_cpm[apply(Krom_R_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Krom_R_Cu_cpm_fil<-Krom_R_Cu_cpm_fil[!duplicated(Krom_R_Cu_cpm_fil$Ath_ID),]
Krom_R_Cu_cpm_Mapman<-cbind(as.character(Krom_R_Cu_cpm_fil$Ath_ID)[!is.na(Krom_R_Cu_cpm_fil$Ath_ID)],Krom_R_Cu_cpm_fil$log2FoldChange[!is.na(Krom_R_Cu_cpm_fil$Ath_ID)])
Krom_R_Cu_cpm_Mapman<-Krom_R_Cu_cpm_Mapman[!duplicated(Krom_R_Cu_cpm_Mapman[,1]),]
write.table(Krom_R_Cu_cpm_Mapman,"Krom_R_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

colnames(Kosi_R_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Kosi_R_Cu_cpm<-merge(Kosi_R_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Kosi_R_Cu_cpm_fil<-Kosi_R_Cu_cpm[apply(Kosi_R_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Kosi_R_Cu_cpm_fil<-Kosi_R_Cu_cpm_fil[!duplicated(Kosi_R_Cu_cpm_fil$Ath_ID),]
Kosi_R_Cu_cpm_Mapman<-cbind(as.character(Kosi_R_Cu_cpm_fil$Ath_ID)[!is.na(Kosi_R_Cu_cpm_fil$Ath_ID)],Kosi_R_Cu_cpm_fil$log2FoldChange[!is.na(Kosi_R_Cu_cpm_fil$Ath_ID)])
Kosi_R_Cu_cpm_Mapman<-Kosi_R_Cu_cpm_Mapman[!duplicated(Kosi_R_Cu_cpm_Mapman[,1]),]
write.table(Kosi_R_Cu_cpm_Mapman,"Kosi_R_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

colnames(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Krom_Kosi_R_high_Cu_cpm<-merge(Krom_Kosi_R_high_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Krom_Kosi_R_high_Cu_cpm_fil<-Krom_Kosi_R_high_Cu_cpm[apply(Krom_Kosi_R_high_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Krom_Kosi_R_high_Cu_cpm_fil<-Krom_Kosi_R_high_Cu_cpm_fil[!duplicated(Krom_Kosi_R_high_Cu_cpm_fil$Ath_ID),]
Krom_Kosi_R_high_Cu_cpm_Mapman<-cbind(as.character(Krom_Kosi_R_high_Cu_cpm_fil$Ath_ID)[!is.na(Krom_Kosi_R_high_Cu_cpm_fil$Ath_ID)],Krom_Kosi_R_high_Cu_cpm_fil$log2FoldChange[!is.na(Krom_Kosi_R_high_Cu_cpm_fil$Ath_ID)])
Krom_Kosi_R_high_Cu_cpm_Mapman<-Krom_Kosi_R_high_Cu_cpm_Mapman[!duplicated(Krom_Kosi_R_high_Cu_cpm_Mapman[,1]),]
write.table(Krom_Kosi_R_high_Cu_cpm_Mapman,"Krom_Kosi_R_high_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

colnames(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Krom_Kosi_R_normal_Cu_cpm<-merge(Krom_Kosi_R_normal_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Krom_Kosi_R_normal_Cu_cpm_fil<-Krom_Kosi_R_normal_Cu_cpm[apply(Krom_Kosi_R_normal_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Krom_Kosi_R_normal_Cu_cpm_fil<-Krom_Kosi_R_normal_Cu_cpm_fil[!duplicated(Krom_Kosi_R_normal_Cu_cpm_fil$Ath_ID),]
Krom_Kosi_R_normal_Cu_cpm_Mapman<-cbind(as.character(Krom_Kosi_R_normal_Cu_cpm_fil$Ath_ID)[!is.na(Krom_Kosi_R_normal_Cu_cpm_fil$Ath_ID)],Krom_Kosi_R_normal_Cu_cpm_fil$log2FoldChange[!is.na(Krom_Kosi_R_normal_Cu_cpm_fil$Ath_ID)])
Krom_Kosi_R_normal_Cu_cpm_Mapman<-Krom_Kosi_R_normal_Cu_cpm_Mapman[!duplicated(Krom_Kosi_R_normal_Cu_cpm_Mapman[,1]),]
write.table(Krom_Kosi_R_normal_Cu_cpm_Mapman,"Krom_Kosi_R_normal_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

colnames(Kosi_S_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Kosi_S_Cu_cpm<-merge(Kosi_S_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Kosi_S_Cu_cpm_fil<-Kosi_S_Cu_cpm[apply(Kosi_S_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Kosi_S_Cu_cpm_fil<-Kosi_S_Cu_cpm_fil[!duplicated(Kosi_S_Cu_cpm_fil$Ath_ID),]
Kosi_S_Cu_cpm_Mapman<-cbind(as.character(Kosi_S_Cu_cpm_fil$Ath_ID)[!is.na(Kosi_S_Cu_cpm_fil$Ath_ID)],Kosi_S_Cu_cpm_fil$log2FoldChange[!is.na(Kosi_S_Cu_cpm_fil$Ath_ID)])
Kosi_S_Cu_cpm_Mapman<-Kosi_S_Cu_cpm_Mapman[!duplicated(Kosi_S_Cu_cpm_Mapman[,1]),]
write.table(Kosi_S_Cu_cpm_Mapman,"Kosi_S_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

colnames(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Krom_Kosi_S_high_Cu_cpm<-merge(Krom_Kosi_S_high_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Krom_Kosi_S_high_Cu_cpm_fil<-Krom_Kosi_S_high_Cu_cpm[apply(Krom_Kosi_S_high_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Krom_Kosi_S_high_Cu_cpm<-Krom_Kosi_S_high_Cu_cpm[!duplicated(Krom_Kosi_S_high_Cu_cpm$Ath_ID),]
Krom_Kosi_S_high_Cu_cpm_Mapman<-cbind(as.character(Krom_Kosi_S_high_Cu_cpm_fil$Ath_ID)[!is.na(Krom_Kosi_S_high_Cu_cpm_fil$Ath_ID)],Krom_Kosi_S_high_Cu_cpm_fil$log2FoldChange[!is.na(Krom_Kosi_S_high_Cu_cpm_fil$Ath_ID)])
Krom_Kosi_S_high_Cu_cpm_Mapman<-Krom_Kosi_S_high_Cu_cpm_Mapman[!duplicated(Krom_Kosi_S_high_Cu_cpm_Mapman[,1]),]
write.table(Krom_Kosi_S_high_Cu_cpm_Mapman,"Krom_Kosi_S_high_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

colnames(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff)[1]<-"Ah_gene"
Krom_Kosi_S_normal_Cu_cpm<-merge(Krom_Kosi_S_normal_Cu_Ath_desc_sorted_gff,Log_cpm,by.y="row.names",by.x="Ah_gene")
Krom_Kosi_S_normal_Cu_cpm_fil<-Krom_Kosi_S_normal_Cu_cpm[apply(Krom_Kosi_S_normal_Cu_cpm[,93:116],1,function(x) any(x>1.2)),]
Krom_Kosi_S_normal_Cu_cpm_fil<-Krom_Kosi_S_normal_Cu_cpm_fil[!duplicated(Krom_Kosi_S_normal_Cu_cpm_fil$Ath_ID),]
Krom_Kosi_S_normal_Cu_cpm_Mapman<-cbind(as.character(Krom_Kosi_S_normal_Cu_cpm_fil$Ath_ID)[!is.na(Krom_Kosi_S_normal_Cu_cpm_fil$Ath_ID)],Krom_Kosi_S_normal_Cu_cpm_fil$log2FoldChange[!is.na(Krom_Kosi_S_normal_Cu_cpm_fil$Ath_ID)])
Krom_Kosi_S_normal_Cu_cpm_Mapman<-Krom_Kosi_S_normal_Cu_cpm_Mapman[!duplicated(Krom_Kosi_S_normal_Cu_cpm_Mapman[,1]),]
write.table(Krom_Kosi_S_normal_Cu_cpm_Mapman,"Krom_Kosi_S_normal_Cu_Mapman_stricter_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

#whole dataset of expressed genes with cpm > 1.2 for Mapman
write.table(Log_cpm,"Log_cpm.table",sep="\t",row.names=T)

OG_Ahalleri<-read.table("Ortho_desc.table",sep="\t",header=T,quote="",fill=T)

#root
All<-merge(OG_Ahalleri,Log_cpm,by.y="row.names",by.x="Ah_ID")
All_cpm_R<-All[apply(All[,10:33],1,function(x) any(x>1.2)),]
All_cpm_R<-All_cpm_R[!duplicated(All_cpm_R$Ath_ID),]
All_cpm_Mapman_R<-cbind(as.character(All_cpm_R$Ath_ID)[!is.na(All_cpm_R$Ath_ID)],All_cpm_R[!is.na(All_cpm_R$Ath_ID),10])
All_cpm_Mapman_R<-All_cpm_Mapman_R[!duplicated(All_cpm_Mapman_R[,1]),]
write.table(All_cpm_Mapman_R,"All_cpm_Mapman_stricter_R_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

#shoot
All<-merge(OG_Ahalleri,Log_cpm,by.y="row.names",by.x="Ah_ID")
All_cpm_S<-All[apply(All[,10:33],1,function(x) any(x>1.2)),]
All_cpm_S<-All_cpm_S[!duplicated(All_cpm_S$Ath_ID),]
All_cpm_Mapman_S<-cbind(as.character(All_cpm_S$Ath_ID)[!is.na(All_cpm_S$Ath_ID)],All_cpm_S[!is.na(All_cpm_S$Ath_ID),10])
All_cpm_Mapman_S<-All_cpm_Mapman_S[!duplicated(All_cpm_Mapman_S[,1]),]
write.table(All_cpm_Mapman_S,"All_cpm_Mapman_stricter_S_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

Krom_Kosi_R_normal_Cu_cpm_Mapman_up<-Krom_Kosi_R_normal_Cu_cpm_Mapman[Krom_Kosi_R_normal_Cu_cpm_Mapman[,2]>0,]
Krom_Kosi_R_normal_Cu_cpm_Mapman_down<-Krom_Kosi_R_normal_Cu_cpm_Mapman[Krom_Kosi_R_normal_Cu_cpm_Mapman[,2]<0,]
write.table(Krom_Kosi_R_normal_Cu_cpm_Mapman_up,"Krom_Kosi_R_normal_Cu_Mapman_stricter_up_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(Krom_Kosi_R_normal_Cu_cpm_Mapman_down,"Krom_Kosi_R_normal_Cu_Mapman_stricter_down_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

Krom_Kosi_S_normal_Cu_cpm_Mapman_up<-Krom_Kosi_S_normal_Cu_cpm_Mapman[Krom_Kosi_S_normal_Cu_cpm_Mapman[,2]>0,]
Krom_Kosi_S_normal_Cu_cpm_Mapman_down<-Krom_Kosi_S_normal_Cu_cpm_Mapman[Krom_Kosi_S_normal_Cu_cpm_Mapman[,2]<0,]
write.table(Krom_Kosi_S_normal_Cu_cpm_Mapman_up,"Krom_Kosi_S_normal_Cu_Mapman_stricter_up_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(Krom_Kosi_S_normal_Cu_cpm_Mapman_down,"Krom_Kosi_S_normal_Cu_Mapman_stricter_down_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

Krom_Kosi_R_high_Cu_cpm_Mapman_up<-Krom_Kosi_R_high_Cu_cpm_Mapman[Krom_Kosi_R_high_Cu_cpm_Mapman[,2]>0,]
Krom_Kosi_R_high_Cu_cpm_Mapman_down<-Krom_Kosi_R_high_Cu_cpm_Mapman[Krom_Kosi_R_high_Cu_cpm_Mapman[,2]<0,]
write.table(Krom_Kosi_R_high_Cu_cpm_Mapman_up,"Krom_Kosi_R_high_Cu_Mapman_stricter_up_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(Krom_Kosi_R_high_Cu_cpm_Mapman_down,"Krom_Kosi_R_high_Cu_Mapman_stricter_down_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)

Krom_Kosi_S_high_Cu_cpm_Mapman_up<-Krom_Kosi_S_high_Cu_cpm_Mapman[Krom_Kosi_S_high_Cu_cpm_Mapman[,2]>0,]
Krom_Kosi_S_high_Cu_cpm_Mapman_down<-Krom_Kosi_S_high_Cu_cpm_Mapman[Krom_Kosi_S_high_Cu_cpm_Mapman[,2]<0,]
write.table(Krom_Kosi_S_high_Cu_cpm_Mapman_up,"Krom_Kosi_S_high_Cu_Mapman_stricter_up_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)
write.table(Krom_Kosi_S_high_Cu_cpm_Mapman_down,"Krom_Kosi_S_high_Cu_Mapman_stricter_down_dedup.txt",sep="\t",row.names=F,col.names=F,quote=F)


#Volcano plots
```{r}
###Enhanced version of volcano plots for pubication

#if (!requireNamespace("BiocManager", quietly = TRUE)) 
#install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")
#BiocManager::install('EnhancedVolcano')
#install.packages("ggrepel")
library(ggrepel)
library(EnhancedVolcano)
library(magrittr)

library(gridExtra)
library(grid)
library(cowplot)
library(gridExtra)
library(gtable)



##Modified version based on source for EnhancedVolcano

#Normalized by all samples
EnhancedVolcano_GJL <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  xlim = c(min(toptable[,x], na.rm=TRUE),
    max(toptable[,x], na.rm=TRUE)),
  ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 12,
  pCutoff = 10e-6,
  pLabellingCutoff = pCutoff,
  FCcutoff = 1,
  title = "",
  titleLabSize = 16,
  transcriptPointSize = 0.8,
  transcriptLabSize = 3.0,
  labhjust = 0,
  labvjust = 1.5,
  col = c("grey30",  "royalblue", "red2"),
  colOverride = NULL,
  colAlpha = 1/2,
  legend = c("NS","DOWN","UP"),
  legendPosition = "top",
  legendLabSize = 10,
  legendIconSize = 3.0,
  legendVisible = TRUE,
  shade = NULL,
  shadeLabel = NULL,
  shadeAlpha = 1/2,
  shadeFill = "grey",
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colConnectors = "black",
  cutoffLineType = "longdash",
  cutoffLineCol = "black",
  cutoffLineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  border = "partial",
  borderWidth = 0.8,
  borderColour = "black")
{
  if(!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call.=FALSE)
  }

  if(!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call.=FALSE)
  }

  if(!is.numeric(toptable[,x])) {
    stop(paste(x, " is not numeric!", sep=""))
  }

  if(!is.numeric(toptable[,y])) {
    stop(paste(y, " is not numeric!", sep=""))
  }

  i <- xvals <- yvals <- Sig <- NULL

  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(toptable[,y]<pCutoff) &
    (toptable[,x] >FCcutoff)] <- "FC_P_UP"
  toptable$Sig[(toptable[,y]<pCutoff) &
    (toptable[,x] < -FCcutoff)] <- "FC_P_DOWN"
  toptable$Sig <- factor(toptable$Sig,
    levels=c("NS","FC_P_DOWN", "FC_P_UP"))

  # some software programs return 0 for very low p-values
  # These throw an error in EnhancedVolcano
  # Detect these, issue warning, and convert these to
  # machine-lowest value possible
  if (min(toptable[,y], na.rm=TRUE) == 0) {
    warning(paste("One or more P values is 0.",
      "Converting to minimum possible value..."),
      call. = FALSE)
    toptable[which(toptable[,y] == 0), y] <- .Machine$double.xmin
  }

  toptable$lab <- lab
  toptable$xvals <- toptable[,x]
  toptable$yvals <- toptable[,y]

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size=14) +

    theme(
      legend.background=element_rect(),
      plot.title=element_text(angle=0, size=16, face="bold", vjust=0.5, hjust=0.5),

      axis.text.x=element_text(angle=0, size=13, vjust=1),
      axis.text.y=element_text(angle=0, size=13, vjust=1),
      axis.title=element_text(size=axisLabSize),

      legend.position=legendPosition,
      legend.key=element_blank(),
      legend.key.size=unit(0.7, "cm"),
      legend.text=element_text(size=legendLabSize),

      title=element_text(size=legendLabSize),
      legend.title=element_blank())

  # Create the plot object differently based on whether colOverride is
  # NULL or not. This helps to avoid messing up the legend.
  if (!is.null(colOverride)) {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      guides(colour = guide_legend(order = 1,
        override.aes=list(size=legendIconSize))) +

      geom_point(aes(color=factor(names(colOverride))),
        alpha=colAlpha,
        size=transcriptPointSize,
        na.rm = TRUE) +

      scale_color_manual(values=colOverride)
  } else {
    plot <- ggplot(toptable, aes(x=xvals, y=-log10(yvals))) + th +

      guides(colour = guide_legend(order = 1,
        override.aes=list(size=legendIconSize))) +

      geom_point(aes(color=factor(Sig)),
        alpha=colAlpha,
        size=transcriptPointSize,
        na.rm = TRUE,
        show.legend = legendVisible) +

      scale_color_manual(values=c(NS=col[1],
        FC_P_DOWN=col[2],
        FC_P_UP = col[3]),
        labels=c(NS=legend[1],
        FC_P_DOWN=paste(legend[2], sep=""),
        FC_P_UP=paste(legend[3], sep="")))
  }

  # add more elements to the plot
  plot <- plot +

    xlab(xlab) +
    ylab(ylab) +

    xlim(xlim[1], xlim[2]) +
    ylim(ylim[1], ylim[2]) +

    ggtitle(title) +

    geom_vline(xintercept=c(-FCcutoff, FCcutoff),
      linetype=cutoffLineType,
      colour=cutoffLineCol,
      size=cutoffLineWidth) +

    geom_hline(yintercept=-log10(pCutoff),
      linetype=cutoffLineType,
      colour=cutoffLineCol,
      size=cutoffLineWidth)

  # Border around plot
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, fill = NA, size = borderWidth))
  } else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, colour = borderColour),
      panel.border = element_blank(),
      panel.background = element_blank())
  } else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }

  # Gridlines
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # For labeling with geom_text_repel (connectors) and
  # geom_text(.., check_overlap = TRUE), 4 possible scenarios
  # can arise
  if (drawConnectors == TRUE && is.null(selectLab)) {
    plot <- plot + geom_text_repel(
      data=subset(toptable,
        toptable[,y]<pLabellingCutoff &
          abs(toptable[,x])>FCcutoff),
        aes(label=subset(toptable,
          toptable[,y]<pLabellingCutoff &
            abs(toptable[,x])>FCcutoff)[,"lab"]),
        size = transcriptLabSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        hjust = labhjust,
        vjust = labvjust,
        na.rm = TRUE)
  } else if (drawConnectors == TRUE && !is.null(selectLab)) {
    plot <- plot + geom_text_repel(
      data=subset(toptable,
        !is.na(toptable[,"lab"])),
        aes(label=subset(toptable,
          !is.na(toptable[,"lab"]))[,"lab"]),
        size = transcriptLabSize,
        segment.color = colConnectors,
        segment.size = widthConnectors,
        hjust = labhjust,
        vjust = labvjust,
        na.rm = TRUE)
  } else if (drawConnectors == FALSE && !is.null(selectLab)) {
    plot <- plot + geom_text(
      data=subset(toptable,
        !is.na(toptable[,"lab"])),
        aes(label=subset(toptable,
          !is.na(toptable[,"lab"]))[,"lab"]),
        size = transcriptLabSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        na.rm = TRUE)
  } else if (drawConnectors == FALSE && is.null(selectLab)) {
    plot <- plot + geom_text(
      data=subset(toptable,
        toptable[,y]<pLabellingCutoff &
          abs(toptable[,x])>FCcutoff),
        aes(label=subset(toptable,
          toptable[,y]<pLabellingCutoff &
            abs(toptable[,x])>FCcutoff)[,"lab"]),
        size = transcriptLabSize,
        check_overlap = TRUE,
        hjust = labhjust,
        vjust = labvjust,
        na.rm = TRUE)
  }

  # shading
  if (!is.null(shade)) {
    plot <- plot + 
      stat_density2d(
        data = subset(toptable,
          rownames(toptable) %in% shade),
        fill = shadeFill,
        alpha = shadeAlpha,
        geom = "polygon",
        contour = TRUE,
        size = shadeSize,
        bins = shadeBins,
        show.legend = FALSE,
        na.rm = TRUE) +

      scale_fill_identity(name = shadeLabel,
        labels = shadeLabel,
        guide = "legend")
  }

  return(plot)
}



KromR <- EnhancedVolcano_GJL(IA_R_Krom, 
                lab = rownames(IA_R_Krom), 
                x = "log2FoldChange",
                y = "padj",
                selectLab = "",
                xlab = "",
                ylab = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10,10),
                ylim = c(0,10),
                title = "Krom Cu root",
                colAlpha = 0.4,
                legend=c("NS","Down-regulated","Up-regulated"),
                legendPosition = "none",
                legendLabSize = 10,
                legendIconSize = 3.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,  
                border = "full",
                borderWidth = 1.2,
                borderColour = "grey30",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0
                )

KosiR <- EnhancedVolcano_GJL(IA_R_Kosi, 
                lab = rownames(IA_R_Kosi), 
                x = "log2FoldChange",
                y = "padj",
                selectLab = "",
                xlab = "",
                ylab = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10,10),
                ylim = c(0,10),
                title = "Kosi Cu root",
                colAlpha = 0.4,
                legend=c("NS","Down-regulated", "Up-regulated"),
                legendPosition = "none",
                legendLabSize = 10,
                legendIconSize = 3.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,  
                border = "full",
                borderWidth = 1.2,
                borderColour = "grey30",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0
                )
KromS <- EnhancedVolcano_GJL(IA_S_Krom, 
                lab = rownames(IA_S_Krom), 
                x = "log2FoldChange",
                y = "padj",
                selectLab = "",
                xlab = "",
                ylab = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10,10),
                ylim = c(0,10),
                title = "Krom Cu shoot",
                colAlpha = 0.4,
                legend=c("NS","Down-regulated", "Up-regulated"),
                legendPosition = "none",
                legendLabSize = 10,
                legendIconSize = 3.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,  
                border = "full",
                borderWidth = 1.2,
                borderColour = "grey30",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0
                )

KosiS <- EnhancedVolcano_GJL(IA_S_Kosi, 
                lab = rownames(IA_S_Kosi), 
                x = "log2FoldChange",
                y = "padj",
                selectLab = "",
                xlab = "",
                ylab = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10,10),
                ylim = c(0,10),
                title = "Kosi Cu shoot",
                colAlpha = 0.4,
                legend=c("NS","Down-regulated", "Up-regulated"),
                legendPosition = "none",
                legendLabSize = 10,
                legendIconSize = 3.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,  
                border = "full",
                borderWidth = 1.2,
                borderColour = "grey30",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0
                )


pdf("Volcano_plots_Cu_reaction_stricter.pdf",width=9.44882,height=9.44882,paper="special")

legend_vol <- get_legend(KosiR + theme(legend.position="bottom"))
prow <- grid.arrange(KromR,KosiR,KromS,KosiS, ncol=2, nrow=2,
                     bottom=textGrob(bquote(~Log[2]~ "fold change"), gp = gpar(cex = 1.3), vjust=0),
                     left=textGrob(bquote(~-Log[10]~adjusted~italic(P)), rot = 90, vjust=1, gp = gpar(cex = 1.3)))
plot_grid( prow, legend_vol, ncol = 1, rel_heights = c(1, .1))             

dev.off()

library("RColorBrewer")
library("pheatmap")
#LogR <- rlog(dds_root, blind = FALSE)
#LogS <- rlog(dds_shoot, blind = FALSE)
Krom_R_Cu_df<-IA_R_Krom[rownames(IA_R_Krom)%in%Krom_R_Cu_Ath_desc_sorted_gff$Ah_gene,]
Krom_R_Cu_mean<-cbind(rownames(Krom_R_Cu_df),mcols(Krom_R_Cu_df)$baseMean)

sampleDistsR <- dist(t(assay(LogR)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- Log$Sample.Name
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("Hierarchical_clustering_Cu_rlog.pdf",width=4.72441,height=4.72441,paper="special")
pheatmap(as.matrix(data.frame()),
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()
#biocLite("genefilter")
library( "genefilter" )
library( "gplots" )
library( "RColorBrewer" )

# Heatmaps are more comparable if the expression data has been transformed. We can transform our data using the function rlog()
rld <- rlog( dds_root )
colnames(rld)<-dds_root$Sample.Name
# We then generate a heatmap by assaying the rld data based on our adjusted p values. Other factors could also be used to subet the gene list.

Krom_Cu_Ah<-merge(IA_R_Krom_sig,countdata_rlog_df,by.x="row.names",by.y="seq_name",all.x=TRUE)
Krom_R_Cu_Ah<-Krom_Cu_Ah[,c(1:19)]
Krom_R_Cu_Ah_ordered<-Krom_R_Cu_Ah[,c("Row.names","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","Krom -Cu R1","Krom -Cu R2","Krom -Cu R3","Krom +Cu R1","Krom +Cu R2","Krom +Cu R3","Kosi -Cu R1","Kosi -Cu R2","Kosi -Cu R3","Kosi +Cu R1","Kosi +Cu R2","Kosi +Cu R3")]
colnames(Krom_R_Cu_Ah_ordered)[8:19]<-c("Krom control R1","Krom control R2","Krom control R3","Krom high Cu R1","Krom high Cu R2","Krom high Cu R3","Kosi control R1","Kosi control R2","Kosi control R3","Kosi high Cu R1","Kosi high Cu R2","Kosi high Cu R3")
pheatmap(Krom_R_Cu_Ah_ordered[,8:19], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(bold(Krom~root~Cu~reaction)))


Cureac1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists_M,Kosi_R_Cu_Ath_desc_sorted_gff_lists_M)
Cureac<-Cureac1[!duplicated(Cureac1$Ah_gene),]
pheatmap(Cureac[,10:21], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(bold(Cu~reaction)))

CureacS<-Kosi_S_Cu_Ath_desc_sorted_gff_lists_M[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists_M$Ah_gene),]
pheatmap(CureacS[,22:33], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(bold(Cu~reaction)),treeheight_row = 0, treeheight_col = 0)




Cureac1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists_M,Kosi_R_Cu_Ath_desc_sorted_gff_lists_M)
Cureac<-Cureac1[!duplicated(Cureac1$Ah_gene),]
rownames(Cureac)<-c("COPT6","FSD1","FRO4_1","FRO4_2","CITF1","CCH","YSL2","ZIP2")
R<-pheatmap(Cureac[,c(19:21,16:18,13:15,10:12)], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Root),treeheight_ro =0,treeheight_col=0,gaps_col=c(3,6,9),show_colnames=F,fontsize_row=10,fontsize_col=10,cellwidth=35,cellheight=35)$gtable

CureacS<-Kosi_S_Cu_Ath_desc_sorted_gff_lists_M[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists_M$Ah_gene),]
rownames(CureacS)<-c("CCS","CV","CSD2","RAD50","AT3G49160")
labels_col<-c("","Kosi -Cu","","","Kosi +Cu","","","Krom -Cu","","","Krom +Cu","")
S<-pheatmap(CureacS[,c(31:33,28:30,25:27,22:24)], scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Shoot),treeheight_ro =0,treeheight_col=0,gaps_col=c(3,6,9),labels_col=labels_col,angle_col=0,fontsize_row=10,fontsize_col=10,cellwidth=35,cellheight=35)$gtable

require(gridExtra)
library(gtable)
library(grid)
g <- rbind(R, S,size="first")
g$widths <- unit.pmax(R$widths, S$widths)
grid.newpage()
pdf("Heatmap_Cu_reaction_stricter.pdf",width=9.44882,height=9.44882,paper="special")
grid.draw(g)
dev.off()


#mean and side by side
Cureac1<-rbind(Krom_R_Cu_Ath_desc_sorted_gff_lists_M,Kosi_R_Cu_Ath_desc_sorted_gff_lists_M)
Cureac2<-Cureac1[!duplicated(Cureac1$Ah_gene),]
Cureac_KosiCtrl<-apply(Cureac2[,c(19:21)],1,mean)
Cureac_KosiCu<-apply(Cureac2[,c(16:18)],1,mean)
Cureac_KromCtrl<-apply(Cureac2[,c(13:15)],1,mean)
Cureac_KromCu<-apply(Cureac2[,c(10:12)],1,mean)
Cureac<-cbind(Cureac_KosiCtrl,Cureac_KosiCu,Cureac_KromCtrl,Cureac_KromCu)
rownames(Cureac)<-c("COPT2","FSD1","FRO4","FRO5","CITF1","CCH","YSL2","ZIP2")
labels_col<-c("",expression(paste("Kosi ",bold(L),"Cu")),"","",expression(paste("Kosi ",bold(H),"Cu")),"","",expression(paste("Krom ",bold(L),"Cu")),"","",expression(paste("Krom ",bold(H),"Cu")),"")
labels_row<-c(expression(italic(COPT2^"*")),expression(italic(FSD1^"*")),expression(italic(FRO4^"*")),expression(italic(FRO5^"*")),expression(italic(CITF1)),expression(italic(CCH)),expression(italic(YSL2)),expression(italic(ZIP2)))

R<-pheatmap(Cureac, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(Root),treeheight_row =0,treeheight_col=0,show_colnames=T,angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,legend=F,labels_col=labels_col,labels_row=labels_row,fontsize=12)$gtable

CureacS1<-Kosi_S_Cu_Ath_desc_sorted_gff_lists_M[!duplicated(Kosi_S_Cu_Ath_desc_sorted_gff_lists_M$Ah_gene),]
CureacS_KosiCtrl<-apply(CureacS1[,c(31:33)],1,mean)
CureacS_KosiCu<-apply(CureacS1[,c(28:30)],1,mean)
CureacS_KromCtrl<-apply(CureacS1[,c(25:27)],1,mean)
CureacS_KromCu<-apply(CureacS1[,c(22:24)],1,mean)
CureacS<-cbind(CureacS_KosiCtrl,CureacS_KosiCu,CureacS_KromCtrl,CureacS_KromCu)
rownames(CureacS)<-c("CCS","CSD2","CV","RAD50","AT3G49160")
labels_row<-c(expression(italic(CCS)),expression(italic(CSD2)),expression(italic(CV)),expression(italic(RAD50)),expression(italic(AT3G49160)))
labels_col<-c("",expression(paste("Kosi ",bold(L),"Cu")),"","",expression(paste("Kosi ",bold(H),"Cu")),"","",expression(paste("Krom ",bold(L),"Cu")),"","",expression(paste("Krom ",bold(H),"Cu")),"")
CureacS<-CureacS[c(1,3,2,4,5),]
S<-pheatmap(CureacS, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(Shoot),treeheight_row=0,treeheight_col=0,labels_col=labels_col,angle_col=45,fontsize_row=12,fontsize_col=12,cellwidth=20,cellheight=20,labels_row=labels_row,fontsize=12)$gtable

require(gridExtra)
library(gtable)
library(grid)
g <- cbind(R, S,size="first")
g$heights <- unit.pmax(R$heights, S$heights)

grid.newpage()
pdf("Heatmap_Cu_reaction_stricter_beside_smaller.pdf",width=6,height=6,paper="special")
grid.arrange(R,S,ncol=2)
dev.off()



LogR <- rlog(dds_root, blind = FALSE)
LogS <- rlog(dds_shoot, blind = FALSE)

sampleDistsR <- assay(LogR)
sampleDistMatrixR <- as.matrix(sampleDistsR)
colnames(sampleDistMatrixR) <-LogR$Sample.Name
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDistMatrixR<-sampleDistMatrixR[,c(11,7,3,12,8,2,10,6,4,9,5,1)]
CandR<-sampleDistMatrixR[rownames(sampleDistMatrixR)%in%Cureac$Ah_gene,]
pheatmap(CandR, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=F, cluster_cols=F,main=expression(bold(Cu~reaction)),treeheight_row = 0, treeheight_col = 0)

par(oma=c(3,1,1,1))
sampleDistsS <- assay(LogS)
sampleDistMatrixS <- as.matrix(sampleDistsS)
colnames(sampleDistMatrixS) <-LogS$Sample.Name
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDistMatrixS<-sampleDistMatrixS[,c(11,7,3,12,8,4,10,6,2,9,5,1)]
CandS<-sampleDistMatrixS[rownames(sampleDistMatrixS)%in%CureacS$Ah_gene,]
pheatmap(CandS, scale="row",col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),cluster_rows=T, cluster_cols=F,main=expression(bold(Cu~reaction)),treeheight_row = 0, treeheight_col = 0)













