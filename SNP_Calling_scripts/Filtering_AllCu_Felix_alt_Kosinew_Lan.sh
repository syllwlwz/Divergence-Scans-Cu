#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_allCu_Kosinew_LanLan
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#Filter 1: select bialleleic SNPs without missing genotypes
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_Felix/All_Cu_Kosinew_Lan.vcf -R Ahalleri_Felix_renamed_sorted.fasta -selectType SNP --excludeNonVariants -restrictAllelesTo BIALLELIC -env -o Filtered_Felix/All_Cu_Kosinew_LanFil1.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil1.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix/All_Cu_Kosinew_LanFil1.vcf -R Ahalleri_Felix_renamed_sorted.fasta -G_filter "DP<2.0" -G_filterName "lowDP" -o Filtered_Felix/All_Cu_Kosinew_LanFil2.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil2.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/All_Cu_Kosinew_LanFil2.vcf -R Ahalleri_Felix_renamed_sorted.fasta --setFilteredGtToNocall -o Filtered_Felix/All_Cu_Kosinew_LanFil3.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil3.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V Filtered_Felix/All_Cu_Kosinew_LanFil3.vcf -R Ahalleri_Felix_renamed_sorted.fasta -select "AN >= 172" -o Filtered_Felix/All_Cu_Kosinew_LanFil4.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil4.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix/All_Cu_Kosinew_LanFil4.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filter "QD<2.0 || FS>40.0 || MQ<50.0 || MQRankSum< -2.5 || ReadPosRankSum < -4.0 || SOR>4.0" -filterName "BP" -o Filtered_Felix/All_Cu_Kosinew_LanFil5a.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil5a.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/All_Cu_Kosinew_LanFil5a.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -o Filtered_Felix/All_Cu_Kosinew_LanFil5.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil5.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/All_Cu_Kosinew_LanFil5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F DP -raw --allowMissingData -o Filtered_Felix/All_Cu_Kosinew_LanFil5.table 2> Filtered_Felix/All_Cu_Kosinew_LanFil5_DPtable.err

grep -v "DP" Filtered_Felix/All_Cu_Kosinew_LanFil5.table > Filtered_Felix/All_Cu_Kosinew_LanFil5_DP_hQD.table
mode=$(./software/R-3.5.3/bin/Rscript -e 'library("modes"); data<-read.table("Filtered_Felix/All_Cu_Kosinew_LanFil5_DP_hQD.table",header=F); cat(modes(data$V1)[1])')
DPmax=2*$mode
DPmin=0.5*$mode

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix/All_Cu_Kosinew_LanFil5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filterName "DP" -filter "DP<"$DPmin" || DP>"$DPmax"" -o Filtered_Felix/All_Cu_Kosinew_LanFil6.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil6.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/All_Cu_Kosinew_LanFil6.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -o Filtered_Felix/All_Cu_Kosinew_LanFil7.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil7.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/All_Cu_Kosinew_LanFil7.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -select "AF<1.0 && AF>0.0" -o Filtered_Felix/All_Cu_Kosinew_LanFil8.vcf 2> Filtered_Felix/All_Cu_Kosinew_LanFil8.err

#Extract tables
#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/All_Cu_Kosinew_LanFil8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/All_Cu_Kosinew_LanFil.table 2> Filtered_Felix/All_Cu_Kosinew_LanFil2.err

