#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_LangBest_UG
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#Filter 1: select bialleleic SNPs without missing genotypes
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V UG_Felix/LangBest_UG.vcf.gz -R Ahalleri_Felix_renamed_sorted.fasta -sf LangBest.list -selectType SNP --excludeNonVariants -restrictAllelesTo BIALLELIC -env -o Filtered_Felix_UG_test/LangBestFil1.vcf 2> Filtered_Felix_UG_test/LangBestFil1.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix_UG_test/LangBestFil1.vcf -R Ahalleri_Felix_renamed_sorted.fasta -G_filter "DP<2.0" -G_filterName "lowDP" -o Filtered_Felix_UG_test/LangBestFil2.vcf 2> Filtered_Felix_UG_test/LangBestFil2.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix_UG_test/LangBestFil2.vcf -R Ahalleri_Felix_renamed_sorted.fasta --setFilteredGtToNocall -o Filtered_Felix_UG_test/LangBestFil3.vcf 2> Filtered_Felix_UG_test/LangBestFil3.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V Filtered_Felix_UG_test/LangBestFil3.vcf -R Ahalleri_Felix_renamed_sorted.fasta -select "AN >= 40" -o Filtered_Felix_UG_test/LangBestFil4.vcf 2> Filtered_Felix_UG_test/LangBestFil4.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix_UG_test/LangBestFil4.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filter "QD<2.0 || FS>40.0 || MQ<50.0 || MQRankSum< -2.5 || ReadPosRankSum < -4.0 || SOR>4.0" -filterName "BP" -o Filtered_Felix_UG_test/LangBestFil5a.vcf 2> Filtered_Felix_UG_test/LangBestFil5a.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix_UG_test/LangBestFil5a.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -o Filtered_Felix_UG_test/LangBestFil5.vcf 2> Filtered_Felix_UG_test/LangBestFil5.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix_UG_test/LangBestFil5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F DP -raw --allowMissingData -o Filtered_Felix_UG_test/LangBestFil5.table 2> Filtered_Felix_UG_test/LangBestFil5_DPtable.err

grep -v "DP" Filtered_Felix_UG_test/LangBestFil5.table > Filtered_Felix_UG_test/LangBestFil5_DP_hQD.table
mode=$(./software/R-3.4.3/bin/Rscript -e 'library("modes"); data<-read.table("Filtered_Felix_UG_test/LangBestFil5_DP_hQD.table",header=F); cat(modes(data$V1)[1])')
DPmax=2*$mode
DPmin=0.5*$mode

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix_UG_test/LangBestFil5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filterName "DP" -filter "DP<"$DPmin" || DP>"$DPmax"" -o Filtered_Felix_UG_test/LangBestFil6.vcf 2> Filtered_Felix_UG_test/LangBestFil6.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix_UG_test/LangBestFil6.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -o Filtered_Felix_UG_test/LangBestFil7.vcf 2> Filtered_Felix_UG_test/LangBestFil7.err

#Split into cohorts
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix_UG_test/LangBestFil7.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Best.list -o Filtered_Felix_UG_test/BestFil.vcf 2> Filtered_Felix_UG_test/BestFil.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix_UG_test/LangBestFil7.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Lang.list -o Filtered_Felix_UG_test/LangFil.vcf 2> Filtered_Felix_UG_test/LangFil.err

#Extract tables
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix_UG_test/BestFil.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix_UG_test/Best.table 2> Filtered_Felix_UG_test/Best2.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix_UG_test/LangFil.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix_UG_test/Lang.table 2> Filtered_Felix_UG_test/Lang2.err

