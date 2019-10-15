#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_KromKosi_indel_hQD_new
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#Filter 1: select bialleleic SNPs without missing genotypes
#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V HC2_Felix/All_Cu_Kosinew.vcf -R Ahalleri_Felix_renamed_sorted.fasta -sf KromKosinew.list -selectType INDEL --excludeNonVariants -restrictAllelesTo BIALLELIC -env -o Filtered_Felix/KromKosinewFil_INDEL1.vcf 2> Filtered_Felix/KromKosinewFil_INDEL1.err

#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix/KromKosinewFil_INDEL1.vcf -R Ahalleri_Felix_renamed_sorted.fasta -G_filter "DP<2.0" -G_filterName "lowDP" -o Filtered_Felix/KromKosinewFil_INDEL2.vcf 2> Filtered_Felix/KromKosinewFil_INDEL2.err

#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil_INDEL2.vcf -R Ahalleri_Felix_renamed_sorted.fasta --setFilteredGtToNocall -o Filtered_Felix/KromKosinewFil_INDEL3.vcf 2> Filtered_Felix/KromKosinewFil_INDEL3.err

#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V Filtered_Felix/KromKosinewFil_INDEL3.vcf -R Ahalleri_Felix_renamed_sorted.fasta -select "AN >= 40" -o Filtered_Felix/KromKosinewFil_INDEL4.vcf 2> Filtered_Felix/KromKosinewFil_INDEL4.err

#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix/KromKosinewFil_INDEL4.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filter "QD<20.0 || FS>40.0 || ReadPosRankSum < -4.0 || SOR>4.0" -filterName "BP" -o Filtered_Felix/KromKosinewFil_INDEL5a.vcf 2> Filtered_Felix/KromKosinewFil_INDEL5a.err

#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil_INDEL5a.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered  -o Filtered_Felix/KromKosinewFil_INDEL5.vcf 2> Filtered_Felix/KromKosinewFil_INDEL5.err

#java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/KromKosinewFil_INDEL5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F DP -raw --allowMissingData -o Filtered_Felix/KromKosinewFil_INDEL5.table 2> Filtered_Felix/KromKosinewFil_INDEL5_DPtable.err

#grep -v "DP" Filtered_Felix/KromKosinewFil_INDEL5.table > Filtered_Felix/KromKosinewFil_INDEL5_DP_hQD.table
mode=$(./software/R-3.5.3/bin/Rscript -e 'library("modes"); data<-read.table("Filtered_Felix/KromKosinewFil_INDEL5_DP_hQD.table",header=F); cat(modes(data$V1)[1])')
DPmax=2*$mode
DPmin=0.5*$mode

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V Filtered_Felix/KromKosinewFil_INDEL5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filterName "DP" -filter "DP<"$DPmin" || DP>"$DPmax"" -o Filtered_Felix/KromKosinewFil_INDEL6.vcf 2> Filtered_Felix/KromKosinewFil_INDEL6.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil_INDEL6.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -o Filtered_Felix/KromKosinewFil_INDEL7.vcf 2> Filtered_Felix/KromKosinewFil_INDEL7.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil_INDEL7.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -select "AF<1.0 && AF>0.0" -o Filtered_Felix/KromKosinewFil_INDEL8.vcf 2> Filtered_Felix/KromKosinewFil_INDEL8.err

#Split into cohorts
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil_INDEL8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Kosinew.list -o Filtered_Felix/KosinewFil_INDEL.vcf 2> Filtered_Felix/KosinewFil_INDEL.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil_INDEL8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Krom.list -o Filtered_Felix/KromnewFil_INDEL.vcf 2> Filtered_Felix/KromnewFil_INDEL.err

#Extract tables
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/KosinewFil_INDEL.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/Kosi_INDEL_new.table 2> Filtered_Felix/Kosi_INDEL2new.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/KromnewFil_INDEL.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/Krom_INDEL_new.table 2> Filtered_Felix/Krom_INDEL2new.err

