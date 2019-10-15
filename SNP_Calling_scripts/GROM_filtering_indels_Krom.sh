#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Filtering_Krom_GROM
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#Filter 1: select bialleleic SNPs without missing genotypes
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V:VCF GROM_analysis/Krom.indels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -sf Krom.list -selectType INDEL --excludeNonVariants -restrictAllelesTo BIALLELIC -env -o GROM_analysis/KromGROMFil_INDEL1.vcf 2> GROM_analysis/KromGROMFil_INDEL1.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V GROM_analysis/KromGROMFil_INDEL1.vcf -R Ahalleri_Felix_renamed_sorted.fasta -G_filter "DP<2.0" -G_filterName "lowDP" -o GROM_analysis/KromGROMFil_INDEL2.vcf 2> GROM_analysis/KromGROMFil_INDEL2.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V GROM_analysis/KromGROMFil_INDEL2.vcf -R Ahalleri_Felix_renamed_sorted.fasta --setFilteredGtToNocall -o GROM_analysis/KromGROMFil_INDEL3.vcf 2> GROM_analysis/KromGROMFil_INDEL3.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 12 -V GROM_analysis/KromGROMFil_INDEL3.vcf -R Ahalleri_Felix_renamed_sorted.fasta -select "AN >= 20" -o GROM_analysis/KromGROMFil_INDEL4.vcf 2> GROM_analysis/KromGROMFil_INDEL4.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V GROM_analysis/KromGROMFil_INDEL4.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filter "QD<20.0 || FS>40.0 || ReadPosRankSum < -4.0 || SOR>4.0" -filterName "BP" -o GROM_analysis/KromGROMFil_INDEL5a.vcf 2> GROM_analysis/KromGROMFil_INDEL5a.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V GROM_analysis/KromGROMFil_INDEL5a.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered  -o GROM_analysis/KromGROMFil_INDEL5.vcf 2> GROM_analysis/KromGROMFil_INDEL5.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V GROM_analysis/KromGROMFil_INDEL5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F DP -raw --allowMissingData -o GROM_analysis/KromGROMFil_INDEL5.table 2> GROM_analysis/KromGROMFil_INDEL5_DPtable.err

grep -v "DP" GROM_analysis/KromGROMFil_INDEL5.table > GROM_analysis/KromGROMFil_INDEL5_DP_hQD.table
mode=$(./software/R-3.4.3/bin/Rscript -e 'library("modes"); data<-read.table("GROM_analysis/KromGROMFil_INDEL5_DP_hQD.table",header=F); cat(modes(data$V1)[1])')
DPmax=2*$mode
DPmin=0.5*$mode

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantFiltration -V GROM_analysis/KromGROMFil_INDEL5.vcf -R Ahalleri_Felix_renamed_sorted.fasta -filterName "DP" -filter "DP<"$DPmin" || DP>"$DPmax"" -o GROM_analysis/KromGROMFil_INDEL6.vcf 2> GROM_analysis/KromGROMFil_INDEL6.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V GROM_analysis/KromGROMFil_INDEL6.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --excludeFiltered -o GROM_analysis/KromGROMFil_INDEL7.vcf 2> GROM_analysis/KromGROMFil_INDEL7.err

#Extract tables
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V GROM_analysis/KromGROMFil_INDEL7.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o GROM_analysis/KromGROM_INDEL.table 2> GROM_analysis/KromGROM_INDEL2.err

