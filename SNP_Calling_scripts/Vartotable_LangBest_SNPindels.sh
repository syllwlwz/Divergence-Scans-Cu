#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Vartotable_SNPindel_LangBest
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#Split into cohorts
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/LangBestSNPIndels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Best.list -o Filtered_Felix/BestSNPIndels.vcf 2> Filtered_Felix/BestSNPIndels.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/LangBestSNPIndels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Lang.list -o Filtered_Felix/LangSNPIndels.vcf 2> Filtered_Felix/LangSNPIndels.err

#Extract tables
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/BestSNPIndels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/BestSNPIndels.table 2> Filtered_Felix/BestSNPIndels2.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/LangSNPIndels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/LangSNPIndels.table 2> Filtered_Felix/LangSNPIndels2.err

