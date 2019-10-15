#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Vartotable_SNPindel_KromKosi
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

#Split into cohorts
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosiSNPIndelsnew.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Kosinew.list -o Filtered_Felix/KosiSNPIndels.vcf 2> Filtered_Felix/KosiSNPIndels.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosiSNPIndelsnew.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 -sf Krom.list -o Filtered_Felix/KromSNPIndels.vcf 2> Filtered_Felix/KromSNPIndels.err

#Extract tables
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/KosiSNPIndels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/KosiSNPIndels.table 2> Filtered_Felix/KosiSNPIndels2.err

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V Filtered_Felix/KromSNPIndels.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F AC -F AN -raw --allowMissingData -o Filtered_Felix/KromSNPIndels.table 2> Filtered_Felix/KromSNPIndels2.err

