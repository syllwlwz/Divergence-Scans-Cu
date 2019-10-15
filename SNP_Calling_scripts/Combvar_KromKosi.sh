#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Variantcombi_KromKosi
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T CombineVariants -nt 4 -V Filtered_Felix/KromKosiFil8.vcf -V Filtered_Felix/KromKosiFil_INDEL8.vcf -genotypeMergeOptions UNSORTED -R Ahalleri_Felix_renamed_sorted.fasta -o Filtered_Felix/KromKosiSNPIndels.vcf 2> Filtered_Felix/KromKosiSNPIndels.err

