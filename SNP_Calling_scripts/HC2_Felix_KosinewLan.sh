#!/bin/bash

#GATK version 3.7
#cohorts separately

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N HaplotypeCaller_Kosinew_Lan
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T GenotypeGVCFs -R Ahalleri_Felix_renamed_sorted.fasta -nt 12 -o HC2_Felix/All_Cu_Kosinew_Lan.vcf -hets 0.02  -indelHeterozygosity 0.01 -newQual -V HC_KosinewLan.list 2> HC2_Felix/All_CuKosinewLan.err

