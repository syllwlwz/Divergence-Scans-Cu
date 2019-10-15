#!/bin/bash

#GATK3.7
#cohorts separately

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Flk_comb
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 10

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/KromKosinewFil8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --concordance Filtered_Felix/KosinewFil3.vcf -o Filtered_Felix/KromKosinew.vcf 2> Filtered_Felix/KromKosinew.err
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/NossPaisnewFil8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --concordance Filtered_Felix/NossnewFil3.vcf -o Filtered_Felix/NossPaisnew.vcf 2> Filtered_Felix/NossPaisnew.err
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/LangBestnewFil8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --concordance Filtered_Felix/LangnewFil3.vcf -o Filtered_Felix/LangBestnew.vcf 2> Filtered_Felix/LangBestnew.err
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T SelectVariants -V Filtered_Felix/WulmHDnewFil8.vcf -R Ahalleri_Felix_renamed_sorted.fasta -nt 10 --concordance Filtered_Felix/WulmnewFil3.vcf -o Filtered_Felix/WulmHDnew.vcf 2> Filtered_Felix/WulmHDnew.err

