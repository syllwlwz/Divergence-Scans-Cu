#!/bin/bash

#Realignment
#GATK3.7

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N HC1_Felix_Lan
#$ -l vf=15G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 4

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T HaplotypeCaller -I realigned_Felix/Lan3_1.realigned.bam --min_base_quality_score 25 --min_mapping_quality_score 25 -rf DuplicateRead -rf BadMate -rf BadCigar -ERC BP_Resolution -R Ahalleri_Felix_renamed_sorted.fasta -o HC1_Felix/Lan3_1.vcf.gz -ploidy 2 --pcr_indel_model NONE -nct 4 --heterozygosity 0.02  --indel_heterozygosity 0.01 -rf NotPrimaryAlignment --output_mode EMIT_ALL_SITES 2> HC1_Felix/Lan3_1.err

