#!/bin/bash

#GATK version 3.7
#cohorts separately

#run from pflaphy-gscan


#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N UG_Felix_LangBest
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T UnifiedGenotyper -I LangBest_realigned.list --min_base_quality_score 25 -rf DuplicateRead -rf BadMate -rf BadCigar -R Ahalleri_Felix_renamed_sorted.fasta -o UG_Felix/LangBest_UG.vcf.gz -ploidy 2 -nct 12 -rf NotPrimaryAlignment -glm SNP --heterozygosity 0.02 -stand_call_conf 25 --output_mode EMIT_VARIANTS_ONLY -dcov 200 2> UG_Felix/LangBest_UG.err



#read filter NotPrimaryAlignment filters out multimapping reads
