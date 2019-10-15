#!/bin/bash

#run from pflaphy-gscan


#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N SNPeff_KromKosi
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#java -Xmx50g -jar software/snpEff/snpEff.jar AHAL_Felix Filtered_Felix/KromKosiFil8.vcf -t > SNPeff_scans/KromKosi.ann.vcf 2> SNPeff_scans/KromKosi.ann.vcf.err

#Extract tables
java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -V SNPeff_scans/KromKosi.ann.vcf -R Ahalleri_Felix_renamed_sorted.fasta -F CHROM -F POS -F ANN -raw --allowMissingData -o SNPeff_scans/KromKosi.ann.table 2> SNPeff_scans/KromKosi.ann.table.err
