#!/bin/bash

#GROM

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N GROM_Krom
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

./software/GROM/dist/GROM -i GROM_analysis/Krom.bam -r Ahalleri_Felix_renamed_sorted.fasta -o GROM_analysis/Krom.predictions.vcf -P 6 -b 25 -q 25 -v 0.02 -e 0.01 -V 0.02 -U 3 -A 4 -p 2 -s 30 2> GROM_analysis/Krom.GROM.err

#Select INDELs and CNVs

#grep "SPR:SEV:SRD:SCO:ECO:SOT:EOT:SSC:HP\|SPR:EPR:SEV:EEV:SRD:ERD:SCO:ECO:SOT:EOT:SSC:ESC:HP"  GROM_analysis/$file.predictions.vcf | grep -v "#" > GROM_analysis/$file.indels.vcf
#grep "SD:Z:CN:CS"  GROM_analysis/$file.predictions.vcf | grep -v "#" > GROM_analysis/$file.cnvs.vcf

