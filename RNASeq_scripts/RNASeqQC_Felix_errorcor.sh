#!/bin/bash


#$ -t 1-24
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N RNASeqqc_Felix_errorcor
#$ -l vf=7G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt.list | head -n $SGE_TASK_ID | tail -n 1)
./software/qualimap_v2.2.1/qualimap rnaseq -a proportional -bam mapped_Felix/$file'_corrected.bam' -gtf AHAL_sorted2.gtf -pe 2> mapped_Felix/$file.errorcor.rnaseqqc.err

