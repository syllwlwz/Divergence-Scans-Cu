#!/bin/bash


#$ -t 1-3
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Compcounts_Felix
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Cutadapt.list | head -n $SGE_TASK_ID | tail -n 1)
./software/qualimap_v2.2.1/qualimap comp-counts -algorithm proportional -bam mapped_Felix/$file'_corrected.bam' -gtf AHAL_sorted2.gtf -id gene_id -out Counts/$file'_counts.txt' -p non-strand-specific -pe -type exon 2> Compcounts_Felix_$file.err
