#!/bin/bash

#GATK version 3.7
#cohorts separately


#$ -t 1-6
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Average_coverage_othersamples
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat Population.list | head -n $SGE_TASK_ID | tail -n 1)

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T DepthOfCoverage -I $file.list -R Ahalleri_Felix_renamed_sorted.fasta -o $file.mean_DOC_all -mbq 25 -mmq 25 -omitBaseOutput -omitIntervals --omitLocusTable 2> $file.mean_DOC.err

