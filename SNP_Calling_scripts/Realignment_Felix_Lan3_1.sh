#!/bin/bash

#Realignment
#GATK3.7

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Realignment_Felix_Lan
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T RealignerTargetCreator -R Ahalleri_Felix_renamed_sorted.fasta -o namefixed_Felix/Lan3_1.IndelRealigner.intervals -I namefixed_Felix/Lan3_1.dedup.bam 2> realigned_Felix/Lan3_1.Interval.err; java -Xmx50g -jar software/GenomeAnalysisTK37.jar -T IndelRealigner -targetIntervals namefixed_Felix/Lan3_1.IndelRealigner.intervals -I namefixed_Felix/Lan3_1.dedup.bam -R Ahalleri_Felix_renamed_sorted.fasta -o realigned_Felix/Lan3_1.realigned.bam 2> realigned_Felix/Lan3_1.err

