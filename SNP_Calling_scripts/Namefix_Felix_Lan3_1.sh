#!/bin/bash

#Namefix
#picard 2.7.1

#run from pflaphy-gscan

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Namefix_Felix_lan
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -Xmx50g -jar software/picard.jar AddOrReplaceReadGroups I=dedup_Felix/Lan3_1.dedup.bam O=namefixed_Felix/Lan3_1.dedup.bam SORT_ORDER=coordinate CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT RGLB=Lan3_1 RGPL=illumina RGPU=AAAAAA RGSM=Lan3_1 2> namefixed_Felix/Lan3_1.err
