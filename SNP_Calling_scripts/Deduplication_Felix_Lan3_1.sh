#!/bin/bash

#Deduplication
#picard tools 2.7.1
#run from pflaphy-cutolgs

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Deduplication_Lan
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi


java -Xmx50g -jar software/picard.jar  MarkDuplicates I=aligned_Felix/Lan3_1.sort.bam O=dedup_Felix/Lan3_1.dedup.bam M=dedup_Felix/duplicateMetricsFile_Lan3_1 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> dedup_Felix/Lan3_1.err
