#!/bin/bash

#Deduplication
#picard tools 2.7.1
#run from pflaphy-cutolgs

#$ -t 1-102
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Deduplication
#$ -l vf=55G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 2

file=$(cat Sample.list | head -n $SGE_TASK_ID | tail -n 1)

java -XX:ParallelGCThreads=2 -Xmx50g -jar software/picard.jar  MarkDuplicates I=aligned_Felix/$file.sort.bam O=dedup_Felix/$file.dedup.bam M=dedup_Felix/duplicateMetricsFile_$file VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> dedup_Felix/$file.err
