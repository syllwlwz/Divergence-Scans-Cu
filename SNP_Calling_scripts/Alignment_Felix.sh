#!/bin/bash

#Alignment
#samtools 1.9 with htslib 1.9
#bwa 0.7.17
#run from pflaphy-gscan

#$ -t 1-56
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Alignment_Felix
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 16

file=$(cat cutadapt_list.txt | head -n $SGE_TASK_ID | tail -n 1)
./software/bwa-0.7.17/bwa mem -t 16 -k 10 Ahalleri_Felix_renamed_sorted.fasta -M cutadapt/$file.R1.cutadapt.fastq.gz cutadapt/$file.R2.cutadapt.fastq.gz > aligned_Felix/$file.sam 2> aligned_Felix/$file.sam.err

#t: number of threads
#k: minimum length of seed region
 
./software/samtools-1.9/samtools view -b aligned_Felix/$file.sam | ./software/samtools-1.9/samtools sort -T $file > aligned_Felix/$file.sort.bam 2> aligned_Felix/$file.sort.bam.err

#b: output in bam format
#T: write temporary files

./software/samtools-1.9/samtools index aligned_Felix/$file.sort.bam

./software/samtools-1.9/samtools idxstats aligned_Felix/$file.sort.bam > aligned_Felix/$file.idxstats

./software/samtools-1.9/samtools flagstat aligned_Felix/$file.sort.bam > aligned_Felix/$file.flagstat

rm aligned_Felix/$file.sam

