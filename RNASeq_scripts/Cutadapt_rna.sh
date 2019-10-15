#!/bin/bash

#cutadapt version 1.18

#$ -t 1-24
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N cutadapt_rna
#$ -l vf=2G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

#ls Fastq/*.fq.gz | sed 's/.fq.gz//g' | sed 's/Fastq\///g' | sed 's/\.1\|\.2//g' | uniq > Cutadapt.list
file=$(cat Cutadapt.list | head -n $SGE_TASK_ID | tail -n 1)
./software/cutadapt-1.18/cutadapt -e 0.15 -O 4 -m 120 --nextseq-trim=20 -q 20 -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' -a "A{100}" -n 2 -o cutadapt/$file.1.cutadapt.fastq.gz -p cutadapt/$file.2.cutadapt.fastq.gz Fastq/$file.1.fq.gz Fastq/$file.2.fq.gz > cutadapt/$file.cutadapt 2> cutadapt/$file.cutadapt.err

#chmod 755 cutadapt.sh
#nohup ./cutadapt.sh &
#option -j number of cores for multithreading available together with pigz installed! 

