#!/bin/bash


#$ -t 1-3
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Hisat2_gemmifera_adjusted_new10
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

file=$(cat Cutadapt.list | head -n $SGE_TASK_ID | tail -n 1)
./software/hisat2-2.1.0/hisat2 -p 12 --no-unal -q -x Halleri.renamed.index -1 cutadapt/$file.1.cutadapt.fastq.gz -2 cutadapt/$file.2.cutadapt.fastq.gz -S mapped_gemmifera/$file.adjusted_new10.out.sam --summary-file mapped_gemmifera/$file.adjusted_new10.summary.txt --max-intronlen 20000 -L 10 --pen-noncansplice 16 --mp 2,0 --score-min L,0,-0.12 --rdg 3,1 --rfg 3,1 --pen-canintronlen G,-0.5,0.1 2> mapped_gemmifera/$file.adjusted_new10.err

