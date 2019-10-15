#!/bin/bash

#samtools 1.9.1 with htslib 1.9.1

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N samtools_merge
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

./software/samtools-1.9/samtools merge -b Kosi_realigned.list GROM_analysis/Kosi.bam 2> GROM_analysis/Kosi.err

./software/samtools-1.9/samtools merge -b Krom_realigned.list GROM_analysis/Krom.bam 2> GROM_analysis/Krom.err

./software/samtools-1.9/samtools merge -b Noss_realigned.list GROM_analysis/Noss.bam 2> GROM_analysis/Noss.err

./software/samtools-1.9/samtools merge -b Pais_realigned.list GROM_analysis/Pais.bam 2> GROM_analysis/Pais.err

./software/samtools-1.9/samtools merge -b Wulm_realigned.list GROM_analysis/Wulm.bam 2> GROM_analysis/Wulm.err

./software/samtools-1.9/samtools merge -b HD_realigned.list GROM_analysis/HD.bam 2> GROM_analysis/HD.err

./software/samtools-1.9/samtools merge -b Lang_realigned.list GROM_analysis/Lang.bam 2> GROM_analysis/Lang.err

./software/samtools-1.9/samtools merge -b Best_realigned.list GROM_analysis/Best.bam 2> GROM_analysis/Best.err

