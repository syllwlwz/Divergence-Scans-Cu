#!/bin/bash

#GATK version 3.7
#cohorts separately

#$ -t 1-2
#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N SweeD_KromKosi_allSNPs
#$ -l vf=70G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

file=$(cat ../GROM.list | head -n $SGE_TASK_ID | tail -n 1)
./../software/sweed-master/SweeD -name $file'_SweeD_allSNPs' -input Sweed_$file.vcf -grid 37500 -folded -noSeparator -reports
