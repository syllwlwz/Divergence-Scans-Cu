#!/bin/bash

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Nielsen_KromKosi_folded
#$ -l vf=5G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 12

./../software/R-3.5.3/bin/Rscript NielsenKromKosi_folded_parallel2.R
