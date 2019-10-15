#!/bin/bash


#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Identify4dgsites
#$ -l vf=4G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

perl software/Identify_4D_Sites.pl AHAL_sorted.gff3 Ahalleri_Felix_renamed_sorted.fasta
