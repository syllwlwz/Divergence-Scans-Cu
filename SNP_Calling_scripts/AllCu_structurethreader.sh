#!/bin/bash

#GATK3.7
#cohorts separately

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Structurethreader
#$ -l vf=25G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi
#$ -pe multislot 8

./software/structure_threader-1.2.14/structure_threader-1_2_14 run -K 15 -R 1 -i Filtered_Felix/AllCucombKKLBNP.bed -t 8 -o Faststructure/ -fs software/fastSTRUCTURE/fastStructure-master/structure.py --ind Ind_allCu.txt --extra_opts prior=logistic --extra_opts full --extra_opts seed=100 2> Faststructure/Structurethreader.err

