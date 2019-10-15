#!/bin/bash

#GATK3.7
#cohorts separately

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Faststructure_plot
#$ -l vf=25G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

python2 software/fastSTRUCTURE/fastStructure-master/distruct.py -K 5 --input=Faststructure/fS_run_K --output=Faststructure/Faststructure_plot_distruct.svg 2> Faststructure/Faststructure_plot.err

