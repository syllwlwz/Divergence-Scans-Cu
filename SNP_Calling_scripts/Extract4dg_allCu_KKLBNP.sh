#!/bin/bash

#GATK version 3.7
#cohorts separately

#$ -cwd
#$ -v PATH=/vol/biotools/bin:/usr/bin
#$ -N Extract4dg_KKLBNP
#$ -l vf=25G
#$ -l arch=lx-amd64
#$ -l idle=1
#$ -P denbi

java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP.vcf -L Ahalleri_Felix_4dg_sites.intervals.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf
#java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf -sf Krom.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Krom.vcf
#java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf -sf Kosi.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Kosi.vcf
#java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf -sf Lang.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Lang.vcf
#java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf -sf Best.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Best.vcf
#java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf -sf Noss.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Noss.vcf
#java -jar software/GenomeAnalysisTK37.jar -T SelectVariants -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites.vcf -sf Pais.list -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Pais.vcf

#java -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites_Krom.vcf -F CHROM -F POS -F AC -F AF -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Krom.table
#java -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites_Kosi.vcf -F CHROM -F POS -F AC -F AF -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Kosi.table
#java -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites_Lang.vcf -F CHROM -F POS -F AC -F AF -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Lang.table
#java -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites_Best.vcf -F CHROM -F POS -F AC -F AF -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Best.table
#java -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites_Noss.vcf -F CHROM -F POS -F AC -F AF -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Noss.table
#java -jar software/GenomeAnalysisTK37.jar -T VariantsToTable -nt 4 -R Ahalleri_Felix_renamed_sorted.fasta -V Filtered_Felix/AllCucombKKLBNP_4dg_sites_Pais.vcf -F CHROM -F POS -F AC -F AF -o Filtered_Felix/AllCucombKKLBNP_4dg_sites_Pais.table

