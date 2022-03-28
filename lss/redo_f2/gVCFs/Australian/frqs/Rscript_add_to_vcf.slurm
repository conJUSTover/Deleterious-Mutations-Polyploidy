#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 9   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "polarize"   # job name
#SBATCH --mem=200G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=merge_frqs.pSONIC.Other_GC.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load r/3.6.3-py3-sxv6dw3
module load r-dplyr/0.7.5-py2-r3.5-4jcgekn
module load r-stringr/1.2.0-py2-r3.5-uoogz7r 
module load r-reshape2/1.4.2-py2-ywy5u5c

#Rscript add_GERP_BAD_vcf.R && echo "add_GERP_BAD_vcf.R completed successfully"
#Rscript vcf_to_models.new.R && echo "vcf_to_models.R completed successfully"
#Rscript polarize_vcf.R && echo "polarization completed successfully"

#sbatch filter_SNPs.slurm

module load vcftools/0.1.14-xp36ajy
module load parallel 

#ls [0-9].lst out.lst | parallel vcftools --gzvcf ../all.GT.vcf.gz --keep {} --freq2 --out {.}

#Rscript merge_frqs.R

#ls [0-9].lst out.lst | parallel vcftools --vcf pSONIC.NO_GC.recode.vcf --keep {} --freq2 --out {.}.pSONIC.noGC
ls [0-9].lst out.lst | parallel vcftools --vcf pSONIC.At_GC.recode.vcf --keep {} --freq2 --out {.}.pSONIC.At_GC
ls [0-9].lst out.lst | parallel vcftools --vcf pSONIC.Dt_GC.recode.vcf --keep {} --freq2 --out {.}.pSONIC.Dt_GC
ls [0-9].lst out.lst | parallel vcftools --vcf pSONIC.Amb_GC.recode.vcf --keep {} --freq2 --out {.}.pSONIC.Amb_GC


ls merge_frqs_pSONIC.[AD]*R | parallel Rscript {}
