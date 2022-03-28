#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 9   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gVCF"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-vcftools.biallele.sitepi.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load vcftools/0.1.14-xp36ajy
module load parallel

#bcftools reheader --samples Australian/header.txt -o pSONIC.reheader.vcf pSONIC.combined.SNPs.recode.vcf
#bcftools annotate -x ^FORMAT/GT -o pSONIC.GT.vcf pSONIC.reheader.vcf
#ls *lst | parallel 'vcftools --vcf all.biallele.vcf --keep {} --freq2 --out {.}.biallele'
#head -100 pSONIC.GT.vcf | grep "#CHROM" | sed 's/#//g' | cut -f 1,2,4,5 > pSONIC.GT.pos.ref.alt.txt
#grep -v "#" pSONIC.GT.vcf | cut -f 1,2,4,5 >> pSONIC.GT.pos.ref.alt.txt
ls frqs/[1-9].lst | parallel vcftools --vcf all.biallele.pSONIC.recode.vcf --keep {} --site-pi --out {.}

