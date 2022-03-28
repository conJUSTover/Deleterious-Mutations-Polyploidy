#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 11:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "site-pi"   # job name
#SBATCH --mem=50G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=site-pi.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever --exclude=legion-[1-8] #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load vcftools/0.1.14-xp36ajy

#vcftools --gzvcf ../all.GT.vcf.gz --weir-fst-pop ../AD1D.lst --weir-fst-pop ../AD2D.lst --fst-window-size 500000 --fst-window-step 50000 --out AD1D.AD2D.fst
ls *lst | parallel -j 30 vcftools --gzvcf ../all.GT.vcf.gz --site-pi --keep {} --out frqs/{.}

