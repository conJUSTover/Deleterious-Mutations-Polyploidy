#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "1.3 1.3.40"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools.1.3.40-1.3.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load vcftools/0.1.14-xp36ajy

out="1.3"
middle="40.0"
vcf="vcf"

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --freq2 --indv "130 " --indv "131 " --indv "132 " --indv "133 " --indv "134 " --indv "135 " --indv "136 " --indv "137 " --indv "138 " --indv "139 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --window-pi 50000 --indv "130 " --indv "131 " --indv "132 " --indv "133 " --indv "134 " --indv "135 " --indv "136 " --indv "137 " --indv "138 " --indv "139 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --site-pi --indv "130 " --indv "131 " --indv "132 " --indv "133 " --indv "134 " --indv "135 " --indv "136 " --indv "137 " --indv "138 " --indv "139 " ; done

#--indv "130 " --indv "131 " --indv "132 " --indv "133 " --indv "134 " --indv "135 " --indv "136 " --indv "137 " --indv "138 " --indv "139 " ; done
