#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "1.2 1.2.50"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools.1.2.50-1.2.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load vcftools/0.1.14-xp36ajy

out="1.2"
middle="50.0"
vcf="vcf"

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --freq2 --indv "120 " --indv "121 " --indv "122 " --indv "123 " --indv "124 " --indv "125 " --indv "126 " --indv "127 " --indv "128 " --indv "129 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --window-pi 50000 --indv "120 " --indv "121 " --indv "122 " --indv "123 " --indv "124 " --indv "125 " --indv "126 " --indv "127 " --indv "128 " --indv "129 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --site-pi --indv "120 " --indv "121 " --indv "122 " --indv "123 " --indv "124 " --indv "125 " --indv "126 " --indv "127 " --indv "128 " --indv "129 " ; done

#--indv "120 " --indv "121 " --indv "122 " --indv "123 " --indv "124 " --indv "125 " --indv "126 " --indv "127 " --indv "128 " --indv "129 " ; done
