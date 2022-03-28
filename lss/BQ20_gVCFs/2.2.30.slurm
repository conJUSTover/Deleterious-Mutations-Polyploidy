#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "2.2 2.2.30"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools.2.2.30-2.2.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load vcftools/0.1.14-xp36ajy

out="2.2"
middle="30.0"
vcf="vcf"

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --freq2 --indv "220 " --indv "221 " --indv "222 " --indv "223 " --indv "224 " --indv "225 " --indv "226 " --indv "227 " --indv "228 " --indv "229 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --window-pi 50000 --indv "220 " --indv "221 " --indv "222 " --indv "223 " --indv "224 " --indv "225 " --indv "226 " --indv "227 " --indv "228 " --indv "229 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --site-pi --indv "220 " --indv "221 " --indv "222 " --indv "223 " --indv "224 " --indv "225 " --indv "226 " --indv "227 " --indv "228 " --indv "229 " ; done

#--indv "220 " --indv "221 " --indv "222 " --indv "223 " --indv "224 " --indv "225 " --indv "226 " --indv "227 " --indv "228 " --indv "229 " ; done
