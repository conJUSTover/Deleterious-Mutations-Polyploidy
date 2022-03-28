#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "1.0 1.0.30"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools.1.0.30-1.0.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load vcftools/0.1.14-xp36ajy

out="1.0"
middle="30.0"
vcf="vcf"

#for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --freq2 --indv "100 " --indv "101 " --indv "102 " --indv "103 " --indv "104 " --indv "105 " --indv "106 " --indv "107 " --indv "108 " --indv "109 " ; done

#for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --window-pi 50000 --indv "100 " --indv "101 " --indv "102 " --indv "103 " --indv "104 " --indv "105 " --indv "106 " --indv "107 " --indv "108 " --indv "109 " ; done

for i in {22..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --site-pi --indv "100 " --indv "101 " --indv "102 " --indv "103 " --indv "104 " --indv "105 " --indv "106 " --indv "107 " --indv "108 " --indv "109 " ; done

#--indv "100 " --indv "101 " --indv "102 " --indv "103 " --indv "104 " --indv "105 " --indv "106 " --indv "107 " --indv "108 " --indv "109 " ; done
