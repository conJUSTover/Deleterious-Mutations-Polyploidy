#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "9 9.40"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools.9.40-9.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load vcftools/0.1.14-xp36ajy

out="9"
middle="40.0"
vcf="vcf"

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --freq2 --indv "900 " --indv "901 " --indv "902 " --indv "903 " --indv "904 " --indv "910 " --indv "911 " --indv "912 " --indv "913 " --indv "914 " --indv "915 " --indv "916 " --indv "917 " --indv "918 " --indv "919 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --window-pi 50000 --indv "900 " --indv "901 " --indv "902 " --indv "903 " --indv "904 " --indv "910 " --indv "911 " --indv "912 " --indv "913 " --indv "914 " --indv "915 " --indv "916 " --indv "917 " --indv "918 " --indv "919 " ; done

for i in {01..26}; do vcftools --vcf $i.biallele.$middle.recode.vcf --min-alleles 2 --max-alleles 2 --out $i.$out.$middle --site-pi --indv "900 " --indv "901 " --indv "902 " --indv "903 " --indv "904 " --indv "910 " --indv "911 " --indv "912 " --indv "913 " --indv "914 " --indv "915 " --indv "916 " --indv "917 " --indv "918 " --indv "919 " ; done

#--indv "900 " --indv "901 " --indv "902 " --indv "903 " --indv "904 " --indv "910 " --indv "911 " --indv "912 " --indv "913 " --indv "914 " --indv "915 " --indv "916 " --indv "917 " --indv "918 " --indv "919 " ; done
