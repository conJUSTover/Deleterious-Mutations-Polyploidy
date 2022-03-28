#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "15 bial"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-15_biallele.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load bcftools/1.9-womp5fh
module load vcftools/0.1.14-xp36ajy
module load vcflib

vcftools --vcf 15.vcf --remove-filtered-geno-all --remove-indels --min-alleles 2 --max-alleles 2 --recode --out 15.biallele;
vcftools --vcf 15.biallele.recode.vcf --depth;
vcffilter -f "MQ > 30 & GQ > 20" 15.biallele.recode.vcf > 15.HighQuality.vcf
