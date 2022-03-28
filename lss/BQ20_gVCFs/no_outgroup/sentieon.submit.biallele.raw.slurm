#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "CHROM bial"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-CHROM_biallele.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load bcftools/1.9-womp5fh
module load vcftools/0.1.14-xp36ajy
module load vcflib

vcftools --vcf CHROM.vcf --remove-filtered-geno-all --remove-indels --min-alleles 2 --max-alleles 2 --recode --out CHROM.biallele;
vcftools --vcf CHROM.biallele.recode.vcf --depth;
vcffilter -f "MQ > 30 & GQ > 20" CHROM.biallele.recode.vcf > CHROM.HighQuality.vcf
