#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "vcf merge"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-merge_biallele.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever # -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load bcftools/1.9-womp5fh
#module load vcftools/0.1.14-xp36ajy
module load sentieon-genomics/201808.01-opfuvzr


#sentieon util vcfindex 912.gVCF
for i in C2-4.gVCF G2-21.gVCF G3-83.gVCF; do bcftools concat -o $i ../../D_genome/$i ../../A_genome/$i && sentieon util vcfindex $i; done && for i in sentieon.submit.[0-9][0-9].slurm; do sbatch $i; done 
#bcftools concat *vcf | vcftools --vcf - --remove-indels --min-alleles 2 --max-alleles 2 --recode --out all.chr.biallele

