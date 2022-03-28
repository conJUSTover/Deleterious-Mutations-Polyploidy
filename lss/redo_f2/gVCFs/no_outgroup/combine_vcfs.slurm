#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "combine"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-combine.vcfs.pSONIC.GT.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=biocrunch,speedy #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load bcftools/1.9-womp5fh
module load vcftools/0.1.14-xp36ajy

bcftools concat -o all.GT.vcf *GT_sorted.vcf
vcftools --vcf all.GT.vcf --bed pSONIC.percentfiltered.375_625.sorted.bed --recode --out pSONIC.GT
#head -1000 pSONIC.GT.recode.vcf | grep "#CHROM" | cut -f 1,2 > pSONIC.seg.pos
grep -v "#" pSONIC.GT.recode.vcf | cut -f 1,2 >> pSONIC.seg.pos

#vcftools --vcf pSONIC.biallele.GT.vcf --keep D5.lst --freq2 --out D5.biallele.pSONIC
#vcftools --vcf pSONIC.biallele.GT.vcf --keep A1.lst --freq2 --out A1.biallele.pSONIC
#vcftools --vcf pSONIC.biallele.GT.vcf --keep AD.lst --freq2 --out AD.biallele.pSONIC

#vcftools all.GT.vcf --keep D5.lst --freq2 --out D5.GT
#vcftools all.GT.vcf --keep A1.lst --freq2 --out A1.GT
#vcftools all.GT.vcf --keep AD.lst --freq2 --out AD.GT

