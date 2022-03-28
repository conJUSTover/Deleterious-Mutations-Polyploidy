#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 10   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Chr CHROM"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-CHROM.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load sentieon-genomics/201808.01-opfuvzr


sentieon util vcfindex 07.vcf
sentieon util vcfindex 11.vcf
#sentieon util vcfindex 227.gVCF
#sentieon util vcfindex 103.gVCF

#sbatch sentieon.submit.09.slurm
#sbatch sentieon.submit.10.slurm
#sbatch sentieon.submit.23.slurm
#time ( sentieon driver -t 10 -r Ghirsutum_458_v1.0.26chr.fa --interval CHROM --algo GVCFtyper CHROM.vcf *.gVCF) &> CHROM.vcf.timelog
