#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 5   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Chr 23"   # job name
#SBATCH --mem=50G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-23.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load sentieon-genomics/201808.01-opfuvzr


time ( sentieon driver -t 5 -r ../Ghirsutum_458_v1.0.26chr.fa --interval 23 --algo GVCFtyper 23.vcf ../[1-9]*.gVCF) &> 23.vcf.timelog
