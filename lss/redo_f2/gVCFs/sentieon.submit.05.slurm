#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 10   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Chr 05"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-05.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load sentieon-genomics/201808.01-opfuvzr


time ( sentieon driver -t 10 -r Ghirsutum_458_v1.0.26chr.fa --interval 05 --algo GVCFtyper 05.vcf [0-9]*.gVCF) &> 05.vcf.timelog
