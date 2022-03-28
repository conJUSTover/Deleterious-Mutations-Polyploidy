#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 23:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 50   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "AG3-83 sort"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-G3-83.sort.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=bigram,biocrunch,speedy

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

bash Sort.sh G3-83 30 Ghirsutum_458_v1.0.26chr.fa && sbatch gVCF.G3-83.slurm

