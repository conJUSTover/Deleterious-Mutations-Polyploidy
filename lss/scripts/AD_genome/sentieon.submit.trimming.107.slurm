#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 2   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "107 trim"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-trim-107.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

bash Sentieon_Trimmomatic.sh 107 2 Ghirsutum_458_v1.0.26chr.fa

