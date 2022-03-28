#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "801 gVCF"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-801.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=legion,speedy,biocrunch

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

bash Sentieon_Mapping_to_gVCF.sh 801 30 Ghirsutum_458_v1.0.26chr.fa

