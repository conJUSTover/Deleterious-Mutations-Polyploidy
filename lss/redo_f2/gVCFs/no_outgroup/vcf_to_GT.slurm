#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 26   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "GT"   # job name
#SBATCH --mem=250G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-GT2.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=biocrunch,speedy #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel 

parallel 'bash vcf_to_GT.sh {} "2>" vcf_to_GT.{}.err' ::: {01..26} 
sbatch combine_vcfs.slurm 
