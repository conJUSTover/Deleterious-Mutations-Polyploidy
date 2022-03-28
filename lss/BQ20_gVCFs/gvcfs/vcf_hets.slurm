#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 2:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 26   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "hets"   # job name
#SBATCH --mem=12G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-hets.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load python/3.8.0-k7w5uj4
module load parallel 

parallel 'python vcf_all_hets.py {}.vcf > {}.hets' ::: {01..26}

cat *hets | gzip > all.hets.gz
