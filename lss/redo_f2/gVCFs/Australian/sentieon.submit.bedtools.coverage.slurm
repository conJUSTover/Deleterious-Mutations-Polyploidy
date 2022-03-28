#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 7   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Chr 01"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-bedtools.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=biocrunch,speedy #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load gcc/7.3.0-xegsmw4
module load bedtools2/2.27.1-opcm3ia
module load parallel 

ls *sort.bam | parallel 'bedtools coverage -sorted -d -a ../../singleton.gff -b {} > F260.{.}.coverage.txt'


