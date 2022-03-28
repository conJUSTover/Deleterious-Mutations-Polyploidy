#!/bin/bash

#Submit this script with: sbatch thefilename
# note: change the memory, threads, wall, etc

#SBATCH -t 24:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores
#SBATCH -J "awk > 0"   # job name
#SBATCH --mem=340G # how much memory you need; each box has ~350G
#SBATCH --output=slurm-awk_greater_0.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

for i in [0-9]*txt; do awk '{ if($NF > 0) {print $0}}' $i >> all.DerivedCount.GERP06.txt  && echo $i; done && gzip all.DerivedCount.GERP06.txt

