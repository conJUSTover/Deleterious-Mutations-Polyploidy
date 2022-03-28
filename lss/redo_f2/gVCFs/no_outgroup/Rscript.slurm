#!/bin/bash

#Submit this script with: sbatch thefilename
# note: change the memory, threads, wall, etc

#SBATCH -t 24:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores
#SBATCH -J "Rscript"   # job name
#SBATCH --mem=340G # how much memory you need; each box has ~350G
#SBATCH --output=slurm-Rscript_Tcacao.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load gcc/10.2.0-zuvaafu
module load r/4.0.4-py3-4khjixy
module load r-stringr/1.4.0-py3-wnm4rst
module load r-dplyr/1.0.3-py3-4pcfkgs
module load r-reshape2/1.4.4-py3-x5mbpmn

Rscript add_outgroup_to_POS.R 
