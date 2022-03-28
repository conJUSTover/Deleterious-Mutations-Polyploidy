#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 2-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Count"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=countRD.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=biocrunch,speedy

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE> 
module load r/3.6.3-py3-sxv6dw3
module load r-dplyr/0.7.5-py2-r3.5-4jcgekn
module load r-stringr/1.2.0-py2-r3.5-uoogz7r
module load r-reshape2/1.4.2-py2-ywy5u5c

Rscript count_read_depth.R && echo "counted reads successfully"


