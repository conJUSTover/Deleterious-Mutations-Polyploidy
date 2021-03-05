#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 5   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "At gerpelem"   # job name
#SBATCH --mem=900G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-gerpelem_At.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load parallel 
module load gcc/6.4.0-vddll7o

cat run_evalues.At.txt | parallel /home/jconover/gerp_recompile/gerpelem -v -f roast.At.msa.in.new.raw.rates -e {} -x ".{}.elems" 

