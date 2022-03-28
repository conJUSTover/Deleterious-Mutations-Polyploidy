#!/bin/bash

#Submit this script with: sbatch thefilename
# note: change the memory, threads, wall, etc

#SBATCH -t 24:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 26   # total number of processor cores in this job; each node has 272 cores
#SBATCH -J "15 run"   # job name
#SBATCH --mem=340G # how much memory you need; each box has ~350G
#SBATCH --output=slurm-count_15_R.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#module load parallel
module load r
module load r-data-table

#ls *[0-9].R | parallel -j 26 'Rscript {}'
#head -1 01.all.frq.txt > all.frq.txt 
#for i in [0-9]*txt; do grep -v "POS" $i >> all.frq.txt; done 
#awk '{if($3 != "./." || $53 != "-") {print $0}}' all.frq.txt > all.frq.Gkirk.txt
Rscript maf_to_daf.R
