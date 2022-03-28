#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "samtools depth"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-samtools_depth.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load samtools
#module load parallel 

for i in 7*sort.bam; do echo $i && samtools depth -a $i | awk '{sum+=$3} END { print "Average = ",sum/NR}'; done > 7.depths.lst

#ls *sort.bam | parallel -j 10 "echo {} && samtools -a depth {} | awk '{sum+=$3} END { print \"Average = \",sum/NR}'" > 7.depths.lst
