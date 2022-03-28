#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gzip ave"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=awk_average.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel 
module load bedtools2/2.27.1-s2mtpsu


#ls *.gVCF | parallel -j 32 'bedtools coverage -a {} -b Ghir.chr.bed -hist > {.}.hist'
#ls *hist | parallel -j 32 "awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' {} > {.}.ave"

#for i in *gVCF; do awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' ${i%.gVCF}.hist > $i.average; done 

ls *hist | parallel -j 30 'gzip {} > {}.gzip'

#for i in *hist; do gzip $i > $i.gzip && echo $i; done 
