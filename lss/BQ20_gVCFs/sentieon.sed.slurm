#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "sed"   # job name
#SBATCH --mem=3G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-exhets.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#for i in {01..26}; do sed "s/CHROM/$i/g" sentieon.submit.biallele.raw.slurm > sentieon.submit.biallele.$i.slurm; done 

for i in no_excessHet/[012][0-9].vcf; do grep -v "#" $i | cut -f 1,2,8 | sed 's/\t[^\t]*ExcessHet=/\t/g' | sed 's/;.*//g' > ${i#no_excessHet/}.ExcessHets.vcf; done 
