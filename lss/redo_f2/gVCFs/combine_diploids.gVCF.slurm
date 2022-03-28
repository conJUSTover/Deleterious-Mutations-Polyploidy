#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 50   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Chr 01"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-typer.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=biocrunch,speedy #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load bcftools/1.9-womp5fh
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
module load parallel 

#cat diploids.lst | parallel 'bcftools concat -o {}.gVCF ../At/{}.gVCF ../Dt/{}.gVCF && gatk IndexFeatureFile -F {}.gVCF'
 
bcftools concat -o 090.gVCF ../At/090.gVCF ../Dt/090.gVCF && gatk IndexFeatureFile -F 090.gVCF

