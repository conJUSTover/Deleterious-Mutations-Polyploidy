#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "090 gVCF"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-007.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

cat 007.1.fq.gz 006.1.fq.gz > 090.1.fq.gz &&
cat 007.2.fq.gz 006.2.fq.gz > 090.2.fq.gz &&
bash Sentieon_Mapping_to_gVCF.sh 090 30 Ghirsutum_458_v1.0.26chr.fa

