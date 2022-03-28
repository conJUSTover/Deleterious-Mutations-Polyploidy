#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "deleasysfs"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=easysfs.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE> 
module load py-pip/9.0.1-py3-dpds55c

python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i all.GT.vcf.gz -p sfs.lst -a -o all.easysfs --proj=20,14,10,8,12,4,8,10,6
