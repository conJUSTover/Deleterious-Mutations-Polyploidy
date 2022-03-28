#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "sbatch raw"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-sbatch_raw.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


module load bwa/0.7.17-zhcbtza

# generate this index before loading Sentieon if you ever want the index to finish
bwa index Ghirsutum_458_v1.0.26chr.fa && for i in sentieon.submit.mapping.[0-9]*; do sbatch $i; done 

#ls *.1.fq.gz | sed 's/.1.fq.gz//g' | while read line; do sed "s/#SAMPLE/$line/g" sentieon.submit.mapping.raw.slurm > sentieon.submit.mapping.$line.slurm; done
