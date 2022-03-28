#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 7-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 30   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "#SAMPLE gVCF"   # job name
#SBATCH --mem=30G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-#SAMPLE.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load bwa/0.7.17-zhcbtza

time ( bwa index AD3.fa )

#module use /opt/rit/local-test/modules
module load sentieon-genomics/201808.01-opfuvzr
module load samtools/1.9-k6deoga
#module load trimmomatic/0.36-lkktrba

samtools faidx AD3.fa
