#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 2   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "SRR3560179 SRA"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-SRR3560179_SRA.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#module load sentieon-genomics/201808.01-opfuvzr
#module load parallel
#module load sra-toolkit

echo "SRR3560179" | while read srr; do wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${srr:0:7}/$srr/${srr}_1.fastq.gz & wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${srr:0:7}/$srr/${srr}_2.fastq.gz; done
