#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 49   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Asort 090.sort.bam"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-090.sort.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=biocrunch,bigram,speedy

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load samtools/1.9-k6deoga

samtools view -f 2 -b -@ 49 ../090.At.sort.bam -o 090.sort.bam && samtools index -@ 49 090.sort.bam && echo "090.sort.bam completed" && sbatch gVCF.090.slurm
#&& for i in sentieon.submit.mapping.30*; do sbatch $i; done 
#samtools index -@ $2 $1.sort.bam
