#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 6:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 32   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "copy"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-copy2.out
#SBATCH --mail-type=FAIL
#SBATCH --constraint=AVX2
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load parallel 

cat done_firstpass.lst | sed 's/_Predictions.txt//g' | parallel -j 32 'mv Align_logs/{}* Align_logs/done/; mv Predictions_Dir/{}* Predictions_Dir/done/; mv Predict_logs/{}* Predict_logs/done/; mv Output_Dir/{}* Output_Dir/done/; mv subs/{}* subs/done/'
