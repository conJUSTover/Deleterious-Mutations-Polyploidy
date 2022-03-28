#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 20:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "BAD concat"   # job name
#SBATCH --mem=40G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-redo_BAD_compile.out
#SBATCH --mail-type=FAIL
#SBATCH --constraint=AVX2
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,bigram

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/4.3.30-qdauveb
module load parallel

source activate bad_mutations
export TMPDIR=/work/LAS/jfw-lab/jconover/deleterious_mutations/BAD_runs/redo_Australian/TEMPB/

/work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations.py compile -P Predictions_Dir/done -S all.GT.subs.longform.txt 2> Compile.log > all.GT.BAD_calls.txt

source deactivate
