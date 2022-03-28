#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "python fasta"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-python_fasta.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=legion -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load python
module load py-biopython/1.70-py3-wos466g

python /work/LAS/jfw-lab/jconover/scripts/acc_to_fasta.py Ghirsutum_458_v1.0.26chr.fa Achr.lst > A_genome/Ghirsutum_458_v1.0.26chr.fa
python /work/LAS/jfw-lab/jconover/scripts/acc_to_fasta.py Ghirsutum_458_v1.0.26chr.fa Dchr.lst > D_genome/Ghirsutum_458_v1.0.26chr.fa
