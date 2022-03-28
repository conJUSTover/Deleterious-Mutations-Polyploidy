#!/bin/bash

#Submit this script with: sbatch thefilename
# note: change the memory, threads, wall, etc

#SBATCH -t 24:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores
#SBATCH -J "15 run"   # job name
#SBATCH --mem=340G # how much memory you need; each box has ~350G
#SBATCH --output=slurm-count_15_R.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load parallel
module load r
module load r-data-table

cd ../
ls [012]*nofilter.recode.vcf | parallel -j 26 'grep -v "#" {} | sed "s/\:[^\t]*//g" | gzip > GERP_calculations/{.}.nocolon.vcf.gz'
cd GERP_calculations
ls count*[0-9].new.R | parallel -j 26 'Rscript {}'
