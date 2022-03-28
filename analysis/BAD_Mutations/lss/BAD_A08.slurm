#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 6-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 32   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "BAD A08"   # job name
#SBATCH --mem=200G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-redo_BAD_A08.out
#SBATCH --mail-type=FAIL
#SBATCH --constraint=AVX2
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram,gpu

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3/4.3.30-qdauveb
module load parallel 
#mkdir Output_Dir
#mkdir Predictions_Dir
#mkdir scratch
#mkdir scratch/BAD_Mutations_Data
#mkdir scratch/BAD_Mutations_Deps

#./BAD_Mutations.py -v DEBUG setup -b scratch/BAD_Mutations_Data -d scratch/BAD_Mutations_Deps -t 'hordeum_vulgare' -e 0.05 -c BAD_Mutations_Config.txt 2> Setup.log

#./BAD_Mutations.py setup --list-species

#./BAD_Mutations.py -v DEBUG fetch -c BAD_Mutations_Config.txt -b scratch/BAD_Mutations_Data --convert-only

source activate bad_mutations

export TMPDIR=/work/LAS/jfw-lab/jconover/deleterious_mutations/redo_Australian/TEMPA/

#/work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations.py align -c /work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations_Config.txt -f fastas/{.}.fasta -o Output_Dir 2> Align_logs/{.}_Alignment.log

ls subs/GohirA08G* | sed 's|subs/||g' | parallel -j 16 '/work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations.py align -c /work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations_Config.txt -f fastas/{.}.fasta -o Output_Dir 2> Align_logs/{.}_Alignment.log; /work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations.py predict -c /work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations_Config.txt -f fastas/{.}.fasta -a Output_Dir/{.}_MSA.fasta -r Output_Dir/{.}.tree -s subs/{.}.subs -o Predictions_Dir 2> Predict_logs/{.}_Predictions.log'

source deactivate

#/work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations.py align -c LRTPredict_Config.txt -f fastas/{.}.fasta -o Output_Dir 2> Align_logs/{.}_Alignment.log; head -1 Output_Dir/{.}.tree > Output_Dir/{.}.tre && mv Output_Dir/{.}.tre Output_Dir/{.}.tree; /work/LAS/jfw-lab/jconover/BAD_Mutations_dev/BAD_Mutations/BAD_Mutations.py predict -c LRTPredict_Config.txt -f fastas/{.}.fasta -a Output_Dir/{.}_MSA.fasta -r Output_Dir/{.}.tree -s subs/{.}.subs -o Predictions_Dir 2> Predict_logs/{.}_Predictions.log'


#./BAD_Mutations.py -v DEBUG \
#                     compile \
#                     -P Predictions_Dir 2> Compile.log
