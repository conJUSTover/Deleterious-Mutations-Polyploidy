#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 2-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 78   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "BAD A#CHROM"   # job name
#SBATCH --mem=400G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-redo_BAD_A#CHROM.out
#SBATCH --mail-type=FAIL
#SBATCH --constraint=AVX2
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load gcc/7.3.0-xegsmw4
module load hyphy/2.3.14-iggx7u3
module load blast-plus/2.7.1-py2-ybam4tg
module load python/2.7.14-h73plf5
module load py-biopython/1.65-py2-fgq2vqv
#module load pasta/1.8.3-py2-alcbobp
module load py-argparse/1.4.0-py2-bfxyi6x
module load pasta/1.8.3-py2-alcbobp
module load parallel

export TMPDIR=/work/LAS/jfw-lab/jconover/deleterious_mutations/BAD_runs/redo_Australian/TEMP/

ls subs/GohirA#CHROMG* | sed 's|subs/||g' | parallel -j 26 '/work/LAS/jfw-lab/jconover/BAD_Mutations/BAD_Mutations.py align -c LRTPredict_Config.txt -f fastas/{.}.fasta -o Output_Dir 2> Align_logs/{.}_Alignment.log; head -1 Output_Dir/{.}.tree > Output_Dir/{.}.tre && mv Output_Dir/{.}.tre Output_Dir/{.}.tree; /work/LAS/jfw-lab/jconover/BAD_Mutations/BAD_Mutations.py predict -c LRTPredict_Config.txt -f fastas/{.}.fasta -a Output_Dir/{.}_MSA.fasta -r Output_Dir/{.}.tree -s subs/{.}.subs -o Predictions_Dir 2> Predict_logs/{.}_Predictions.log'


