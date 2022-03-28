#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 1-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 5   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "Dpython"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-acc_pythonDt.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load gcc/7.3.0-xegsmw4
module load py-biopython/1.73-py3-7gocqag
module load parallel

python /work/LAS/jfw-lab/jconover/scripts/acc_to_fasta.py ../map/Ghirsutum_458_v1.0.26chr.fa chr.lst > Ghirsutum_458_v1.0.26chr.fa

module load bwa/0.7.17-zhcbtza

time ( bwa index Ghirsutum_458_v1.0.26chr.fa )

#module use /opt/rit/local-test/modules
module load sentieon-genomics/201808.01-opfuvzr
module load samtools/1.9-k6deoga
#module load trimmomatic/0.36-lkktrba

samtools faidx Ghirsutum_458_v1.0.26chr.fa



