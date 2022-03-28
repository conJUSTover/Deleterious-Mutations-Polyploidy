#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "SRR6311686"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-redo_SRR6311686.dump.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=legion -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#module load sentieon-genomics/201808.01-opfuvzr
module load parallel
module load sra-toolkit


#grep -v "don" samples.code.lst | grep -v "^8" | cut -f 3 | parallel 'fastq-dump --split-files --gzip {}'
#prefetch -O . --max-size 200GB SRR6311686 
fasterq-dump --split-files /work/LAS/jfw-lab/del_mut_gvcf/SRR6311650.sra
#fastq-dump --split-files SRR6311686


#time ( sentieon driver -t 10 -r Ghirsutum_458_v1.0.26chr.fa --interval CHROM --algo GVCFtyper CHROM.vcf *.reID.gz) &> CHROM.vcf.timelog
