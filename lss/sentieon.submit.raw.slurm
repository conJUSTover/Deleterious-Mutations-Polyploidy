#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 72:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 10   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "wget SRA"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-remap3.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#module load sentieon-genomics/201808.01-opfuvzr
module load parallel
#module load sra-toolkit


#grep -v "don" samples.code.lst | grep -v "^8" | cut -f 3 | parallel 'fastq-dump --split-files --gzip {}'
#cat samples.remap3.lst | parallel 'prefetch -O . --max-size 200GB {}; fasterq-dump --split-files {}.sra'
#cat samples.remap3.lst | parallel 'fastq-dump --split-files /work/LAS/jfw-lab/del_mut_gvcf/{}.sra'
ls *fastq | parallel 'gzip $i'

#time ( sentieon driver -t 10 -r Ghirsutum_458_v1.0.26chr.fa --interval CHROM --algo GVCFtyper CHROM.vcf *.reID.gz) &> CHROM.vcf.timelog
