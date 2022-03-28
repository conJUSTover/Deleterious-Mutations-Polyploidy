#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 2-00:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 21   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "easysfs"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=easysfs.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch 

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE> 
module load py-pip/9.0.1-py3-dpds55c
module load vcftools/0.1.14-xp36ajy
module load parallel

ls *sites | parallel vcftools --vcf ../all.GT.polarized.vcf --positions {} --recode --recode-INFO-all --out {.} 
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 1.pop.lst --unfolded -a -o {.}.AD1.sfs --proj 20
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 3.pop.lst --unfolded -a -o {.}.AD3.sfs --proj 14
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 4.pop.lst --unfolded -a -o {.}.AD4.sfs --proj 10
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 5.pop.lst --unfolded -a -o {.}.AD5.sfs --proj 8
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 6.pop.lst --unfolded -a -o {.}.AD6.sfs --proj 12
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 7.pop.lst --unfolded -a -o {.}.AD7.sfs --proj 4
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 8.pop.lst --unfolded -a -o {.}.D5.sfs --proj 8
ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 9.pop.lst --unfolded -a -o {.}.A1.sfs --proj 10


#python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i all.GT.vcf.gz -p sfs.lst --unfolded -a -o all.easysfs.unfold --proj=20,14,10,8,12,4,8,10,6

# 10 1.pop.lst
#  7 3.pop.lst
#  5 4.pop.lst
#  4 5.pop.lst
#  6 6.pop.lst
#  2 7.pop.lst
#  4 8.pop.lst
#  5 9.pop.lst
