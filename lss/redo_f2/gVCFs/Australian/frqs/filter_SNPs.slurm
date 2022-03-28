#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gVCF"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-NO_GC_sfs.polar4.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#module load bcftools/1.9-womp5fh
module load vcftools/0.1.14-xp36ajy
#module load parallel
#module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
#module load sentieon-genomics/201808.01-opfuvzr

#vcftools --gzvcf ../all.GT.vcf.gz --positions NO_GC.homoeoSNPs.pos --recode --out pSONIC.NO_GC
#vcftools --vcf all.GT.polarized.vcf --positions NO_GC.homoeoSNPs.pos --recode --out pSONIC.NO_GC.polarized 

for i in [AD]*homoeoSNPs.pos; do vcftools --gzvcf ../all.GT.vcf.gz --positions $i --recode --out pSONIC.${i%.homoeoSNPs.pos}; done 
sbatch Rscript_add_to_vcf.slurm

#module load py-pip/9.0.1-py3-dpds55c
#module load vcftools/0.1.14-xp36ajy
#module load parallel

#python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i pSONIC.NO_GC.recode.vcf -p AD3.pop.lst -a -o pSONIC.NO_GC.AD3.sfs --proj 14
#python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i pSONIC.NO_GC.polarized.recode.vcf -p AD3.pop.split.lst --unfolded -a -o pSONIC.NO_GC.AD3.split.sfs --proj 8,6


#ls *sites | parallel vcftools --vcf ../all.GT.polarized.vcf --positions {} --recode --recode-INFO-all --out {.}
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 1.pop.lst --unfolded -a -o {.}.AD1.sfs --proj 20
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 3.pop.lst --unfolded -a -o {.}.AD3.sfs --proj 14
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 4.pop.lst --unfolded -a -o {.}.AD4.sfs --proj 10
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 5.pop.lst --unfolded -a -o {.}.AD5.sfs --proj 8
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 6.pop.lst --unfolded -a -o {.}.AD6.sfs --proj 12
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 7.pop.lst --unfolded -a -o {.}.AD7.sfs --proj 4
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 8.pop.lst --unfolded -a -o {.}.D5.sfs --proj 8
#ls *vcf | parallel python3 /work/LAS/jfw-lab/jconover/Selfing_hetero/Mapping_Page_et_al/easySFS/easySFS.py -i {} -p 9.pop.lst --unfolded -a -o {.}.A1.sfs --proj 10


