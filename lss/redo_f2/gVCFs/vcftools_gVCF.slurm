#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 12:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 10   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gVCF"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-vcftools.sentieon2.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch,bigram #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load bcftools/1.9-womp5fh
module load vcftools/0.1.14-xp36ajy
module load parallel
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
module load sentieon-genomics/201808.01-opfuvzr

#ls [CG1-9]*gVCF | parallel -j 46 'vcftools --vcf {} --positions Australian/pSONIC.homoeologs.GT.pos --recode --out {.}.pSONIC && gatk IndexFeatureFile -F {.}.pSONIC.recode.vcf'
#time ( sentieon driver -t 10 -r Ghirsutum_458_v1.0.26chr.fa --interval Australian/pSONIC.percentfiltered.375_625.sorted.bed --algo GVCFtyper --emit_mode all pSONIC.homoeoSNPs.combined.allsites.vcf [CG1-9]*.gVCF) &> pSONIC.vcf.timelog
vcftools --vcf pSONIC.homoeoSNPs.combined.allsites.vcf --remove-indels --positions Australian/pSONIC.homoeoSNPs.GT.pos --max-alleles 2 --recode --recode-INFO-all --out pSONIC.combined.SNPs
bcftools reheader --samples Australian/header.txt -o pSONIC.reheader.vcf pSONIC.combined.SNPs.recode.vcf
bcftools annotate -x ^FORMAT/GT -o pSONIC.GT.vcf pSONIC.reheader.vcf
ls Australian/*lst | parallel 'vcftools --vcf pSONIC.GT.vcf --keep {} --freq2 --out {.}.pSONIC.GT'
head -100 pSONIC.GT.vcf | grep "#CHROM" | sed 's/#//g' | cut -f 1,2,4,5 > pSONIC.GT.pos.ref.alt.txt
grep -v "#" pSONIC.GT.vcf | cut -f 1,2,4,5 >> pSONIC.GT.pos.ref.alt.txt
