#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 11:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 1   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gatk4"   # job name
#SBATCH --mem=150G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools10.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel
module load vcftools/0.1.14-xp36ajy
#module load bcftools/1.9-womp5fh
module load r/3.6.3-py3-sxv6dw3
#module load python/3.8.0-k7w5uj4
module load r-dplyr/0.7.5-py2-r3.5-4jcgekn
module load r-stringr/1.2.0-py2-r3.5-uoogz7r 
module load r-reshape2/1.4.2-py2-ywy5u5c


########
#parallel vcftools --vcf {}.Nocall.vcf --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 --recode --out {}.biallele ::: {01..26}

#for i in {01..26}; do vcftools --vcf $i.biallele.recode.vcf --indv "C2-4 " --indv "G2-21 " --indv "G3-83" --freq2 --out $i.Australians; done  

#for i in *Australians.frq; do grep -v "CHROM" $i | awk 'OFS="\t" {if($4 > 4 && $6 == 1) {print $1,$2}}' > ${i%Australians.frq}keep.pos; grep -v "CHROM" $i | awk 'OFS="\t" {if($4 > 4 && $6 == 0) {print $1,$2}}' >> ${i%Australians.frq}keep.pos; done

#parallel vcftools --vcf {}.biallele.recode.vcf --positions {}.keep.pos --recode --recode-INFO-all --out {}.Australian ::: {01..26}

#for i in {01..26}; do vcftools --vcf $i.Australian.recode.vcf --remove-indv "C2-4 " --remove-indv "G2-21 " --remove-indv "G3-83" --freq2 --out $i.Australians.not; done

#for i in *Australians.not.frq; do grep -v "CHROM" $i | awk 'OFS="\t" {if($6 < 1 && $6 > 0) {print $1,$2}}' > ${i%Australians.not.frq}keep.AD.pos; done 

#parallel vcftools --vcf {}.Australian.recode.vcf --positions {}.keep.AD.pos --recode --recode-INFO-all --out {}.ancestral ::: {01..26}

#parallel bcftools annotate -x ^FORMAT/GT -o {}.GT.vcf {}.ancestral.recode.vcf ::: {01..26}

#bcftools concat -o all.GT.vcf *GT.vcf

#bcftools reheader --samples header.txt -o all.GT.reheader.vcf all.GT.vcf && mv all.GT.reheader.vcf all.GT.vcf && gzip all.GT.vcf

#mkdir frqs

#ls *lst | parallel vcftools --gzvcf all.GT.vcf.gz --keep {} --freq2 --out frqs/{.}

#Rscript merge_frqs.R

#parallel 'python vcf_all_hets.py {}.vcf > {}.hets' ::: {01..26}

#cat *hets | gzip > all.hets.gz && rm *hets

#tar -cvf - all_gVCFs/ | lz4 - all_BQ20_Australianoutgroup.lz4

Rscript filter_hets.R
