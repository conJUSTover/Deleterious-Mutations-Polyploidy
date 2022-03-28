#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 23:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 26   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gatk4"   # job name
#SBATCH --mem=450G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=vcftools2.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=whatever --exclude=legion-[1-8] #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel
module load vcftools/0.1.14-xp36ajy
module load bcftools/1.9-womp5fh
#module load r/3.6.3-py3-sxv6dw3
#module load python/3.8.0-k7w5uj4
#module load r-dplyr/0.7.5-py2-r3.5-4jcgekn
#module load r-stringr/1.2.0-py2-r3.5-uoogz7r 
#module load r-reshape2/1.4.2-py2-ywy5u5c
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx

########
#parallel 'bcftools reheader --samples header.txt -o {}.vcf ../{}.vcf' ::: {01..26}  
#
#parallel vcftools --vcf {}.vcf --recode-INFO-all --remove-indels --min-alleles 2 --max-alleles 2 --recode --out {}.prelim ::: {01..26}
#echo "step 1"
#parallel gatk VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V {}.prelim.recode.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --output {}.filtered_snp.vcf "2> {}.VariantFiltration_snp.err" ::: {01..26}
#echo "step2" 
#
#module load python/3.8.0-k7w5uj4
#parallel 'python ExcessHets_bypop.py {}.filtered_snp.vcf pops.lsts > {}.filtered_snp_exhet.vcf' ::: {01..26}
#echo "step3" 
#
#module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
#
#parallel gatk IndexFeatureFile -F {}.filtered_snp_exhet.vcf ::: {01..26}
#parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp_exhet.vcf --set-filtered-gt-to-nocall --output {}.biallele.recode.vcf "2> {}.Nocall.err" ::: {01..26}
#echo "step4"
#parallel vcftools --vcf {}.biallele.recode.vcf --keep out.lst --freq2 --out {}.Australians ::: {01..26}
#echo "step5"
#for i in *Australians.frq; do grep -v "CHROM" $i | awk 'OFS="\t" {if($4 > 4 && $6 == 1) {print $1,$2}}' > ${i%Australians.frq}keep.pos; grep -v "CHROM" $i | awk 'OFS="\t" {if($4 > 4 && $6 == 0) {print $1,$2}}' >> ${i%Australians.frq}keep.pos; done
#echo "step6"
#parallel vcftools --vcf {}.biallele.recode.vcf --positions {}.keep.pos --recode --recode-INFO-all --out {}.Australian ::: {01..26}
#echo "step7"
#parallel vcftools --vcf {}.biallele.recode.vcf --remove out.lst --freq2 --out {}.Australians.not ::: {01..26}
#echo "step8"
#for i in *Australians.not.frq; do grep -v "CHROM" $i | awk 'OFS="\t" {if($6 < 1 && $6 > 0) {print $1,$2}}' > ${i%Australians.not.frq}keep.AD.pos; done 
#echo "step9"
#parallel vcftools --vcf {}.Australian.recode.vcf --positions {}.keep.AD.pos --recode --recode-INFO-all --out {}.ancestral ::: {01..26}
#echo "step10"
#parallel bcftools annotate -x ^FORMAT/GT -o {}.GT.vcf {}.ancestral.recode.vcf ::: {01..26}
#echo "step11"
#for i in *GT.vcf; do cat $i | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > ${i%.vcf}_sorted.vcf; done
#echo "step 11.5"
#bcftools concat -o all.GT.vcf *GT_sorted.vcf && gzip all.GT.vcf
#echo "step12"
#mkdir frqs

#ls *lst | parallel vcftools --gzvcf all.GT.vcf.gz --keep {} --freq2 --out frqs/{.}
#echo "step13"

module load r/3.6.3-py3-sxv6dw3
module load r-dplyr/0.7.5-py2-r3.5-4jcgekn
module load r-stringr/1.2.0-py2-r3.5-uoogz7r 
module load r-reshape2/1.4.2-py2-ywy5u5c

Rscript merge_frqs.R
echo "step14"
#parallel 'python vcf_all_hets.py {}.vcf > {}.hets' ::: {01..26}

#cat *hets | gzip > all.hets.gz && rm *hets

#tar -cvf - all_gVCFs/ | lz4 - all_BQ20_Australianoutgroup.lz4

#Rscript filter_hets.R
