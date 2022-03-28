#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 36:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 13   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "gatk4"   # job name
#SBATCH --mem=450G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=GATK4-exhet30.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jconover@iastate.edu   # email address
#SBATCH --partition=speedy,biocrunch #-- use sinfo to get queue names

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load parallel
module load gatk/4.0.4.0-py2-r3.5-2rrbjyx
module load vcftools/0.1.14-xp36ajy
#module load picard/2.17.0-ft5qztz

#picard CreateSequenceDictionary R=Ghirsutum_458_v1.0.26chr.fa O=Ghirsutum_458_v1.0.26chr.dict

#ExHet=30.0

#parallel -j 13 gatk VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V no_excessHet/{}.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --filter-expression "ExcessHet\>$ExHet" --filter-name "ExHet$ExHet" --output {}.filtered_snp.$ExHet.vcf "2> {}.VariantFiltration_snp.$ExHet.err" ::: {01..26}


#parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.$ExHet.vcf --exclude-filtered --output {}.PASS_snp.$ExhHet.vcf "2> {}.PASS_snp.$ExHet.err" ::: {01..26}

#parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.$ExHet.vcf --set-filtered-gt-to-nocall --output {}.Nocall.$ExHet.vcf "2> {}.Nocall.$ExHet.err" ::: {01..26}

#parallel vcftools --vcf {}.Nocall.$ExHet.vcf --recode-INFO-all --remove-indels --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode --out {}.biallele.$ExHet ::: {01..26}

#mv [012]*.biallele.recode.vcf ../gerp_files/
#cd ../gerp_files
#sbatch Rscript.slurm

ExHet=40.0

parallel -j 13 gatk VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V no_excessHet/{}.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --filter-expression "ExcessHet\>$ExHet" --filter-name "ExHet$ExHet" --output {}.filtered_snp.$ExHet.vcf "2> {}.VariantFiltration_snp.$ExHet.err" ::: {01..26}


#parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.$ExHet.vcf --exclude-filtered --output {}.PASS_snp.$ExhHet.vcf "2> {}.PASS_snp.$ExHet.err" ::: {01..26}

parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.$ExHet.vcf --set-filtered-gt-to-nocall --output {}.Nocall.$ExHet.vcf "2> {}.Nocall.$ExHet.err" ::: {01..26}

parallel vcftools --vcf {}.Nocall.$ExHet.vcf --recode-INFO-all --remove-indels --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode --out {}.biallele.$ExHet ::: {01..26}

ExHet=60.0

parallel -j 13 gatk VariantFiltration -R Ghirsutum_458_v1.0.26chr.fa -V no_excessHet/{}.vcf --genotype-filter-expression "QD\<2.0" --genotype-filter-name "QD2" --genotype-filter-expression "FS\>60.0" --genotype-filter-name "FS60" --genotype-filter-expression "MQ\<40.0" --genotype-filter-name "MW40" --genotype-filter-expression "SOR\>4.0" --genotype-filter-name "SOR4" --genotype-filter-expression "MQRankSum\<-12.5" --genotype-filter-name "MQRankSum-12.5" --genotype-filter-expression "ReadPosRankSum\<-8.0" --genotype-filter-name "ReadPosRankSum" --genotype-filter-expression "DP\>60" --genotype-filter-name "DP60" --filter-expression "ExcessHet\>$ExHet" --filter-name "ExHet$ExHet" --output {}.filtered_snp.$ExHet.vcf "2> {}.VariantFiltration_snp.$ExHet.err" ::: {01..26}


#parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.$ExHet.vcf --exclude-filtered --output {}.PASS_snp.$ExhHet.vcf "2> {}.PASS_snp.$ExHet.err" ::: {01..26}

parallel gatk SelectVariants -R Ghirsutum_458_v1.0.26chr.fa --variant {}.filtered_snp.$ExHet.vcf --set-filtered-gt-to-nocall --output {}.Nocall.$ExHet.vcf "2> {}.Nocall.$ExHet.err" ::: {01..26}

parallel vcftools --vcf {}.Nocall.$ExHet.vcf --recode-INFO-all --remove-indels --remove-filtered-all --min-alleles 2 --max-alleles 2 --recode --out {}.biallele.$ExHet ::: {01..26}

